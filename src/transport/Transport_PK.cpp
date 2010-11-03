#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"

#include "simple_mesh/Mesh_maps_base.hh"
#include "Transport_PK.hpp"
#include "flow/cell_geometry.hpp"

using namespace std;
using namespace Teuchos;
using namespace cell_geometry;



/* ************************************************************* */
/*  Constructor for initializing the transport PK.               */
/*  Its call is made usually at time T=0 by the MPC.             */
/* ************************************************************* */
Transport_PK::Transport_PK( ParameterList &parameter_list_MPC,
			    RCP<Transport_State> TS_MPC )

{ 
  parameter_list = parameter_list_MPC;

  number_components = TS_MPC->get_total_component_concentration()->NumVectors();

  /* make a copy of the transport state object */
  TS = rcp( new Transport_State() );
  TS->copy_constant_state( *TS_MPC );

  /* allocate memory for internal (next) transport state */
  TS_next = rcp( new Transport_State () );
  TS_next->create_internal_state( *TS_MPC );

  /* set null/zero values to all internal parameters */
  dT = 0.0;
  status = TRANSPORT_NULL;


  /* frequently used data */
  Epetra_Map cell_map = TS->get_mesh_maps()->cell_map(false);
  Epetra_Map face_map = TS->get_mesh_maps()->face_map(false);

  cmin = cell_map.MinLID();
  cmax = cell_map.MaxLID();
  number_cells = cmax - cmin + 1;

  fmin = face_map.MinLID();
  fmax = face_map.MaxLID(); 
  number_faces = fmax - fmin + 1;


  /* process parameter list */
  process_parameter_list();

  /* future geometry package */
  cell_volume.resize( cmax + 1 );
  face_area.resize( fmax + 1 );

  upwind_cell.resize( fmax + 1 );
  downwind_cell.resize( fmax + 1 );

  geometry_package();

  /* other preliminaries for Demo I */
  darcy_flux.resize( fmax + 1 );
};




/* ************************************************************* */
/* process parameter list: needs to be called only once on each  */ 
/* processor                                                     */
/* ************************************************************* */
void Transport_PK::process_parameter_list()
{
  RCP<Mesh_maps_base> mesh = TS->get_mesh_maps();

  /* global transport parameters */
  cfl = parameter_list.get<double>( "CFL", 1.0 );
 
  // the number of components is given by the state and not read
  // from the parameter list - M.B.
  // number_components = parameter_list.get<int>( "number of components" );

  cout << "Transport PK: CFL = " << cfl << endl;
  cout << "              Total number of components = " << number_components << endl;


  /* read number of boundary consitions */ 
  ParameterList  BC_list;
  int  i, k, nBCs;

  BC_list = parameter_list.get<ParameterList>("Transport BCs");
  nBCs = BC_list.get<int>("number of BCs");


  /* create list of boundary data */
  bcs.resize( nBCs );
  for( i=0; i<nBCs; i++ ) bcs[i] = Transport_BCs( 0, number_components );
  
  for( i=0; i<nBCs; i++ ) {
     char bc_char_name[10];
    
     sprintf(bc_char_name, "BC %d", i);
     string bc_name(bc_char_name);

     if ( ! BC_list.isSublist( bc_name ) ) throw exception();

     ParameterList bc_ss = BC_list.sublist( bc_name );

     unsigned int  ssid, ntcc;
     string  type;
     double  value;

     ssid = bc_ss.get<int>("Side set ID");
     ntcc = bc_ss.get<int>("number of components");
     type = bc_ss.get<string>("Type");

     /* check all existing components */
     for( k=0; k<number_components; k++ ) {
        char tcc_char_name[10];

        sprintf( tcc_char_name, "Component %d", k );
        string tcc_name( tcc_char_name );

        if ( bc_ss.isParameter( tcc_name ) ) { value = bc_ss.get<double>( tcc_name ); }
        else                                 { value = 0.0; }

        bcs[i].values[k] = value;
     }

     bcs[i].side_set_id = ssid;
     if ( !mesh->valid_set_id( ssid, Mesh_data::FACE ) ) throw exception();

     /* populate list of n boundary faces: it could be empty */
     int  n;

     n = mesh->get_set_size( ssid, Mesh_data::FACE, OWNED );
     bcs[i].faces.resize( n );

     mesh->get_set( ssid, Mesh_data::FACE, OWNED, bcs[i].faces.begin(), bcs[i].faces.end() );
  }
}




/* ************************************************************* */
/* Estimation of the time step based on T.Barth (Lecture Notes   */
/* presented at VKI Lecture Series 1994-05, Theorem 4.2.2.       */
/* Routine must be called every time we update a flow field      */
/* ************************************************************* */
double Transport_PK::calculate_transport_dT()
{
  if ( status == TRANSPORT_NULL ) {
     extract_darcy_flux();
     identify_upwind_cells();
  }

  RCP<Mesh_maps_base> mesh = TS->get_mesh_maps();

  vector<double>  total_influx(number_cells, 0.0);

  /* loop over faces and accumulate upwinding fluxes */
  int  i, f, c, c1;
  double  area, u;
  Epetra_Map face_map = mesh->face_map(false);

  for( f=fmin; f<=fmax; f++ ) {
     c = downwind_cell[f];

     area = face_area[f];
     if( c >= 0 ) total_influx[c] += area * fabs( darcy_flux[f] ); 
  }


  /* loop over cells and calculate minimal dT */
  double  influx, dT_cell; 

  RCP<const Epetra_Vector>  ws  = TS->get_water_saturation();
  RCP<const Epetra_Vector>  phi = TS->get_porosity();

  dT = dT_cell = TRANSPORT_LARGE_TIME_STEP;
  for( c=cmin; c<=cmax; c++ ) {
     influx = total_influx[c];
     if( influx ) dT_cell = cell_volume[c] * (*phi)[c] * (*ws)[c] / influx;

     dT = min( dT, dT_cell);
  }


  /* parallel garther and scatter of dT */ 
  //cout << "Transport time step dT = " << dT << endl;
}




/* ************************************************************* */
/* MPC will call this function to advance the transport state    */
/* ************************************************************* */
void Transport_PK::advance()
{
  /* this should be moved to MPC */
  calculate_transport_dT();

  status = TRANSPORT_STATE_BEGIN;


  /* Step 1: Loop over internal faces: update concentrations */
  int i, c1, c2;
  unsigned int f;
  double u, phi_ws1, phi_ws2, tcc_flux;

  /* access raw data */
  RCP<Epetra_MultiVector>   tcc      = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector>   tcc_next = TS_next->get_total_component_concentration();

  RCP<const Epetra_Vector>  ws  = TS->get_water_saturation();
  RCP<const Epetra_Vector>  phi = TS->get_porosity();

  /* copy previous state */
  *tcc_next = *tcc; 


  /* advance each component */ 
  int num_components = tcc->NumVectors();

  for( f=fmin; f<=fmax; f++ ) {
     c1 = upwind_cell[f]; 
     c2 = downwind_cell[f]; 

     if ( c1 >=0 && c2 >= 0 ) {
        u = fabs(darcy_flux[f]);

        phi_ws1 = (*phi)[c1] * (*ws)[c1]; 
        phi_ws2 = (*phi)[c2] * (*ws)[c2]; 

        for( i=0; i<num_components; i++ ) {
           tcc_flux = cfl * dT * u * (*tcc)[i][c1];

           (*tcc_next)[i][c1] -= tcc_flux / phi_ws1;
           (*tcc_next)[i][c2] += tcc_flux / phi_ws2;
        }
     }
  }


  /* Step 2: Create an interface map */  

  /* Step 3: Parallel communication */

  /* Step 4: Loop over boundary sets */
  int  k, n;
  for( n=0; n<bcs.size(); n++ ) {
     for( k=0; k<bcs[n].faces.size(); k++ ) {
        f = bcs[n].faces[k];

        c1 = upwind_cell[f];
        c2 = downwind_cell[f]; 

        u = fabs(darcy_flux[f]);

        if ( c1 >= 0 ) {
           phi_ws1  = (*phi)[c1] * (*ws)[c1]; 

           for( i=0; i<num_components; i++ ) {
              tcc_flux = cfl * dT * u * (*tcc)[i][c1];
              (*tcc_next)[i][c1] -= tcc_flux / phi_ws1;
           }
        }

        if ( c2 >= 0 ) {
           phi_ws2 = (*phi)[c2] * (*ws)[c2]; 

           for( i=0; i<num_components; i++ ) {
              tcc_flux = cfl * dT * u * bcs[n].values[i];
              (*tcc_next)[i][c2] += tcc_flux / phi_ws2;
           }
        }
     }
  }


  /* Step 5: Paralell communication */

  status = TRANSPORT_STATE_COMPLETE;
};




/* ************************************************************* */
/*  MPC will call this function to indicate to the transport PK  */
/*  that it can commit the advanced state it has created.        */
/*  This  call indicates that the MPC has accepted the new state */
/* ************************************************************* */
void Transport_PK::commit_state( RCP<Transport_State> TS )
{
  /* nothing is done her since a pointer to the state is kept */ 
};




/* ************************************************************* */
/* calculate Darcy velocity from mass flux:                      */
/* serial implementaition                                        */
/* ************************************************************* */
void Transport_PK::extract_darcy_flux()
{
  int  f, c1, c2;
  double  density;

  RCP<Mesh_maps_base>  mesh = TS->get_mesh_maps();
  RCP<const Epetra_Vector>  rho = TS->get_water_density();

  double *u;
  TS->get_darcy_flux()->ExtractView( &u );

  for( f=fmin; f<=fmax; f++ ) {
     c1 = upwind_cell[f]; 
     c2 = downwind_cell[f]; 

     if ( c1 >= cmin && c2 >= cmin ) { density = ((*rho)[c1] + (*rho)[c2]) / 2; }
     else if ( c1 >= cmin )          { density = (*rho)[c1]; }
     else if ( c2 >= cmin )          { density = (*rho)[c2]; }

     darcy_flux[f] = u[f] / density;
  }
}




/* ************************************************************* */
/* Identify flux direction based on direction of the face normal */
/* andsign of the  Darcy velocity.                               */
/* ************************************************************* */
void Transport_PK::identify_upwind_cells()
{
  RCP<Mesh_maps_base> mesh = TS->get_mesh_maps();

  /* negative value is indicator of a boundary  */
  int  f;
  for( f=fmin; f<=fmax; f++ ) {
     upwind_cell[f]   = -1;
     downwind_cell[f] = -1;
  }

  /* populate upwind and downwind cells ids */
  int  i, c;

  vector<unsigned int>  c2f(6);
  vector<int>           dirs(6);

  Epetra_Map  face_map = mesh->face_map(false);

  double *darcy_flux;
  TS->get_darcy_flux()->ExtractView( &darcy_flux );

  for( c=cmin; c<=cmax; c++) {
     mesh->cell_to_faces( c, c2f.begin(), c2f.end() );
     mesh->cell_to_face_dirs( c, dirs.begin(), dirs.end() );

     for ( i=0; i<6; i++ ) {
        f = face_map.LID(c2f[i]);
        if ( darcy_flux[f] * dirs[i] >= 0 ) { upwind_cell[f] = c; }
        else                                { downwind_cell[f] = c; }
     }
  }
}




/* ************************************************************* */
/*  This is part of the future geometry package                  */
/* ************************************************************* */
void Transport_PK::geometry_package()
{
  RCP<Mesh_maps_base> mesh = TS->get_mesh_maps();

  /* loop over faces and calculate areas */
  unsigned int f;
  double x[8][3], v1[3], v2[3];

  for( f=fmin; f<=fmax; f++ ) {
     mesh->face_to_coordinates( f, (double*) x, (double*) x+12 );

     face_area[f] = quad_face_area(x[0], x[1], x[2], x[3]);
  }


  /* loop over cells, then faces, to calculate cell volume */
  int  i, j, c;
  double  center[3], normal[3], volume;

  vector<unsigned int>  c2f(6);
  vector<int>           dirs(6);

  Epetra_Map  face_map = mesh->face_map(false);

  for( c=cmin; c<=cmax; c++ ) {
     TS->get_mesh_maps()->cell_to_coordinates( c, (double*) x, (double*) x+24 );

     for( i=0; i<3; i++ ) { 
        center[i] = 0;
        for( j=0; j<8; j++ ) center[i] += x[j][i];
        center[i] /= 8;
     }

     mesh->cell_to_faces( c, c2f.begin(), c2f.end() );
     mesh->cell_to_face_dirs( c, dirs.begin(), dirs.end() );
     
     /* assume flat faces */
     volume = 0.0;
     for ( j=0; j<6; j++ ) {
        f = face_map.LID(c2f[j]);
        mesh->face_to_coordinates( f, (double*) x, (double*) x+12 );

        quad_face_normal(normal, x[0], x[1], x[2], x[3]);

        for( i=0; i<3; i++ ) v1[i] = x[0][i] - center[i];

        volume += dirs[j] * dot_product(normal, v1, 3);
     }
     cell_volume[c] = volume / 3;
  }
}


