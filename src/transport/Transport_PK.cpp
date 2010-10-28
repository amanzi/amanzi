#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"

#include "simple_mesh/Mesh_maps_base.hh"
#include "Transport_PK.hpp"
#include "flow/cell_geometry.hpp"

using namespace std;
using namespace Teuchos;
using namespace cell_geometry;


/*  Constructor for initializing the transport PK.    */
/*  Its call is made usually at time T=0 by the MPC.  */
Transport_PK::Transport_PK ( ParameterList &parameter_list_MPC,
			     RCP<Transport_State> TS_MPC )

{ 
  parameter_list = parameter_list_MPC;

  /* make a copy of the transport state object */
  TS = rcp( new Transport_State() );
  TS->copy_constant_state( *TS_MPC );

  /* allocate memory for internal (next) transport state */
  TS_next = rcp( new Transport_State () );
  TS_next->create_internal_state( *TS_MPC );

  /* set null/zero values to all internal parameters */
  dT = 0.01;
  status = 0;


  /* future geometry package */
  Epetra_Map face_map = TS->get_mesh_maps()->face_map(false);

  face_area  = new Epetra_Vector(face_map);
  face_to_cell_upwind   = new Epetra_IntVector(face_map); 
  face_to_cell_downwind = new Epetra_IntVector(face_map); 

  geometry_package();


  /* read the CFL number from the parameter list with a default of 1.0 */
  cfl = parameter_list.get<double>("CFL", 1.0);
};



/* destructor */
Transport_PK::~Transport_PK()
{
  delete face_area;
  delete face_to_cell_upwind;
  delete face_to_cell_downwind;
}




/*  This is part of the future geometry package  */
void Transport_PK::geometry_package()
{
  RCP<Mesh_maps_base> mesh = TS->get_mesh_maps();


  /* loop over faces and calculate areas */
  unsigned int i, f, c;
  double x[4][3], v1[3], v2[3];
  Epetra_Map face_map = mesh->face_map(false);

  for ( f=face_map.MinLID(); f<face_map.MaxLID(); f++ ) {
     mesh->face_to_coordinates( f, (double*) x, (double*) x+12 );

     (*face_area)[f] = quad_face_area(x[0], x[1], x[2], x[3]);
  }


  /* loop over cells */
  vector<unsigned int> c2f(6);
  Epetra_Map cell_map = mesh->cell_map(false);

  /* clean the reserce map f->c */
  for ( f=face_map.MinLID(); f<face_map.MaxLID(); f++ ) {
     (*face_to_cell_upwind)[f]   = -1;
     (*face_to_cell_downwind)[f] = -1;
  }

  /* populate upwind and downwind cells LID */
  for ( c=cell_map.MinLID(); c<cell_map.MaxLID(); c++) {
     mesh->cell_to_faces( c, c2f.begin(), c2f.end() );

     for ( i=0; i<6; i++ ) {
        f = face_map.LID(c2f[i]);
        if ( (*face_to_cell_upwind)[f] == -1 ) { (*face_to_cell_upwind)[f] = c; }
        else                                   { (*face_to_cell_downwind)[f] = c; }
     }
  }

  /* select the upwinding cell */
  int c1, c2;
  double *darcy_flux;

  TS->get_darcy_flux()->ExtractView( &darcy_flux );

  for ( f=face_map.MinLID(); f<face_map.MaxLID(); f++ ) {
     c1 = (*face_to_cell_upwind)[f]; 
     c2 = (*face_to_cell_downwind)[f]; 

     if ( darcy_flux[f] >= 0 ) {
        (*face_to_cell_upwind)[f]   = max(c1, c2); 
        (*face_to_cell_downwind)[f] = min(c1, c2); 
     } else {
        (*face_to_cell_upwind)[f]   = min(c1, c2); 
        (*face_to_cell_downwind)[f] = max(c1, c2); 
     }
  }
}



/* MPC will call this function to advance the transport state  */
void Transport_PK::advance()
{
  /* this should be moved to MPC */
  //calculate_transport_dT();

  /* Step 1: Loop over internal faces: update concentrations */
  int i, c1, c2;
  unsigned int f;
  double u, tcc_mass_flux;

  /* access raw data */
  double **tcc_data, **tcc_next_data;

  RCP<Epetra_MultiVector>  tcc      = TS->get_total_component_concentration ();
  RCP<Epetra_MultiVector>  tcc_next = TS_next->get_total_component_concentration ();

  tcc->ExtractView( &tcc_data );
  tcc_next->ExtractView( &tcc_next_data );

  double *darcy_flux;
  TS->get_darcy_flux()->ExtractView( &darcy_flux );

  RCP<Mesh_maps_base> mesh = TS->get_mesh_maps();
  Epetra_Map face_map = mesh->face_map(false);

  /* advance each component */ 
  int num_components = tcc->NumVectors();

  for ( f=face_map.MinLID(); f<face_map.MaxLID(); f++ ) {
     c1 = (*face_to_cell_upwind)[f]; 
     c2 = (*face_to_cell_downwind)[f]; 

     if ( c1 >=0 && c2 >= 0 ) {
        u = darcy_flux[f];

        for ( i=0; i<num_components; i++ ) {
            tcc_mass_flux = cfl * dT * u * tcc_next_data[i][c2];

            tcc_next_data[i][c1] = tcc_data[i][c1] + tcc_mass_flux;
            tcc_next_data[i][c2] = tcc_data[i][c2] - tcc_mass_flux;
        }
     }
  }


  /* Step 2: Create an interface map */  

  /* Step 3: Parallel communication */

  /* Step 4: Loop over boundary faces */
};



/*  MPC will call this function to indicate to the transport PK  */
/*  that it can commit the advanced state it has created.        */
/*  This  call indicates that the MPC has accepted the new state */
void Transport_PK::commit_state ( RCP<Transport_State> TS )
{
  /* nothing is done her since a pointer to the state is kept */ 
};



/*  DEBUGing routines  */
double Transport_PK::get_face_area( int f ) {
  double *internal_data;

  face_area->ExtractView( &internal_data );

  return internal_data[f];
}
