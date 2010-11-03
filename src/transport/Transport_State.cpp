
#include "mpc/State.hpp"
#include "flow/cell_geometry.hpp"
#include "Transport_State.hpp"


using namespace Teuchos;
using namespace cell_geometry;


/* ************************************************************* */
/* at the moment a transport state is a copy of the global state */  
/* ************************************************************* */
Transport_State::Transport_State( State S )
{
  total_component_concentration = S.get_total_component_concentration();
  porosity                      = S.get_porosity();
  darcy_flux                    = S.get_darcy_flux();
  water_saturation              = S.get_water_saturation();
  water_density                 = S.get_water_density();
  mesh_maps                     = S.get_mesh_maps();
}




/* ************************************************************* */
/* trivial (at the moment) copy of a constant transport state    */  
/* ************************************************************* */
void Transport_State::copy_constant_state( Transport_State & S )
{
  total_component_concentration = S.get_total_component_concentration();
  porosity                      = S.get_porosity();
  darcy_flux                    = S.get_darcy_flux();
  water_saturation              = S.get_water_saturation();
  water_density                 = S.get_water_density();
  mesh_maps                     = S.get_mesh_maps();
}




/* ************************************************************* */
/* internal transport state uses internal variable for the total */
/* component concentration                                       */
/* ************************************************************* */
void Transport_State::create_internal_state( Transport_State & S )
{
  porosity         = S.get_porosity(); 
  darcy_flux       = S.get_darcy_flux(); 
  water_saturation = S.get_water_saturation(); 
  water_density    = S.get_water_density();

  RCP<Epetra_MultiVector> tcc = S.get_total_component_concentration();
  total_component_concentration = rcp( new Epetra_MultiVector( *tcc ) );
}




/* ************************************************************* */
/* DEBUG: create constant analytical Darcy velocity field:       */
/* u = (1,2,3)                                                   */
/* ************************************************************* */
void Transport_State::analytic_darcy_flux()
{
  int  i, f;
  double x[4][3], normal[3], length;

  Epetra_Map face_map = mesh_maps->face_map(false);

  for( f=face_map.MinLID(); f<=face_map.MaxLID(); f++ ) { 
     mesh_maps->face_to_coordinates( f, (double*) x, (double*) x+12 );

     quad_face_normal(normal, x[0], x[1], x[2], x[3]);
     length = vector_length( normal, 3 );

     (*darcy_flux)[f] = (normal[0] + 2 * normal[1] + 3 * normal[2]) / length;
  }
}




/* ************************************************************* */
/* DEBUG: create constant analytical concentration C_0 = x       */
/* ************************************************************* */
void Transport_State::analytic_total_component_concentration()
{
  int  i, j, c;
  double x[8][3], center[3];

  Epetra_Map cell_map = mesh_maps->cell_map(false);

  for( c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++ ) { 
     mesh_maps->cell_to_coordinates( c, (double*) x, (double*) x+24);

     for( i=0; i<3; i++ ) { 
        center[i] = 0;
        for( j=0; j<8; j++ ) center[i] += x[j][i];
        center[i] /= 8;
     }

     (*total_component_concentration)[0][c] = center[0] / 100;
  }
}




/* ************************************************************* */
/* DEBUG: create constant analytical porosity                    */
/* ************************************************************* */
void Transport_State::analytic_porosity( double phi )
{
  int  c;
  Epetra_Map cell_map = mesh_maps->cell_map(false);

  for( c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++ ) { 
     (*porosity)[c] = phi;  /* default is 0.2 */
  }
}




/* ************************************************************* */
/* DEBUG: create constant analytical water saturation            */
/* ************************************************************* */
void Transport_State::analytic_water_saturation( double ws )
{
  int  c;
  Epetra_Map cell_map = mesh_maps->cell_map(false);

  for( c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++ ) { 
     (*water_saturation)[c] = ws;  /* default is 1.0 */
  }
}




/* ************************************************************* */
/* DEBUG: create constant analytical water density               */
/* ************************************************************* */
void Transport_State::analytic_water_density( double wd )
{
  int  c;
  Epetra_Map cell_map = mesh_maps->cell_map(false);

  for( c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++ ) { 
     (*water_density)[c] = wd;  /* default is 1000.0 */
  }
}
