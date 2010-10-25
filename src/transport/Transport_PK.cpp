#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"

#include "simple_mesh/Mesh_maps_base.hh"
#include "Transport_PK.hpp"
#include "flow/cell_geometry.hpp"

using namespace std;
using namespace Teuchos;
using namespace cell_geometry;


/* 
   Constructor for initializing the transport PK.
   Its call is made usually at time T=0 by the MPC.
 */
Transport_PK::Transport_PK ( RCP<Transport_State> TS_MPC )
{ 
  /* make a copy of the transport state object */
  TS = TS_MPC;

  /* copy pointers for state variables that will remain unchanged */
  TS_next = rcp( new Transport_State () );

  TS_next->get_porosity () = TS->get_porosity (); 
  TS_next->get_water_saturation () = TS->get_water_saturation (); 
  TS_next->get_darcy_flux () = TS->get_darcy_flux (); 

  /* allocate memory for state variables that will be changed */
  RCP<Epetra_MultiVector>  tcc      = TS->get_total_component_concentration ();
  RCP<Epetra_MultiVector>  tcc_next = TS_next->get_total_component_concentration ();

  tcc_next = rcp( new Epetra_MultiVector( *tcc ) );


  /* set null/zero values to all internal parameters */
  dT = 0.0;
  status = 0;


  /* future geometry package */
  Epetra_Map face_map = TS->get_mesh_maps()->face_map(false);

  face_area  = new Epetra_Vector(face_map);
  face_to_cell_upwind   = new Epetra_IntVector(face_map); 
  face_to_cell_downwind = new Epetra_IntVector(face_map); 

  geometry_package();
};



/* this is part of the future geometry package */
void Transport_PK::geometry_package()
{
  RCP<const Mesh_maps_simple> mesh = TS->get_mesh_maps();


  /* loop over faces and calculate areas */
  unsigned int i, f, c;
  double x[4][3], v1[3], v2[3];
  Epetra_Map face_map = mesh->face_map(false);

  for ( f=face_map.MinLID(); f<face_map.MaxLID(); f++ ) {
     //mesh->face_to_coordinates( f, (double*) x, (double*) x+12 );

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
     //mesh->cell_to_faces( c, c2f.begin(), c2f.end() );

     for ( i=0; i<6; i++ ) {
        f = face_map.LID(c2f[i]);
     }
  }
}



/* MPC will call this function to advance the state  with this 
   particular process kernel */
void Transport_PK::advance_transport_state ()
{
  /* Step 1: Create reverse map: face -> cells  */

  /* Step 2: Loop over internal faces: update concentrations */

  /* Step 3: Create an interface map */  

  /* Step 4: Parallel communication */

  /* Step 5: Loop over interface faces */
};



/* 
  DEBUGing routines
*/
double Transport_PK::get_face_area( int f ) {
  double *internal_data;

  face_area->ExtractView( &internal_data );

  return internal_data[f];
}
