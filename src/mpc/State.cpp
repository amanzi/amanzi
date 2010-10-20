#include "State.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Mesh_maps_simple.hh"

State::State( int number_of_components_,
	      Teuchos::RCP<Mesh_maps_simple> mesh_maps_):
  number_of_components(number_of_components_),
  mesh_maps(mesh_maps_)
{
  // create the Eptera_Vector objects

  water_density =    Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  pressure =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  darcy_flux =       Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  porosity =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  water_saturation = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) ); 
  total_component_concentration 
     = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), number_of_components ) );  

};


void State::set_time ( double new_time ) {

  if ( status == UPDATING ) {
    
    time = new_time;

  } else {
    
    // throw an error

  }

}
