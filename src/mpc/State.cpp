#include "State.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "MeshWrapper.hpp"

State::State( Teuchos::RCP<DataLayout> data_layout_ , int number_of_components_,
	      Teuchos::RCP<MeshWrapper> mesh_wrapper_):
  data_layout(data_layout_),
  number_of_components(number_of_components_),
  mesh_wrapper(mesh_wrapper_)
{
  // create the Eptera_Vector objects

  water_density =    Teuchos::rcp( new Epetra_Vector( *data_layout->get_element_map() ) );
  pressure =         Teuchos::rcp( new Epetra_Vector( *data_layout->get_element_map() ) );
  darcy_flux =       Teuchos::rcp( new Epetra_Vector( *data_layout->get_element_map() ) );
  porosity =         Teuchos::rcp( new Epetra_Vector( *data_layout->get_element_map()   ) );
  water_saturation = Teuchos::rcp( new Epetra_Vector( *data_layout->get_element_map() ) ); 
  total_component_concentration 
    = Teuchos::rcp( new Epetra_MultiVector( *data_layout->get_element_map(), number_of_components ) );  

};




void State::set_time ( double new_time ) {

  if ( status == UPDATING ) {
    
    time = new_time;

  } else {
    
    // throw an error

  }

}
