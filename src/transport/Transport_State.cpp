
#include "../mpc/State.hpp"
#include "Transport_State.hpp"


// this method creates a deep copy of the Transport_State
// object that is passed in
 
void Transport_State::copy ( Teuchos::RCP<Transport_State> TS )
{
  mesh_maps = TS->mesh_maps;

  total_component_concentration = Teuchos::rcp( new Epetra_MultiVector( * TS->total_component_concentration ) );
  water_saturation = Teuchos::rcp( new Epetra_Vector( * TS->water_saturation ) );
  darcy_flux = Teuchos::rcp (new Epetra_Vector( * TS->darcy_flux ) );
  porosity = Teuchos::rcp( new Epetra_Vector( * TS->porosity ) );
}


