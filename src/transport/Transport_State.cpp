
#include "../mpc/State.hpp"
#include "Transport_State.hpp"


using namespace Teuchos;


/* at the moment the transport state is a copy of the global state */  
Transport_State::Transport_State ( RCP<State> S )
{
  total_component_concentration = S->get_total_component_concentration();
  porosity                      = S->get_porosity();
  darcy_flux                    = S->get_darcy_flux();
  water_saturation              = S->get_water_saturation();
  mesh_maps                     = S->get_mesh_maps();
}



Transport_State::Transport_State ( State S )
{
  total_component_concentration = S.get_total_component_concentration();
  porosity                      = S.get_porosity();
  darcy_flux                    = S.get_darcy_flux();
  water_saturation              = S.get_water_saturation();
  mesh_maps                     = S.get_mesh_maps();
}



/*
  mesh_maps = TS->mesh_maps;

  total_component_concentration = Teuchos::rcp( new Epetra_MultiVector( * TS->total_component_concentration ) );
  water_saturation = Teuchos::rcp( new Epetra_Vector( * TS->water_saturation ) );
  darcy_flux = Teuchos::rcp (new Epetra_Vector( * TS->darcy_flux ) );
  porosity = Teuchos::rcp( new Epetra_Vector( * TS->porosity ) );
*/


