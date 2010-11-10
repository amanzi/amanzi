#include "Flow_State.hpp"

Flow_State::Flow_State(Teuchos::RCP<State> S): 
  mesh_maps_(S->get_mesh_maps()),
  gravity_(S->get_gravity()),
  fluid_density_(S->get_density()),
  fluid_viscosity_(S->get_viscosity()),
  permeability_Epetra_(S->get_permeability()),
  pressure_(S->get_pressure())
{

};
