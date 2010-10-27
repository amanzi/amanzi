
#include "../mpc/State.hpp"
#include "Transport_State.hpp"


using namespace Teuchos;


/* at the moment the transport state is a copy of the global state */  
Transport_State::Transport_State ( State S )
{
  total_component_concentration = S.get_total_component_concentration();
  porosity                      = S.get_porosity();
  darcy_flux                    = S.get_darcy_flux();
  water_saturation              = S.get_water_saturation();
  mesh_maps                     = S.get_mesh_maps();
}



/* trivial (at the moment) copy of a constant transport state */  
void Transport_State::copy_constant_state ( Transport_State & S )
{
  total_component_concentration = S.get_total_component_concentration();
  porosity                      = S.get_porosity();
  darcy_flux                    = S.get_darcy_flux();
  water_saturation              = S.get_water_saturation();
  mesh_maps                     = S.get_mesh_maps();
}



/* internal transport state uses internal variable for the total component concentration */
void Transport_State::create_internal_state ( Transport_State & S )
{
  porosity         = S.get_porosity(); 
  water_saturation = S.get_water_saturation(); 
  darcy_flux       = S.get_darcy_flux(); 

  RCP<Epetra_MultiVector> tcc = S.get_total_component_concentration();
  total_component_concentration = rcp( new Epetra_MultiVector( *tcc ) );
}




