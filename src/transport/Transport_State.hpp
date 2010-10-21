#ifndef __Transport_State_hpp__
#define __Transport_State_hpp__

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh_maps_simple.hh"


using namespace Teuchos;

class Transport_State {

public:
  Transport_State (RCP<State> S):
     total_component_concentration(S->get_total_component_concentration()),
     porosity(S->get_porosity()),
     darcy_flux(S->get_darcy_flux()),
     water_saturation(S->get_water_saturation()),
     mesh_maps(S->get_mesh_maps()) {};

  Transport_State (State S):
     total_component_concentration(S.get_total_component_concentration()),
     porosity(S.get_porosity()),
     darcy_flux(S.get_darcy_flux()),
     water_saturation(S.get_water_saturation()),
     mesh_maps(S.get_mesh_maps()) {};

  Transport_State () {};

  ~Transport_State () {};

  /* access methods for state variables */
  RCP<Epetra_MultiVector> get_total_component_concentration()       { return total_component_concentration; };
  RCP<Epetra_MultiVector> get_total_component_concentration() const { return total_component_concentration; };
  RCP<const Epetra_Vector>      get_porosity ()               const { return porosity; };
  RCP<const Epetra_Vector>      get_water_saturation ()       const { return water_saturation; };
  RCP<const Epetra_Vector>      get_darcy_flux ()             const { return darcy_flux; };
  
  void copy (RCP<Transport_State> TS);

private:
  /* state variables that are relevant to transport */
  RCP<Epetra_MultiVector>   total_component_concentration;
  RCP<const Epetra_Vector>  water_saturation;
  RCP<const Epetra_Vector>  darcy_flux;
  RCP<const Epetra_Vector>  porosity;

  /* mesh infranstructure */
  RCP<const Mesh_maps_simple> mesh_maps;
};


#endif

