#ifndef __Transport_State_hpp__
#define __Transport_State_hpp__

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh_maps_simple.hh"



class Transport_State {

public:
  Transport_State (Teuchos::RCP<State> S):
     total_component_concentration(S->get_total_component_concentration()),
     porosity(S->get_porosity()),
     darcy_flux(S->get_darcy_flux()),
     water_saturation(S->get_water_saturation()),
     mesh_maps(S->get_mesh_maps()) {};

  ~Transport_State () {};

  /* access methods for state variables */
  Teuchos::RCP<const Epetra_MultiVector> get_total_component_concentration()       { return total_component_concentration; };
  Teuchos::RCP<const Epetra_MultiVector> get_total_component_concentration() const { return total_component_concentration; };
  Teuchos::RCP<const Epetra_Vector>      get_porosity ()                     const { return porosity; };
  Teuchos::RCP<const Epetra_Vector>      get_water_saturation ()             const { return water_saturation; };
  Teuchos::RCP<const Epetra_Vector>      get_darcy_flux ()                   const { return darcy_flux; };
  
  void copy (Teuchos::RCP<Transport_State> TS);

private:
  /* state variables that are relevant to transport */
  Teuchos::RCP<const Epetra_MultiVector>  total_component_concentration;
  Teuchos::RCP<const Epetra_Vector>       water_saturation;
  Teuchos::RCP<const Epetra_Vector>       darcy_flux;
  Teuchos::RCP<const Epetra_Vector>       porosity;

  /* mesh infranstructure */
  Teuchos::RCP<const Mesh_maps_simple> mesh_maps;
};


#endif

