#ifndef __Flow_State_hpp__
#define __Flow_State_hpp__

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Mesh_maps_simple.hh"
#include "State.hpp"

class Flow_State {

public:
  Flow_State (Teuchos::RCP<State> S):
    total_component_concentration(S->get_total_component_concentration()),
    porosity(S->get_porosity()),
    water_density(S->get_water_density()),
    water_saturation(S->get_water_saturation()),
    mesh_maps(S->get_mesh_maps())
  { };

  ~Flow_State () {};

  // access methods
  Teuchos::RCP<const Epetra_MultiVector> get_total_component_concentration() 
  { return total_component_concentration; }; 
  
  Teuchos::RCP<const Epetra_Vector> get_porosity () const { return porosity; };
  Teuchos::RCP<const Epetra_Vector> get_water_saturation () const { return water_saturation; };
  Teuchos::RCP<const Epetra_Vector> get_water_density () const { return water_density; };

  const Teuchos::RCP<const Mesh_maps_simple> get_mesh_maps() const { return mesh_maps;};

private:
  // variables that are relevant to chemistry
  Teuchos::RCP<const Epetra_MultiVector> total_component_concentration;
  Teuchos::RCP<const Epetra_Vector> porosity;
  Teuchos::RCP<const Epetra_Vector> water_saturation;
  Teuchos::RCP<const Epetra_Vector> water_density;
  
  const Teuchos::RCP<const Mesh_maps_simple> mesh_maps;
};



#endif
