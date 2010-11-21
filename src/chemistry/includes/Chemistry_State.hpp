/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Chemistry_State_hpp__
#define __Chemistry_State_hpp__

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"

class Chemistry_State {

 public:
  Chemistry_State (Teuchos::RCP<State> S);

  ~Chemistry_State();

  // access methods
  Teuchos::RCP<const Epetra_MultiVector> get_total_component_concentration()
  { return total_component_concentration_; };

  Teuchos::RCP<const Epetra_Vector> get_porosity() const { return porosity_; };
  Teuchos::RCP<const Epetra_Vector> get_water_saturation() const { return water_saturation_; };
  Teuchos::RCP<const Epetra_Vector> get_water_density() const { return water_density_; };

  Teuchos::RCP<const Mesh_maps_base> get_mesh_maps() const { return mesh_maps_; };

  Teuchos::RCP<const Epetra_SerialDenseVector> get_volume() const { return volume_; };

 private:
  // variables that point to main State object
  Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_;
  Teuchos::RCP<const Epetra_Vector> porosity_;
  Teuchos::RCP<const Epetra_Vector> water_saturation_;
  Teuchos::RCP<const Epetra_Vector> water_density_;

  Teuchos::RCP<Mesh_maps_base> mesh_maps_;

  // local variable
  Teuchos::RCP<Epetra_SerialDenseVector> volume_;

  void ExtractVolumeFromMesh(void);

};



#endif
