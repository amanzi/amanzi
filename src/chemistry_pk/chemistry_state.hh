/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_CHEMISTRY_STATE_HH_
#define AMANZI_CHEMISTRY_CHEMISTRY_STATE_HH_

#include "State.hpp"
#include "Mesh.hh"

#include "Teuchos_RCP.hpp"

// forward declarations
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_SerialDenseVector;
class ParameterList;

namespace amanzi {
namespace chemistry {

using Amanzi::AmanziMesh::Mesh;

class Chemistry_State {
 public:
  explicit Chemistry_State(Teuchos::RCP<State> S);

  virtual ~Chemistry_State();

  void AllocateMemory(const int num_aqueous,
                      const int num_free_ion,
                      const int num_minerals,
                      const int num_ion_exchange_sites,
                      const int num_total_sorbed,
                      const int num_sorption_sites);

  // access methods
  Teuchos::RCP<const Epetra_MultiVector> total_component_concentration() {
    return simulation_state_->get_total_component_concentration();
  }

  Teuchos::RCP<const Epetra_Vector> porosity() const {
    return simulation_state_->get_porosity();
  }
  Teuchos::RCP<const Epetra_Vector> water_saturation() const {
    return simulation_state_->get_water_saturation();
  }
  Teuchos::RCP<const Epetra_Vector> water_density() const {
    return simulation_state_->get_water_density();
  }

  Teuchos::RCP<const Mesh> mesh_maps() const {
    return simulation_state_->get_mesh_maps();
  }

  Teuchos::RCP<const Epetra_Vector> volume() const {
    return simulation_state_->volume();
  }

  Teuchos::RCP<Epetra_MultiVector> aqueous_components() const {
    return aqueous_components_;
  }

  Teuchos::RCP<Epetra_MultiVector> free_ion_species() const {
    return free_ion_species_;
  }

  Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions() const {
    return mineral_volume_fractions_;
  }

  Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area() const {
    return mineral_specific_surface_area_;
  }

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites() const {
    return ion_exchange_sites_;
  }

  Teuchos::RCP<Epetra_MultiVector> sorption_sites() const {
    return sorption_sites_;
  }

  Teuchos::RCP<Epetra_MultiVector> total_sorbed() const {
    return total_sorbed_;
  }

 private:
  Teuchos::RCP<State> simulation_state_;

  // local variable
  Teuchos::RCP<Epetra_Vector> volume_;
  Teuchos::RCP<Epetra_MultiVector> aqueous_components_;
  Teuchos::RCP<Epetra_MultiVector> free_ion_species_;
  Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions_;
  Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area_;
  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites_;
  Teuchos::RCP<Epetra_MultiVector> sorption_sites_;
  Teuchos::RCP<Epetra_MultiVector> total_sorbed_;

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_CHEMISTRY_STATE_HH_
