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

  // stuff chemistry can't/shouldn't change
  Teuchos::RCP<const Mesh> mesh_maps() const {
    return simulation_state_->get_mesh_maps();
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

  Teuchos::RCP<const Epetra_Vector> volume() const {
    return simulation_state_->volume();
  }

  int number_of_aqueous_components(void) const {
    return simulation_state_->get_number_of_components();
  }

  // stuff chemistry can change
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration() {
    return simulation_state_->get_total_component_concentration();
  }

  Teuchos::RCP<Epetra_MultiVector> free_ion_species() const {
    return simulation_state_->free_ion_concentrations();
  }

  int number_of_minerals(void) const {
    return simulation_state_->number_of_minerals();
  }

  Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions() const {
    return simulation_state_->mineral_volume_fractions();
  }

  Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area() const {
    return simulation_state_->mineral_specific_surface_area();
  }

  Teuchos::RCP<Epetra_MultiVector> total_sorbed() const {
    return simulation_state_->total_sorbed();
  }

  int number_of_ion_exchange_sites(void) const {
    return simulation_state_->number_of_ion_exchange_sites();
  }

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites() const {
    return simulation_state_->ion_exchange_sites();
  }

  int number_of_sorption_sites(void) const {
    return simulation_state_->number_of_sorption_sites();
  }

  Teuchos::RCP<Epetra_MultiVector> sorption_sites() const {
    return simulation_state_->sorption_sites();
  }

  bool using_sorption(void) const {
    return simulation_state_->using_sorption();
  }

 private:
  Teuchos::RCP<State> simulation_state_;

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_CHEMISTRY_STATE_HH_
