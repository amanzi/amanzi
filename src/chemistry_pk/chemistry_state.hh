/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_CHEMISTRY_STATE_HH_
#define AMANZI_CHEMISTRY_CHEMISTRY_STATE_HH_

#include "State.hh"
#include "Mesh.hh"

#include "Teuchos_RCP.hpp"

#include "beaker.hh"

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

  void AllocateAdditionalChemistryStorage(
      const amanzi::chemistry::Beaker::BeakerComponents&);

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

  Teuchos::RCP<Epetra_MultiVector> primary_activity_coeff() const {
    return simulation_state_->primary_activity_coeff();
  }

  Teuchos::RCP<Epetra_MultiVector> secondary_activity_coeff() const {
    return simulation_state_->secondary_activity_coeff();
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

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_ref_cation_conc() const {
    return simulation_state_->ion_exchange_ref_cation_conc();
  }

  int number_of_sorption_sites(void) const {
    return simulation_state_->number_of_sorption_sites();
  }

  Teuchos::RCP<Epetra_MultiVector> sorption_sites() const {
    return simulation_state_->sorption_sites();
  }

  Teuchos::RCP<Epetra_MultiVector> surface_complex_free_site_conc() const {
    return simulation_state_->surface_complex_free_site_conc();
  }

  bool using_sorption(void) const {
    return simulation_state_->using_sorption();
  }

  Teuchos::RCP<Epetra_MultiVector> isotherm_kd() const {
    return simulation_state_->isotherm_kd();
  }

  Teuchos::RCP<Epetra_MultiVector> isotherm_freundlich_n() const {
    return simulation_state_->isotherm_freundlich_n();
  }

  Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b() const {
    return simulation_state_->isotherm_langmuir_b();
  }

  bool using_sorption_isotherms(void) const {
    return simulation_state_->use_sorption_isotherms();
  }

 private:
  Teuchos::RCP<State> simulation_state_;

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_CHEMISTRY_STATE_HH_
