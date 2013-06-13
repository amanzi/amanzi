/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Chemistry

License: see COPYRIGHT
Author: Ethan Coon

Interface layer between Chemistry_PK and State, this is a harness for
accessing the new state-dev from the old Flow PK.

 ------------------------------------------------------------------------- */

#ifndef AMANZI_CHEMISTRY_STATE_NEW_HH_
#define AMANZI_CHEMISTRY_STATE_NEW_HH_

#include "PK_State.hh"
#include "beaker.hh"

namespace Amanzi {
namespace AmanziChemistry {


class Chemistry_State : public PK_State {

 public:

  Chemistry_State(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<State>& S);

  Chemistry_State(const Teuchos::RCP<State>& S,
                  int number_of_aqueous_components,
                  int number_of_minerals,
                  int number_of_ion_exchange_sites,
                  int number_of_sorption_sites,
                  bool using_sorption,
                  bool using_sorption_isotherms);

  virtual ~Chemistry_State() {}

  void AllocateAdditionalChemistryStorage(const Beaker::BeakerComponents&);

  void Initialize();

  // accessors
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_maps() const {
    return S_->GetMesh();
  }

  // access methods for state variables
  // const methods
  Teuchos::RCP<const Epetra_MultiVector> total_component_concentration() const {
    return S_->GetFieldData("total_component_concentration")->ViewComponent("cell", true);
  }

  Teuchos::RCP<const Epetra_Vector> porosity() const {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("porosity")->ViewComponent("cell", ghosted_))(0));
  }

  Teuchos::RCP<const Epetra_Vector> water_saturation() const {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("water_saturation")->ViewComponent("cell", ghosted_))(0));
  }

  Teuchos::RCP<const Epetra_Vector> water_density() const {
    if (!water_density_initialized_) {
      water_density_->PutScalar(*S_->GetScalarData("fluid_density"));
      water_density_initialized_ = true;
    }
    return water_density_;
  }

  Teuchos::RCP<const Epetra_Vector> volume() const {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("cell_volume")->ViewComponent("cell", ghosted_))(0));
  }

  int number_of_aqueous_components(void) const {
    return number_of_aqueous_components_;
  }


  // non-const accessors
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration() {
    return S_->GetFieldData("total_component_concentration", name_)
        ->ViewComponent("cell", true);
  }

  Teuchos::RCP<Epetra_MultiVector> free_ion_species() {
    return S_->GetFieldData("free_ion_species", name_)
        ->ViewComponent("cell", true);
  }

  Teuchos::RCP<Epetra_MultiVector> primary_activity_coeff() {
    return S_->GetFieldData("primary_activity_coeff", name_)
        ->ViewComponent("cell", true);
  }

  Teuchos::RCP<Epetra_MultiVector> secondary_activity_coeff() {
    return S_->GetFieldData("secondary_activity_coeff", name_)
        ->ViewComponent("cell", true);
  }


  int number_of_minerals(void) const {
    return number_of_minerals_;
  }

  Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions() {
    return S_->GetFieldData("mineral_volume_fractions", name_)
        ->ViewComponent("cell", true);
  }

  Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area() {
    return S_->GetFieldData("mineral_specific_surface_area", name_)
        ->ViewComponent("cell", true);
  }

  Teuchos::RCP<Epetra_MultiVector> total_sorbed() {
    return S_->GetFieldData("total_sorbed", name_)
        ->ViewComponent("cell", true);
  }

  int number_of_ion_exchange_sites(void) const {
    return number_of_ion_exchange_sites_;
  }

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites() {
    return S_->GetFieldData("ion_exchange_sites", name_)
        ->ViewComponent("cell", true);
  }

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_ref_cation_conc() {
    return S_->GetFieldData("ion_exchange_ref_cation_conc", name_)
        ->ViewComponent("cell", true);
  }

  int number_of_sorption_sites(void) const {
    return number_of_sorption_sites_;
  }

  Teuchos::RCP<Epetra_MultiVector> sorption_sites() {
    return S_->GetFieldData("sorption_sites", name_)
        ->ViewComponent("cell", true);
  }

  Teuchos::RCP<Epetra_MultiVector> surface_complex_free_site_conc() {
    return S_->GetFieldData("surface_complex_free_site_conc", name_)
        ->ViewComponent("cell", true);
  }


  bool using_sorption() const {
    return using_sorption_;
  }


  Teuchos::RCP<Epetra_MultiVector> isotherm_kd() {
    return S_->GetFieldData("isotherm_kd", name_)
        ->ViewComponent("cell", true);
  }

  Teuchos::RCP<Epetra_MultiVector> isotherm_freundlich_n() {
    return S_->GetFieldData("isotherm_freundlich_n", name_)
        ->ViewComponent("cell", true);
  }

    Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b() {
    return S_->GetFieldData("isotherm_langmuir_b", name_)
        ->ViewComponent("cell", true);
  }


  bool using_sorption_isotherms() const {
    return using_sorption_isotherms_;
  }

 protected:

  void InitializeField_(Teuchos::ParameterList& ic_plist, std::string fieldname,
                        bool sane_default, double default_val);

  void SetupSoluteNames_();
  void SetupMineralNames_();
  void SetupSorptionSiteNames_();
  void ParseMeshBlocks_();
  void VerifyMineralogy_(const std::string& region_name,
                         const Teuchos::ParameterList& minerals_list);
  void VerifySorptionIsotherms_(const std::string& region_name,
          const Teuchos::ParameterList& isotherms_list);
  void VerifySorptionSites_(const std::string& region_name,
                            const Teuchos::ParameterList& sorption_site_list);
  void RequireData_();

 private:
  // not implemented
  Chemistry_State(const Chemistry_State& other);

 protected:
  Teuchos::ParameterList plist_;

  int number_of_aqueous_components_;
  int number_of_minerals_;
  int number_of_ion_exchange_sites_;
  int number_of_sorption_sites_;
  bool using_sorption_;
  bool using_sorption_isotherms_;

  std::vector<std::string> compnames_;
  std::map<std::string,int> comp_name_id_map_;
  std::vector<std::string> mineral_names_;
  std::map<std::string, int> mineral_name_id_map_;
  std::vector<std::string> sorption_site_names_;
  std::map<std::string, int> sorption_site_name_id_map_;

  mutable bool water_density_initialized_;
  Teuchos::RCP<Epetra_Vector> water_density_;

};
} // namespace
} // namespace

#endif
