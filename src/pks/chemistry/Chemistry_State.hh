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

#include "State.hh"
#include "beaker.hh"

#ifdef ALQUIMIA_ENABLED
#include "Teuchos_RCP.hpp"
#include "ChemistryEngine.hh"
#endif

namespace Amanzi {
namespace AmanziChemistry {


class Chemistry_State {

 public:

  Chemistry_State(Teuchos::ParameterList& plist,
                  const std::vector<std::string>& component_names,
                  const Teuchos::RCP<State>& S);

  virtual ~Chemistry_State() {};

  void Setup();
  void AllocateAdditionalChemistryStorage(const Beaker::BeakerComponents&);
  void AllocateAdditionalChemistryStorage(int num_aqueous_components);
  void SetAuxDataNames(const std::vector<std::string>& aux_data_names);

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
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("saturation_liquid")->ViewComponent("cell", ghosted_))(0));
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
    try {
      return S_->GetFieldData("free_ion_species", name_)
        ->ViewComponent("cell", true);
    } catch (...) {
      return Teuchos::null;
    }
  }

  Teuchos::RCP<Epetra_MultiVector> primary_activity_coeff() {
    try {
      return S_->GetFieldData("primary_activity_coeff", name_)
        ->ViewComponent("cell", true);
    } catch (...) {
      return Teuchos::null;
    } 
  }

  Teuchos::RCP<Epetra_MultiVector> secondary_activity_coeff() {
    try {
      return S_->GetFieldData("secondary_activity_coeff", name_)
        ->ViewComponent("cell", true);
    } catch (...) {
      return Teuchos::null;
    }
  }


  int number_of_minerals(void) const {
    return number_of_minerals_;
  }

  Teuchos::RCP<Epetra_MultiVector> mineral_volume_fractions() {
    if (number_of_minerals_ > 0) {
      return S_->GetFieldData("mineral_volume_fractions", name_)
        ->ViewComponent("cell", true);
    }  else {
      return Teuchos::null;
    }
  }

  Teuchos::RCP<Epetra_MultiVector> mineral_specific_surface_area() {
    if (number_of_minerals_ > 0) {
      return S_->GetFieldData("mineral_specific_surface_area", name_)
        ->ViewComponent("cell", true);
    } else {
      return Teuchos::null;
    }
  }

  Teuchos::RCP<Epetra_MultiVector> total_sorbed() {
    if (using_sorption_ > 0) {    
      return S_->GetFieldData("total_sorbed", name_)
        ->ViewComponent("cell", true);
    } else {
      return Teuchos::null;
    }
  }

  int number_of_ion_exchange_sites(void) const {
    return number_of_ion_exchange_sites_;
  }

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites() {
    if (number_of_ion_exchange_sites_ > 0) {
      return S_->GetFieldData("ion_exchange_sites", name_)
        ->ViewComponent("cell", true);
    } else {
      return Teuchos::null;
    }
  }

  Teuchos::RCP<Epetra_MultiVector> ion_exchange_ref_cation_conc() {
    if (number_of_ion_exchange_sites_ > 0) {
      return S_->GetFieldData("ion_exchange_ref_cation_conc", name_)
        ->ViewComponent("cell", true);
    } else {
      return Teuchos::null;
    } 
  }

  int number_of_sorption_sites(void) const {
    return number_of_sorption_sites_;
  }

  Teuchos::RCP<Epetra_MultiVector> sorption_sites() {
    if (number_of_sorption_sites_ > 0) {
      return S_->GetFieldData("sorption_sites", name_)
        ->ViewComponent("cell", true);
    } else {
      return Teuchos::null;
    }
  }

  Teuchos::RCP<Epetra_MultiVector> surface_complex_free_site_conc() {
    if (number_of_sorption_sites_ > 0) {    
      return S_->GetFieldData("surface_complex_free_site_conc", name_)
        ->ViewComponent("cell", true);
    } else {
      return Teuchos::null;
    } 
  }

  // accessor of aux data
  Teuchos::RCP<Epetra_MultiVector> aux_data(std::string &auxname) {
    return S_->GetFieldData(auxname, name_)
      ->ViewComponent("cell", true);
  }  



  bool using_sorption() const {
    return using_sorption_;
  }


  Teuchos::RCP<Epetra_MultiVector> isotherm_kd() {
    if (using_sorption_isotherms_) {
      return S_->GetFieldData("isotherm_kd", name_)
        ->ViewComponent("cell", true);
    } else {
      return Teuchos::null;
    }
  }

  Teuchos::RCP<Epetra_MultiVector> isotherm_freundlich_n() {
    if (using_sorption_isotherms_) {
      return S_->GetFieldData("isotherm_freundlich_n", name_)
        ->ViewComponent("cell", true);
    } else {
      return Teuchos::null;
    }
  }

  Teuchos::RCP<Epetra_MultiVector> isotherm_langmuir_b() {
    if (using_sorption_isotherms_) {
      return S_->GetFieldData("isotherm_langmuir_b", name_)
        ->ViewComponent("cell", true);
    } else {
      return Teuchos::null;
    }
  }


  bool using_sorption_isotherms() const {
    return using_sorption_isotherms_;
  }


  
  // Set the solution component names.
//  void SetComponentNames(const std::vector<std::string>& comp_names);

#ifdef ALQUIMIA_ENABLED
  // The following methods are designed to help the Chemistry State interact with 
  // Alquimia via the ChemistryEngine.

  // Copies the chemistry state in the given cell to the given Alquimia containers.
  void CopyToAlquimia(const int cell_id,
                      AlquimiaMaterialProperties& mat_props,
                      AlquimiaState& state,
                      AlquimiaAuxiliaryData& aux_data);

  // Copies the chemistry state in the given cell to the given Alquimia containers, 
  // taking the aqueous components from the given multivector.
  void CopyToAlquimia(const int cell_id,
                      Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                      AlquimiaMaterialProperties& mat_props,
                      AlquimiaState& state,
                      AlquimiaAuxiliaryData& aux_data);

  // Copies the data in the given Alquimia containers to the given cell within the 
  // chemistry state. The aqueous component concentrations are placed into 
  // the aqueous_components multivector.
  void CopyFromAlquimia(const int cell_id,
                        const AlquimiaMaterialProperties& mat_props,
                        const AlquimiaState& state,
                        const AlquimiaAuxiliaryData& aux_data,
                        const AlquimiaAuxiliaryOutputData& aux_output,
                        Teuchos::RCP<const Epetra_MultiVector> aqueous_components);

  // Copies the data in the given Alquimia containers to the given cell within the 
  // chemistry state. Modifies only uninitialize fields

  void InitFromAlquimia(const int cell_id,
                        const AlquimiaMaterialProperties& mat_props,
                        const AlquimiaState& state,
                        const AlquimiaAuxiliaryData& aux_data,
                        const AlquimiaAuxiliaryOutputData& aux_output);

  void SetAllFieldsInitialized();

#endif

 protected:
  void InitializeField_(Teuchos::ParameterList& ic_plist, std::string fieldname,
                        bool sane_default, double default_val);

  void SetupMineralNames_();
//  void SetupSoluteNames_();
  void SetupSorptionSiteNames_();
  void ParseMeshBlocks_();
  void VerifyMineralogy_(const std::string& region_name,
                         const Teuchos::ParameterList& minerals_list);
  void VerifySorptionIsotherms_(const std::string& region_name,
          const Teuchos::ParameterList& isotherms_list);
  void VerifySorptionSites_(const std::string& region_name,
                            const Teuchos::ParameterList& sorption_site_list);
  void RequireData_();
  void RequireAuxData_();
  

 private:
  // not implemented
  Chemistry_State(const Chemistry_State& other);
  Chemistry_State& operator=(const Chemistry_State& other);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<State> S_;
  bool ghosted_;
  std::string name_;

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

  // Auxiliary data, maintained by Amanzi and updated
  int num_aux_data_;
  Teuchos::RCP<Epetra_MultiVector> aux_data_;


};
} // namespace
} // namespace

#endif
