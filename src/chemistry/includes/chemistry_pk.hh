/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_PK_HH_
#define AMANZI_CHEMISTRY_PK_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "beaker.hh"
#include "chemistry_exception.hh"
#include "verbosity.hh"

// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_SerialDenseVector;

namespace amanzi {
namespace chemistry {

// forward declarations from chemistry
class Chemistry_State;

// Trilinos based chemistry process kernel for the unstructured mesh
class Chemistry_PK {
 public:


  Chemistry_PK(const Teuchos::ParameterList& param_list,
               Teuchos::RCP<Chemistry_State> chem_state);

  ~Chemistry_PK();

  void InitializeChemistry(void);

  void advance(const double& delta_time,
               Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star);
  void commit_state(Teuchos::RCP<Chemistry_State> chem_state, const double& delta_time);
  Teuchos::RCP<Epetra_MultiVector> get_total_component_concentration(void) const;


  Verbosity verbosity(void) const {
    return this->verbosity_;
  }
  void set_verbosity(const Verbosity verbosity) {
    this->verbosity_ = verbosity;
  }

  void set_max_time_step(const double mts) {
    this->max_time_step_ = mts;
  }
  double max_time_step(void) const {
    return this->max_time_step_;
  }

  int number_aqueous_components(void) const {
    return this->number_aqueous_components_;
  }
  void set_number_aqueous_components(const int nac) {
    this->number_aqueous_components_ = nac;
  }

  int have_free_ion_guess(void) const {
    return this->have_free_ion_guess_;
  }
  void set_have_free_ion_guess(const int hfi) {
    this->have_free_ion_guess_ = hfi;
  }

  int number_free_ion(void) const {
    return this->number_free_ion_;
  }
  void set_number_free_ion(const int nfi) {
    this->number_free_ion_ = nfi;
  }

  int number_total_sorbed(void) const {
    return this->number_total_sorbed_;
  }
  void set_number_total_sorbed(const int nts) {
    this->number_total_sorbed_ = nts;
  }

  int number_minerals(void) const {
    return this->number_minerals_;
  }
  void set_number_minerals(const int nm) {
    this->number_minerals_ = nm;
  }

  int number_ion_exchange_sites(void) const {
    return this->number_ion_exchange_sites_;
  }
  void set_number_ion_exchange_sites(const int ies) {
    this->number_ion_exchange_sites_ = ies;
  }

  int number_sorption_sites(void) const {
    return this->number_sorption_sites_;
  }
  void set_number_sorption_sites(const int scs) {
    this->number_sorption_sites_ = scs;
  }

  int using_sorption(void) const {
    return this->using_sorption_;
  }
  void set_using_sorption(const int us) {
    this->using_sorption_ = us;
  }

  // Ben: the following two routines provide the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> get_extra_chemistry_output_data();
  void set_chemistry_output_names(std::vector<std::string>* names);

  // Ben: this routine should set the strings that will be
  // appended to the component_x tag in the cgns output
  void set_component_names(std::vector<std::string>* names);

 protected:

 private:
  Verbosity verbosity_;
  double max_time_step_;
  // auxilary state for process kernel
  Teuchos::RCP<Chemistry_State> chemistry_state_;

  // parameter list
  Teuchos::ParameterList parameter_list_;

  Beaker* chem_;
  Beaker::BeakerParameters beaker_parameters_;
  Beaker::BeakerComponents beaker_components_;
  Beaker::BeakerComponents beaker_components_copy_;

  double current_time_;
  double saved_time_;
  int number_aqueous_components_;
  int number_free_ion_;
  int number_total_sorbed_;
  int number_minerals_;
  int number_ion_exchange_sites_;
  int number_sorption_sites_;
  int using_sorption_;
  int have_free_ion_guess_;
  std::vector<std::string> aux_names_;
  std::vector<int> aux_index_;

  struct InternalStorage {
    // things we don't change, just point to State object
    Teuchos::RCP<const Epetra_Vector> porosity;
    Teuchos::RCP<const Epetra_Vector> water_saturation;
    Teuchos::RCP<const Epetra_Vector> water_density;
    Teuchos::RCP<const Epetra_SerialDenseVector> volume;
    // things we need internal copies of
    Teuchos::RCP<Epetra_MultiVector> aqueous_components;
    Teuchos::RCP<Epetra_MultiVector> free_ion_species;
    Teuchos::RCP<Epetra_MultiVector> minerals;
    Teuchos::RCP<Epetra_MultiVector> ion_exchange_sites;
    Teuchos::RCP<Epetra_MultiVector> sorption_sites;
    Teuchos::RCP<Epetra_MultiVector> total_sorbed;
    // geh can do without for now.    Teuchos::RCP<Epetra_MultiVector> free_site_concentrations;
  };

  Teuchos::RCP<Epetra_MultiVector> aux_data_;

  void InitializeInternalStorage(InternalStorage* storage);
  void SwapCurrentAndSavedStorage(void);

  InternalStorage current_state_;
  InternalStorage saved_state_;

  void XMLParameters(void);
  void LocalPhysicalState(void);
  void LocalInitialConditions(void);
  void SetupAuxiliaryOutput(void);
  void ExtractInitialCondition(const std::string& type,
                               const std::string& keyword,
                               const int number_to_find,
                               const int block,
                               const Teuchos::ParameterList& mesh_block_list,
                               const int mesh_block_ID,
                               Teuchos::RCP<Epetra_MultiVector> data);
  void set_const_values_for_block(const std::vector<double>& values,
                                  const int num_values,
                                  Teuchos::RCP<Epetra_MultiVector>& vec,
                                  const int mesh_block_id);
  void set_cell_value_in_mesh_block(const double value,
                                    Epetra_Vector& vec,
                                    const int mesh_block_id);
  void SizeBeakerComponents(void);
  void CopyCellToBeakerComponents(
      const int cell_id,
      Teuchos::RCP<const Epetra_MultiVector> aqueous_components);
  void CopyBeakerComponentsToCell(const int cell_id);
  void CopyStateToBeakerParameters(const int cell_id);
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_PK_HH_
