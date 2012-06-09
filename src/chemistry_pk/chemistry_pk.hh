/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_PK_HH_
#define AMANZI_CHEMISTRY_PK_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "beaker.hh"
#include "chemistry_exception.hh"
#include "chemistry_verbosity.hh"

#include "chemistry_state.hh"

// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_SerialDenseVector;

namespace amanzi {
namespace chemistry {

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

  void set_max_time_step(const double mts) {
    this->max_time_step_ = mts;
  }
  double max_time_step(void) const {
    return this->max_time_step_;
  }

  int number_aqueous_components(void) const {
    return chemistry_state_->number_of_aqueous_components();
  }

  int number_free_ion(void) const {
    return chemistry_state_->number_of_aqueous_components();
  }

  int number_total_sorbed(void) const {
    return chemistry_state_->number_of_aqueous_components();
  }

  int number_minerals(void) const {
    return chemistry_state_->number_of_minerals();
  }

  int number_ion_exchange_sites(void) const {
    return chemistry_state_->number_of_ion_exchange_sites();
  }

  int number_sorption_sites(void) const {
    return chemistry_state_->number_of_sorption_sites();
  }

  int using_sorption(void) const {
    return chemistry_state_->using_sorption();
  }

  int using_sorption_isotherms(void) const {
    return chemistry_state_->using_sorption_isotherms();
  }

  bool debug(void) const {
    return debug_;
  }

  void set_debug(const bool value) {
    debug_ = value;
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
  bool debug_;
  bool display_free_columns_;
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

  std::vector<std::string> aux_names_;
  std::vector<int> aux_index_;

  Teuchos::RCP<Epetra_MultiVector> aux_data_;

  void UpdateChemistryStateStorage(void);

  void XMLParameters(void);
  void SetupAuxiliaryOutput(void);
  void SizeBeakerStructures(void);
  void CopyCellStateToBeakerStructures(
      const int cell_id,
      Teuchos::RCP<const Epetra_MultiVector> aqueous_components);
  void CopyBeakerStructuresToCellState(const int cell_id);
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_PK_HH_
