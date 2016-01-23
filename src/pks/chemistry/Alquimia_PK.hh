/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_ALQUIMIA_CHEMISTRY_PK_HH_
#define AMANZI_ALQUIMIA_CHEMISTRY_PK_HH_

#include <string>
#include <vector>
#include <map>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Chemistry_PK_Base.hh"
#include "ChemistryEngine.hh"
#include "VerboseObject.hh"

// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_SerialDenseVector;

namespace Amanzi {
namespace AmanziChemistry {

// Trilinos based chemistry process kernel for the unstructured mesh
class Alquimia_PK: public Chemistry_PK_Base {
 public:

  // Constructor. Note that we must pass the "Main" parameter list
  // to this PK so that it has access to all information about the 
  // problem.
  Alquimia_PK(const Teuchos::ParameterList& param_list,
              Teuchos::RCP<Chemistry_State> chem_state,
              Teuchos::RCP<ChemistryEngine> chem_engine);

  ~Alquimia_PK();

  void InitializeChemistry(void);

  void Advance(const double& delta_time,
               Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star);
  void CommitState(Teuchos::RCP<Chemistry_State> chem_state, const double& time);
  Teuchos::RCP<Epetra_MultiVector> get_total_component_concentration(void) const;

  double time_step(void) const {
    return this->time_step_;
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

  // Ben: the following routine provides the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> get_extra_chemistry_output_data();

 protected:

 private:

  // Timestepping controls.
  double time_step_, max_time_step_, min_time_step_, prev_time_step_;
  std::string time_step_control_method_;
  int num_iterations_for_time_step_cut_, num_steps_before_time_step_increase_;
  double time_step_cut_factor_, time_step_increase_factor_;
  int num_iterations_, num_successful_steps_;
  void ComputeNextTimeStep();

  // auxilary state for process kernel
  Teuchos::RCP<Chemistry_State> chemistry_state_;

  // parameter lists
  Teuchos::ParameterList main_param_list_, chem_param_list_;

  bool chem_initialized_;

  // Chemistry engine.
  Teuchos::RCP<ChemistryEngine> chem_engine_;

  // Alquimia data structures.
  AlquimiaState alq_state_;
  AlquimiaMaterialProperties alq_mat_props_;
  AlquimiaAuxiliaryData alq_aux_data_;
  AlquimiaAuxiliaryOutputData alq_aux_output_;

  // Mapping of region names to geochemical condition names. A region is identified 
  // by a string, and all cells within a region will have all geochemical 
  // conditions in the corresponding condition vector applied to them. 
  std::map<std::string, std::string> chem_initial_conditions_;
  
  double current_time_;
  double saved_time_;

  // Auxiliary output data, requested by and stored within Amanzi.
  std::vector<std::string> aux_names_;
  Teuchos::RCP<Epetra_MultiVector> aux_output_;

  // For printing diagnostic information.
  Teuchos::RCP<VerboseObject> vo_;

  void UpdateChemistryStateStorage(void);
  int InitializeSingleCell(int cell_index, const std::string& condition);
  int AdvanceSingleCell(double delta_time, 
                        Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star,
                        int cell_index);

  void ParseChemicalConditionRegions(const Teuchos::ParameterList& param_list,
                                     std::map<std::string, std::string>& conditions);
  void XMLParameters(void);

  // These helpers copy data back and forth between a set of buffers and the chemistry state.
  // given cell.
  void CopyAmanziStateToAlquimia(const int cell_id,
                                 Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                                 AlquimiaMaterialProperties& mat_props,
                                 AlquimiaState& state,
                                 AlquimiaAuxiliaryData& aux_data);

  void CopyAlquimiaStateToAmanzi(const int cell_id,
                                 const AlquimiaMaterialProperties& mat_props,
                                 const AlquimiaState& state,
                                 const AlquimiaAuxiliaryData& aux_data,
                                 const AlquimiaAuxiliaryOutputData& aux_output);

  void InitAmanziStateFromAlquimia(const int cell_id,
                                   const AlquimiaMaterialProperties& mat_props,
                                   const AlquimiaState& state,
                                   const AlquimiaAuxiliaryData& aux_data,
                                   const AlquimiaAuxiliaryOutputData& aux_output);

};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif  // AMANZI_ALQUIMIA_CHEMISTRY_PK_HH_
