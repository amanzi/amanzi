/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_ALQUIMIA_CHEMISTRY_PK_HH_
#define AMANZI_ALQUIMIA_CHEMISTRY_PK_HH_

#include <string>
#include <vector>
#include <map>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "chemistry_pk_base.hh"
#include "ChemistryEngine.hh"

// forward declarations
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_SerialDenseVector;

namespace Amanzi {
namespace AmanziChemistry {

// Trilinos based chemistry process kernel for the unstructured mesh
class Alquimia_Chemistry_PK: public Chemistry_PK_Base {
 public:

  // Constructor. Note that we must pass the "Main" parameter list
  // to this PK so that it has access to all information about the 
  // problem.
  Alquimia_Chemistry_PK(const Teuchos::ParameterList& param_list,
                        Teuchos::RCP<Chemistry_State> chem_state,
                        Teuchos::RCP<ChemistryEngine> chem_engine);

  ~Alquimia_Chemistry_PK();

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

  // Ben: the following routine provides the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> get_extra_chemistry_output_data();

 protected:

 private:
  bool debug_;
  bool display_free_columns_;
  double max_time_step_;
  // auxilary state for process kernel
  Teuchos::RCP<Chemistry_State> chemistry_state_;

  // parameter lists
  Teuchos::ParameterList main_param_list_, chem_param_list_;

  bool chem_initialized_;

  // Chemistry engine.
  Teuchos::RCP<ChemistryEngine> chem_engine_;

  // Alquimia data structures.
  AlquimiaState* alq_state_;
  AlquimiaMaterialProperties* alq_mat_props_;
  AlquimiaAuxiliaryData* alq_aux_data_;
  AlquimiaAuxiliaryOutputData* alq_aux_output_;

  // Mapping of region names to geochemical condition names. A region is identified 
  // by a string, and all cells within a region will have all geochemical 
  // conditions in the corresponding condition vector applied to them. 
  std::map<std::string, std::string> chem_initial_conditions_;
  
  double current_time_;
  double saved_time_;

  // Auxiliary output data, requested by and stored within Amanzi.
  std::vector<std::string> aux_names_;
  Teuchos::RCP<Epetra_MultiVector> aux_output_;

  // Auxiliary data, maintained by Amanzi and updated
  Teuchos::RCP<Epetra_MultiVector> aux_data_;

  void UpdateChemistryStateStorage(void);
  int InitializeSingleCell(int cellIndex, const std::string& condition);
  int AdvanceSingleCell(double delta_time, 
                        Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star,
                        int cellIndex);

  void ParseChemicalConditions(const Teuchos::ParameterList& param_list,
                               std::map<std::string, std::string>& conditions);
  void XMLParameters(void);

  // These helpers copy data back and forth between a set of buffers and the chemistry state.
  // given cell.
  void CopyAmanziStateToAlquimia(const int cell_id,
                                 Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                                 AlquimiaState* state,
                                 AlquimiaMaterialProperties* mat_props,
                                 AlquimiaAuxiliaryData* aux_data);

  void CopyAlquimiaStateToAmanzi(const int cell_id,
                                 const AlquimiaState* state,
                                 const AlquimiaMaterialProperties* mat_props,
                                 const AlquimiaAuxiliaryData* aux_data,
                                 const AlquimiaAuxiliaryOutputData* aux_output);

};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif  // AMANZI_ALQUIMIA_CHEMISTRY_PK_HH_
