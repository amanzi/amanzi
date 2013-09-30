/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_ALQUIMIA_CHEMISTRY_PK_HH_
#define AMANZI_ALQUIMIA_CHEMISTRY_PK_HH_

#include <string>
#include <vector>
#include <map>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "chemistry_pk_base.hh"

#include "alquimia_memory.h"
#include "alquimia_util.h"
#include "alquimia_constants.h"
#include "alquimia_containers.h"
#include "alquimia_interface.h"

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
                        Teuchos::RCP<Chemistry_State> chem_state);

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

  // Alquimia data structures.
  bool chem_initialized_;
  AlquimiaInterface chem_;
  AlquimiaEngineStatus chem_status_;
  AlquimiaData chem_data_;

  // Mapping of region names to geochemical conditions. A region is identified 
  // by a string, and all cells within a region will have all geochemical 
  // conditions in the corresponding condition vector applied to them. NOTE
  // that these maps do not own the geochemical conditions--they only hold 
  // pointers to the objects.
  std::map<std::string, AlquimiaGeochemicalCondition*> chem_initial_conditions_;
  std::map<std::string, AlquimiaGeochemicalCondition*> chem_boundary_conditions_;
  
  // Vector that takes responsibility for ownership of geochemical conditions.
  std::vector<AlquimiaGeochemicalCondition*> all_chem_conditions_;

  // Back-end engine name and input file.
  std::string chem_engine_inputfile_;
  std::string chem_engine_name_;

  double current_time_;
  double saved_time_;

  // Auxiliary data, requested by and stored within Amanzi.
  std::vector<std::string> aux_names_;
  Teuchos::RCP<Epetra_MultiVector> aux_data_;

  void UpdateChemistryStateStorage(void);
  int InitializeSingleCell(int cellIndex, AlquimiaGeochemicalCondition* condition);
  int AdvanceSingleCell(double delta_time, 
                        Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star,
                        int cellIndex, AlquimiaGeochemicalCondition* condition);

  void ParseChemicalConditions(const Teuchos::ParameterList& param_list,
                               std::map<std::string, AlquimiaGeochemicalCondition*>& conditions);
  void XMLParameters(void);
  void SetupAuxiliaryOutput(void);

  // These helpers copy data at the given cell from Amanzi's chemistry state to 
  // their corresponding locations within Alquimia.
  void CopyAmanziStateToAlquimia(
      const int cell_id,
      Teuchos::RCP<const Epetra_MultiVector> aqueous_components);
  void CopyAmanziMaterialPropertiesToAlquimia(
      const int cell_id,
      Teuchos::RCP<const Epetra_MultiVector> aqueous_components);
  void CopyAmanziGeochemicalConditionsToAlquimia(
      const int cell_id,
      Teuchos::RCP<const Epetra_MultiVector> aqueous_components);

  // These helpers copy Alquimia's data to Amanzi's chemistry state at the 
  // given cell.
  void CopyAlquimiaStateToAmanzi(const int cell_id);
  void CopyAlquimiaMaterialPropertiesToAmanzi(const int cell_id);

};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif  // AMANZI_ALQUIMIA_CHEMISTRY_PK_HH_
