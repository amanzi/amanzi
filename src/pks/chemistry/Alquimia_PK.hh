/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/
 
#ifndef AMANZI_CHEMISTRY_ALQUIMIA_PK_HH_
#define AMANZI_CHEMISTRY_ALQUIMIA_PK_HH_

#include <map>
#include <string>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Chemistry
#include "Chemistry_PK.hh"
#include "ChemistryEngine.hh"

namespace Amanzi {
namespace AmanziChemistry {

// Trilinos based chemistry process kernel for the unstructured mesh
class Alquimia_PK: public Chemistry_PK {
 public:
  // Constructor. Note that we must pass the "Main" parameter list
  // to this PK so that it has access to all information about the 
  // problem.
  Alquimia_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
              Teuchos::RCP<ChemistryEngine> chem_engine,
              Teuchos::RCP<State> S,
              Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  ~Alquimia_PK();

  // members required by PK interface
  virtual void Setup();
  virtual void Initialize();

  void InitializeChemistry(void);

  void Advance(const double& delta_time,
               Teuchos::RCP<Epetra_MultiVector> total_component_concentration);
  void CommitState(const double& time);

  double time_step(void) const { return this->time_step_; }

  // Ben: the following routine provides the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> get_extra_chemistry_output_data();

  // Copies the chemistry state in the given cell to the given Alquimia containers.
  void CopyToAlquimia(const int cell_id,
                      AlquimiaMaterialProperties& mat_props,
                      AlquimiaState& state,
                      AlquimiaAuxiliaryData& aux_data);
  
 private:
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

  void UpdateChemistryStateStorage(void);
  int InitializeSingleCell(int cell_index, const std::string& condition);
  int AdvanceSingleCell(double delta_time, 
                        Teuchos::RCP<Epetra_MultiVector> total_component_concentration,
                        int cell_index);

  void ParseChemicalConditionRegions(const Teuchos::ParameterList& param_list,
                                     std::map<std::string, std::string>& conditions);
  void XMLParameters(void);

  void CopyAlquimiaStateToAmanzi(const int cell_id,
                                 const AlquimiaMaterialProperties& mat_props,
                                 const AlquimiaState& state,
                                 const AlquimiaAuxiliaryData& aux_data,
                                 const AlquimiaAuxiliaryOutputData& aux_output,
                                 Teuchos::RCP<Epetra_MultiVector> total_component_concentration);

 private:
  Teuchos::RCP<Teuchos::ParameterList> glist_, cp_list_;

  // Time stepping controls.
  double time_step_, max_time_step_, min_time_step_, prev_time_step_;
  std::string time_step_control_method_;
  int num_iterations_for_time_step_cut_, num_steps_before_time_step_increase_;
  double time_step_cut_factor_, time_step_increase_factor_;
  int num_iterations_, num_successful_steps_;
  void ComputeNextTimeStep();

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

  int num_aux_data_;
  Teuchos::RCP<Epetra_MultiVector> aux_data_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif

