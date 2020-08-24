/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Trilinos based chemistry process kernel for the unstructured mesh.
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

// Amanzi
#include "ChemistryEngine.hh"
#include "PK_Factory.hh"
#include "TreeVector.hh"

// Chemistry PK
#include "Chemistry_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

#ifdef ALQUIMIA_ENABLED
class Alquimia_PK: public Chemistry_PK {
 public:
  Alquimia_PK(Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& soln);

  ~Alquimia_PK();

  // members required by PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual void set_dt(double dt) {};
  virtual double get_dt() { return this->time_step_; }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) { extra_chemistry_output_data(); }

  virtual std::string name() { return "chemistry alquimia"; }

  void CopyFieldstoNewState(const Teuchos::RCP<State>& S_next);

  // Ben: the following routine provides the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> extra_chemistry_output_data();

  // Copies the chemistry state in the given cell to the given Alquimia containers.
  void CopyToAlquimia(int cell_id,
                      AlquimiaProperties& mat_props,
                      AlquimiaState& state,
                      AlquimiaAuxiliaryData& aux_data);
  
 private:
  // Copy cell state to the given Alquimia containers taking 
  // the aqueous components from the given multivector.
  void CopyToAlquimia(int cell_id,
                      Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                      AlquimiaProperties& mat_props,
                      AlquimiaState& state,
                      AlquimiaAuxiliaryData& aux_data);

  // Copy the data in the given Alquimia containers to the given cell state.
  // The aqueous components are placed into the given multivector.
  void CopyFromAlquimia(const int cell,
                        const AlquimiaProperties& mat_props,
                        const AlquimiaState& state,
                        const AlquimiaAuxiliaryData& aux_data,
                        const AlquimiaAuxiliaryOutputData& aux_output,
                        Teuchos::RCP<Epetra_MultiVector> aqueous_components);

  int InitializeSingleCell(int cell, const std::string& condition);
  int AdvanceSingleCell(double dt, Teuchos::RCP<Epetra_MultiVector>& aquesous_components,
                        int cell);

  void ParseChemicalConditionRegions(const Teuchos::ParameterList& param_list,
                                     std::map<std::string, std::string>& conditions);
  void XMLParameters();

  void CopyAlquimiaStateToAmanzi(const int cell,
                                 const AlquimiaProperties& mat_props,
                                 const AlquimiaState& state,
                                 const AlquimiaAuxiliaryData& aux_data,
                                 const AlquimiaAuxiliaryOutputData& aux_output,
                                 Teuchos::RCP<Epetra_MultiVector> aquesous_components);

  void ComputeNextTimeStep();

  // maps
  void InitializeAuxNamesMap_();

 protected:
  Teuchos::RCP<TreeVector> soln_;

 private:
  // Time stepping controls. Some parameters are defined in the base class
  double time_step_, max_time_step_, min_time_step_, prev_time_step_;
  std::string time_step_control_method_;
  int num_iterations_for_time_step_cut_, num_steps_before_time_step_increase_;
  double time_step_cut_factor_, time_step_increase_factor_;

  bool chem_initialized_;

  // Alquimia data structures for interface with Amanzi.
  AlquimiaState alq_state_;
  AlquimiaProperties alq_mat_props_;
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
  Teuchos::RCP<Epetra_MultiVector> aux_data_;

  std::vector<std::vector<int> > map_;
  std::vector<std::string> mineral_names_, primary_names_;

 private:
  // factory registration
  static RegisteredPKFactory<Alquimia_PK> reg_;
};
#endif

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif

