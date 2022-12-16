/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry PK

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
class Alquimia_PK : public Chemistry_PK {
 public:
  Alquimia_PK(Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& soln);

  ~Alquimia_PK();

  // members required by PK interface
  virtual void Setup() final;
  virtual void Initialize() final;

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) final;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;
  virtual void CalculateDiagnostics(const Tag& tag) final { extra_chemistry_output_data(); }

  // Ben: the following routine provides the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> extra_chemistry_output_data();

  // Copies the chemistry state in the given cell to the given Alquimia containers.
  void CopyToAlquimia(int cell_id,
                      AlquimiaProperties& mat_props,
                      AlquimiaState& state,
                      AlquimiaAuxiliaryData& aux_data,
                      const Tag& water_tag = Tags::DEFAULT);

 private:
  // Copy cell state to the given Alquimia containers taking
  // the aqueous components from the given multivector.
  void CopyToAlquimia(int cell_id,
                      Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                      AlquimiaProperties& mat_props,
                      AlquimiaState& state,
                      AlquimiaAuxiliaryData& aux_data,
                      const Tag& water_tag = Tags::DEFAULT);

  // Copy the data in the given Alquimia containers to the given cell state.
  // The aqueous components are placed into the given multivector.
  void CopyFromAlquimia(const int cell,
                        const AlquimiaProperties& mat_props,
                        const AlquimiaState& state,
                        const AlquimiaAuxiliaryData& aux_data,
                        const AlquimiaAuxiliaryOutputData& aux_output,
                        Teuchos::RCP<Epetra_MultiVector> aqueous_components);

  int InitializeSingleCell(int cell, const std::string& condition);
  int AdvanceSingleCell(double dt, Teuchos::RCP<Epetra_MultiVector>& aquesous_components, int cell);

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

 private:
  // Time stepping controls. Some parameters are defined in the base class
  std::string dt_control_method_;

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
  std::vector<std::vector<std::string>> aux_subfield_names_;
  Teuchos::RCP<Epetra_MultiVector> aux_output_;
  Teuchos::RCP<Epetra_MultiVector> aux_data_;

  std::vector<std::vector<int>> map_;
  std::vector<std::string> mineral_names_, primary_names_;

 private:
  // factory registration
  static RegisteredPKFactory<Alquimia_PK> reg_;
};
#endif

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
