/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/
 
#ifndef CHEMISTRY_AMANZI_PK_HH_
#define CHEMISTRY_AMANZI_PK_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Beaker.hh"
#include "BeakerState.hh"
#include "Chemistry_PK.hh"
#include "PK_Factory.hh"
#include "Mesh.hh"
#include "TreeVector.hh"

namespace Amanzi {
namespace AmanziChemistry {

// Trilinos based chemistry process kernel for the unstructured mesh
class Amanzi_PK : public Chemistry_PK {
 public:
  Amanzi_PK(Teuchos::ParameterList& pk_tree,
            const Teuchos::RCP<Teuchos::ParameterList>& glist,
            const Teuchos::RCP<State>& S,
            const Teuchos::RCP<TreeVector>& soln);

  ~Amanzi_PK();

  // members required by PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual void set_dt(double dt) {};
  virtual double get_dt() { return dt_next_; }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) { extra_chemistry_output_data(); }

  virtual std::string name() { return "chemistry amanzi"; }

  // The following two routines provide the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> extra_chemistry_output_data();
  void set_chemistry_output_names(std::vector<std::string>* names);

  // functions used in Rransport PK
  void CopyCellStateToBeakerState(
      int c, Teuchos::RCP<Epetra_MultiVector> aqueous_components);

  // access
  Beaker* get_engine() { return chem_; }
  const BeakerParameters& beaker_parameters() const { return beaker_parameters_; }
  BeakerState beaker_state() { return beaker_state_; }

 private:
  void AllocateAdditionalChemistryStorage_();

  void XMLParameters();
  void SetupAuxiliaryOutput();
  void SizeBeakerState_();

  void CopyBeakerStructuresToCellState(
      int c, Teuchos::RCP<Epetra_MultiVector> aqueous_components);

 protected:
  Teuchos::RCP<TreeVector> soln_;

 private:
  Beaker* chem_;
  BeakerParameters beaker_parameters_;
  BeakerState beaker_state_, beaker_state_copy_;

  std::string dt_control_method_;
  double current_time_, saved_time_;
  double dt_max_, dt_next_, dt_cut_factor_, dt_increase_factor_;
  int dt_cut_threshold_, dt_increase_threshold_;

  std::vector<std::string> aux_names_;
  std::vector<int> aux_index_;
  Teuchos::RCP<Epetra_MultiVector> aux_data_;

 private:
  // factory registration
  static RegisteredPKFactory<Amanzi_PK> reg_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
