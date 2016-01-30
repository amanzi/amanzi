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
#include "beaker.hh"
#include "chemistry_exception.hh"
#include "chemistry_verbosity.hh"
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
  virtual void Setup();
  virtual void Initialize();

  virtual void set_dt(double dt) {};
  virtual double get_dt() { return this->max_time_step_; }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);
  virtual void CommitStep(double t_old, double t_new);
  virtual void CalculateDiagnostics() { extra_chemistry_output_data(); }

  virtual std::string name() { return "chemistry amanzi"; }

  // The following two routines provide the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> extra_chemistry_output_data();
  void set_chemistry_output_names(std::vector<std::string>* names);

 private:
  void AllocateAdditionalChemistryStorage_(const Beaker::BeakerComponents& components);

  void XMLParameters();
  void SetupAuxiliaryOutput();
  void SizeBeakerStructures_();

  void CopyCellStateToBeakerStructures(
      int cell_id, Teuchos::RCP<Epetra_MultiVector> aqueous_components);
  void CopyBeakerStructuresToCellState(
      int cell_id, Teuchos::RCP<Epetra_MultiVector> aqueous_components);

 protected:
  Teuchos::RCP<TreeVector> soln_;

 private:
  Teuchos::RCP<Teuchos::ParameterList> cp_list_;

  Beaker* chem_;
  Beaker::BeakerParameters beaker_parameters_;
  Beaker::BeakerComponents beaker_components_;
  Beaker::BeakerComponents beaker_components_copy_;

  double max_time_step_;
  double current_time_;
  double saved_time_;

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
