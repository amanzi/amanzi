/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors:: Daniil Svyatskiy

  Temporary wrapper converting the Chemistry_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.

*/

#ifndef AMANZI_CHEMISTRY_PK_WRAPPER_HH_
#define AMANZI_CHEMISTRY_PK_WRAPPER_HH_

#include "Teuchos_RCP.hpp"

#include "../pks/PK.hh"
#include "chemistry_pk.hh"
#include "../pks/PK_Factory.hh"
#include "Chemistry_State.hh"

namespace Amanzi {
namespace AmanziChemistry{

class Chemistry_PK_Wrapper : public PK {
 public:
  Chemistry_PK_Wrapper(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln);
  // Setup
  virtual void Setup() {};

  // Initialize owned (dependent) variables.
  virtual void Initialize() {  
    dt_ = -1;
    CS->Initialize();
    pk_->InitializeChemistry();
  }

  // Choose a time step compatible with physics.
  virtual double get_dt() { return pk_->time_step(); }
  virtual void set_dt(double dt) { dt_ = dt; }

  // Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool AdvanceStep(double t_old, double t_new){
    pk_->Advance(t_new - t_old, total_component_concentration);
    return true;
  }

  // Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new) {
    pk_->CommitState(CS, t_new);
  }

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics() {
    // get the auxillary data
    Teuchos::RCP<Epetra_MultiVector> aux = pk_->get_extra_chemistry_output_data();
  }

  virtual std::string name() { return pk_->name(); }

  Teuchos::RCP<Epetra_MultiVector> get_total_component_concentration() { return total_component_concentration; }

  void set_total_component_concentration(Teuchos::RCP<Epetra_MultiVector> tcc) {
    total_component_concentration = tcc;
  }

 protected:
  // storage for the component concentration intermediate values
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration;
  Teuchos::RCP<AmanziChemistry::Chemistry_State> CS;
  std::vector<std::string> comp_names_;
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Chemistry_PK> pk_;
  Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<State> S_;
  double dt_;

 private:
  // factory registration
  static RegisteredPKFactory<Chemistry_PK_Wrapper> reg_;
};

}  // namespace Chemistry
}  // namespace Amanzi

#endif
