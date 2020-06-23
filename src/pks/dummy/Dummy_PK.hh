/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Daniil Svyatskiy

  Dummy PK which demonstrates the require interface for PK
  BDFFnBase<CompositeVector>, to use TreeVectors.
*/

#ifndef AMANZI_DUMMY_PK_HH_
#define AMANZI_DUMMY_PK_HH_

#include "Teuchos_RCP.hpp"

#include "TreeVector.hh"
#include "PK_PhysicalBDF.hh"
#include "PK_Factory.hh"

namespace Amanzi {

class Dummy_PK : public PK_PhysicalBDF {
 public:
  Dummy_PK(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln);

  // Setup
  virtual void Setup(const Teuchos::Ptr<State>& S) {dummy_dt = 1000; step_count = 1;}

  // Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S) {}
  
  // Choose a time step compatible with physics.
  virtual double get_dt() {
    return dummy_dt;
  }

  virtual void set_dt(double dt) {
    dummy_dt = dt;
  };

  // Advance PK by step size dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false);

  // Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {};

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {};

  virtual std::string name() { return "dummy_pk"; }

  // Time integration interface
  // computes the non-linear functional f = f(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {};

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) { return 1; }

  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) { return 0.0; }

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {};

  // check the admissibility of a solution
  // override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) { return true; }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) { return true; }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> du) {
      return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED; }

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() {};

 protected:
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::ParameterList ti_list_;
  Teuchos::RCP<TreeVector> soln_;

  double dummy_dt;
  int step_count;

 private:
  // factory registration
  static RegisteredPKFactory<Dummy_PK> reg_;
};

}  // namespace Amanzi

#endif
