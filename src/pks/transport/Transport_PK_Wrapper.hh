/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Ethan Coon

  Temporary wrapper converting the Transport_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.

*/

#ifndef AMANZI_TRANSPORT_PK_WRAPPER_HH_
#define AMANZI_TRANSPORT_PK_WRAPPER_HH_

#include "Teuchos_RCP.hpp"

#include "FnTimeIntegratorPK.hh"
#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

class Transport_PK_Wrapper : public FnTimeIntegratorPK {

 public:
  Transport_PK_Wrapper(const Teuchos::RCP<Teuchos::ParameterList>& plist,
		  Teuchos::ParameterList& FElist,
		  const Teuchos::RCP<TreeVector>& soln) {
    pk_ = Teuchos::rcp(new Transport_PK(plist, FElist, soln));
  }

  // Setup
  virtual void Setup(const Teuchos::Ptr<State>& S) {
    pk_->Setup(S);
  }

  // Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S) {
    pk_->Initialize(S);
  }

  // Choose a time step compatible with physics.
  virtual double get_dt() {
    return pk_->get_dt();
  }

  // Advance from state S0 to state S1 at time S0.time + dt.
  // Due to Transport PK / MPC conflict (FIXME when MPC will be upgraded)
  //  virtual int Advance(double dt, double& dt_actual) = 0;
  virtual bool Advance(double dt) {
    return pk_->Advance(dt);
  }

  // Commit any secondary (dependent) variables.
  virtual void CommitState(double dt, const Teuchos::Ptr<State>& S) {
    pk_->CommitState(dt, S);
  }

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::Ptr<State>& S) {
    pk_->CalculateDiagnostics(S);
  }

  virtual void SetState(const Teuchos::RCP<State>& S) {
    pk_->SetState(S);
  }

  virtual std::string name() {
    return pk_->name();
  }

  // Time integration interface
  // computes the non-linear functional f = f(t,u,udot)
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {
    pk_->Functional(t_old, t_new, u_old->Data(),
		    u_new->Data(), f->Data());    
  }

  // applies preconditioner to u and returns the result in Pu
  virtual void ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
    pk_->ApplyPreconditioner(u->Data(), Pu->Data());
  }

  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u, 
			   Teuchos::RCP<const TreeVector> du) {
    return pk_->ErrorNorm(u->Data(), du->Data());
  }

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, 
				    double h) {
    pk_->UpdatePreconditioner(t, up->Data(), h);
  }

  // check the admissibility of a solution
  // override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) {
    return pk_->IsAdmissible(up->Data());
  }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, 
			       Teuchos::RCP<TreeVector> u) {
    return pk_->ModifyPredictor(h, u0->Data(), u->Data());
  }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, 
		       Teuchos::RCP<TreeVector> du) {
    return pk_->ModifyCorrection(h, res->Data(), u->Data(), du->Data());
  }

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() {
    pk_->ChangedSolution();
  }

 protected:
  Teuchos::RCP<Transport_PK> pk_;

};

} // namespace
} // namespace

#endif
