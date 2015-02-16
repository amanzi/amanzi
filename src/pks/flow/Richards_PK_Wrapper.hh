/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon

  Temporary wrapper converting the Richards_PK, which inherits from
  BDFFnBase<CompositeVector>, to use TreeVectors.
*/

#ifndef AMANZI_RICHARDS_PK_WRAPPER_HH_
#define AMANZI_RICHARDS_PK_WRAPPER_HH_

#include "Teuchos_RCP.hpp"

#include "TreeVector.hh"
#include "FnTimeIntegratorPK.hh"
#include "Richards_PK.hh"
#include "PK_Factory.hh"

namespace Amanzi {

class State;

namespace Flow {

class Richards_PK_Wrapper : public FnTimeIntegratorPK {
 public:
  Richards_PK_Wrapper(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln);

  ~Richards_PK_Wrapper() {};

  // Setup
  virtual void Setup() {
    dt_ = -1.0;
    pk_->InitializeFields();
  }

  // Initialize owned (dependent) variables.
  virtual void Initialize() {
    pk_->Initialize(S_.ptr());
    pk_->InitializeAuxiliaryData(); 
    pk_->InitTimeInterval();
    pk_->UpdateAuxilliaryData();
  }

  // Choose a time step compatible with physics.
  virtual double get_dt() {    
    return pk_->get_dt();    
  }

  //  Set a time step
  virtual void set_dt(double dt){
    dt_ = dt;
    pk_->set_dt(dt);
  }

  // Advance PK by step size dt.
  virtual bool AdvanceStep(double t_old, double t_new);

  // Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new) {
    pk_->CommitState(t_new-t_old, S_.ptr());
  }

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics() {
    pk_->CalculateDiagnostics(S_.ptr());
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
  // Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::ParameterList ti_list_;
  Teuchos::RCP<Richards_PK> pk_; 
  Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<State> S_;
  double dt_;

 private:
  // factory registration
  static RegisteredPKFactory<Richards_PK_Wrapper> reg_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
