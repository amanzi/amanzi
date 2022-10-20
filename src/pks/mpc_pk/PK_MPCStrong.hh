/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
          Daniil Svyatskiy
          Konstantin Lipnikov

  Interface for the derived PK_MPCStrong class.  Is both a PK and a Model
  Evalulator, providing needed methods for BDF time integration of the
  coupled system.

  Completely automated and generic to any sub PKs, this uses a block diagonal
  preconditioner.

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#ifndef AMANZI_PK_MPC_STRONG_HH_
#define AMANZI_PK_MPC_STRONG_HH_

#include <vector>
#include "Teuchos_ParameterList.hpp"

#include "BDF1_TI.hh"
#include "PK_BDF.hh"

#include "PK_MPC.hh"
#include "State.hh"
#include "PK_Factory.hh"
#include "TreeOperator.hh"


namespace Amanzi {

template<class PK_Base>
class PK_MPCStrong : virtual public PK_MPC<PK_Base>, public PK_BDF
{
 public:
  PK_MPCStrong(Teuchos::ParameterList& pk_tree,
            const Teuchos::RCP<Teuchos::ParameterList>& global_list,
            const Teuchos::RCP<State>& S,
            const Teuchos::RCP<TreeVector>& soln);

  // PK_MPCStrong is a PK
  virtual void Setup();
  virtual void Initialize();

  // -- dt is the minimum of the sub pks
  virtual double get_dt() { return dt_; }
  virtual void set_dt(double dt) { dt_ = dt; }

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  // PK_MPCStrong is an PK_Implicit
  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk FunctionalResidual().
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // -- enorm for the coupled system
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du,
                           Teuchos::RCP<const TreeVector> res,
                           const AmanziSolvers::ConvergenceMonitor& monitor);

  // PK_MPCStrong's preconditioner is, by default, just the block-diagonal
  // operator formed by placing the sub PK's preconditioners on the diagonal.
  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void ChangedSolution();

  // -- Admissibility of the solution.
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> u);

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> du);

  // access
  Teuchos::RCP<Operators::TreeOperator> op_tree_matrix() { return op_tree_matrix_; }
  Teuchos::RCP<Operators::TreeOperator> op_tree_pc() { return op_tree_pc_; }
  Teuchos::RCP<TreeVector> op_tree_rhs() { return op_tree_rhs_; }
  
 protected:
  using PK_MPC<PK_Base>::sub_pks_;
  using PK_MPC<PK_Base>::S_;

  using PK_MPC<PK_Base>::global_list_;
  using PK_MPC<PK_Base>::pk_tree_;
  using PK_MPC<PK_Base>::my_list_;

  using PK_MPC<PK_Base>::solution_;

  // timestep control
  double dt_;
  Teuchos::RCP<Amanzi::BDF1_TI<TreeVector, TreeVectorSpace> > time_stepper_;

  Teuchos::RCP<Operators::TreeOperator> op_tree_matrix_, op_tree_pc_;
  Teuchos::RCP<TreeVector> op_tree_rhs_;
  
 private:
  // factory registration
  static RegisteredPKFactory<PK_MPCStrong<PK_Base> > reg_;
};


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template<class PK_Base>
PK_MPCStrong<PK_Base>::PK_MPCStrong(Teuchos::ParameterList& pk_tree,
                              const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                              const Teuchos::RCP<State>& S,
                              const Teuchos::RCP<TreeVector>& soln) :
  PK_MPC<PK_Base>(pk_tree, global_list, S, soln), dt_(0.0) {};


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
template<class PK_Base>
void PK_MPCStrong<PK_Base>::Setup()
{
  // Tweak the sub-PK parameter lists. This allows the PK to
  // potentially not assemble things.
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(global_list_, "PKs");

  for (auto param = pk_tree_.begin(); param != pk_tree_.end(); ++param) {
    std::string pname = param->first;
    if (pks_list->isSublist(pname)) {
      pks_list->sublist(pname).set("strongly coupled PK", true);
    }
  }

  // call each sub-PKs Setup()
  PK_MPC<PK_Base>::Setup();

  // Set the initial timestep as the min of the sub-pk sizes.
  dt_ = get_dt();
}


// -----------------------------------------------------------------------------
// Initialize each sub-PK and the time integrator.
// -----------------------------------------------------------------------------
template<class PK_Base>
void PK_MPCStrong<PK_Base>::Initialize()
{
  // Just calls both subclass's initialize.  NOTE - order is important
  // here -- MPC<PK_Base> grabs the primary variables from each sub-PK
  // and stuffs them into the solution, which must be done prior to
  // initializing the timestepper.

  // Initialize all sub PKs.
  PK_MPC<PK_Base>::Initialize();

  // set up the timestepping algorithm if this is not strongly coupled
  if (!my_list_->template get<bool>("strongly coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::ParameterList& ts_plist = my_list_->sublist("time integrator").sublist("BDF1");
    ts_plist.set("initial time", S_->get_time());
    time_stepper_ = Teuchos::rcp(new Amanzi::BDF1_TI<TreeVector,
        TreeVectorSpace>(*this, ts_plist, solution_));

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->SetInitialState(S_->get_time(), solution_, solution_dot);
  }
}


// -----------------------------------------------------------------------------
// Make one time step 
// -----------------------------------------------------------------------------
template<class PK_Base>
bool PK_MPCStrong<PK_Base>::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  // save a copy of solution, i.e. primary variables
  TreeVector solution_copy(*solution_);

  // take a bdf timestep
  double dt_solver;
  bool fail;
  if (true) { // this is here simply to create a context for timer,
              // which stops the clock when it is destroyed at the
              // closing brace.
    fail = time_stepper_->TimeStep(dt_, dt_solver, solution_);
  }

  if (!fail) {
    // commit the step as successful
    time_stepper_->CommitSolution(dt_, solution_);
    dt_ = dt_solver;
  } else {
    // take the decreased timestep size
    dt_ = dt_solver;

    // recover the original solution
    *solution_ = solution_copy;
    ChangedSolution();
  }

  return fail;
}


// -----------------------------------------------------------------------------
// Compute the non-linear functional g = g(t,u,udot).
// -----------------------------------------------------------------------------
template<class PK_Base>
void PK_MPCStrong<PK_Base>::FunctionalResidual(
    double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // loop over sub-PKs
  for (unsigned int i = 0; i != sub_pks_.size(); ++i) {
    // pull out the old solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_old(Teuchos::null);
    if (u_old != Teuchos::null) {
      pk_u_old = u_old->SubVector(i);
      if (pk_u_old == Teuchos::null) {
        Errors::Message message("MPC: vector structure does not match PK structure");
        Exceptions::amanzi_throw(message);
      }
    }

    // pull out the new solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_new = u_new->SubVector(i);
    if (pk_u_new == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the residual sub-vector
    Teuchos::RCP<TreeVector> pk_g = g->SubVector(i);
    if (pk_g == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // fill the nonlinear function with each sub-PKs contribution
    sub_pks_[i]->FunctionalResidual(t_old, t_new, pk_u_old, pk_u_new, pk_g);
  }
}


// -----------------------------------------------------------------------------
// Applies preconditioner to u and returns the result in Pu.
// -----------------------------------------------------------------------------
template<class PK_Base>
int PK_MPCStrong<PK_Base>::ApplyPreconditioner(
    Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  int ierr(0);
  // loop over sub-PKs
  for (unsigned int i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the preconditioned u sub-vector
    Teuchos::RCP<TreeVector> pk_Pu = Pu->SubVector(i);
    if (pk_Pu == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // Fill the preconditioned u as the block-diagonal product using each sub-PK.
    sub_pks_[i]->ApplyPreconditioner(pk_u, pk_Pu);
  }
  return ierr;
}


// -----------------------------------------------------------------------------
// Compute a norm on u-du-res and returns the result.
// For a Strong MPC, the enorm is just the max of the sub PKs enorms.
// -----------------------------------------------------------------------------
template<class PK_Base>
double PK_MPCStrong<PK_Base>::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                        Teuchos::RCP<const TreeVector> du)
{
  double norm = 0.0;

  // loop over sub-PKs
  for (unsigned int i = 0; i != sub_pks_.size(); ++i) {
    // pull out sub-vectors
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    Teuchos::RCP<const TreeVector> pk_du = du->SubVector(i);

    if (pk_u == Teuchos::null || pk_du == Teuchos::null) {
      Errors::Message message("PK_MPCStrong: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // norm is the max of the sub-PK norms
    double tmp_norm = sub_pks_[i]->ErrorNorm(pk_u, pk_du);
    norm = std::max(norm, tmp_norm);
  }
  return norm;
}


// -----------------------------------------------------------------------------
// Compute a norm on u-du-res and returns the result.
// For a Strong MPC, the enorm is just the max of the sub PKs enorms.
// -----------------------------------------------------------------------------
template<class PK_Base>
double PK_MPCStrong<PK_Base>::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                        Teuchos::RCP<const TreeVector> du,
                                        Teuchos::RCP<const TreeVector> res,
                                        const AmanziSolvers::ConvergenceMonitor& monitor)
{
  double norm = 0.0;

  // loop over sub-PKs
  for (unsigned int i = 0; i != sub_pks_.size(); ++i) {
    // pull out sub-vectors
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    Teuchos::RCP<const TreeVector> pk_du = du->SubVector(i);
    Teuchos::RCP<const TreeVector> pk_res = res->SubVector(i);

    if (pk_u == Teuchos::null || pk_du == Teuchos::null || pk_res == Teuchos::null) {
      Errors::Message message("PK_MPCStrong: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // norm is the max of the sub-PK norms
    double tmp_norm = sub_pks_[i]->ErrorNorm(pk_u, pk_du, pk_res, monitor);
    norm = std::max(norm, tmp_norm);
  }
  return norm;
}


// -----------------------------------------------------------------------------
// Update the preconditioner.
// -----------------------------------------------------------------------------
template<class PK_Base>
void PK_MPCStrong<PK_Base>::UpdatePreconditioner(
    double t, Teuchos::RCP<const TreeVector> up, double h)
{
  for (unsigned int i = 0; i != sub_pks_.size(); ++i) {
    // pull out the up sub-vector
    Teuchos::RCP<const TreeVector> pk_up = up->SubVector(i);
    if (pk_up == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // update precons of each of the sub-PKs
    sub_pks_[i]->UpdatePreconditioner(t, pk_up, h);
  }
}


// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time integration
// scheme is changing the value of the solution in state.
// -----------------------------------------------------------------------------
template<class PK_Base>
void PK_MPCStrong<PK_Base>::ChangedSolution() {
  // loop over sub-PKs
  for (typename PK_MPC<PK_Base>::SubPKList::iterator pk = PK_MPC<PK_Base>::sub_pks_.begin();
      pk != PK_MPC<PK_Base>::sub_pks_.end(); ++pk) {
    (*pk)->ChangedSolution();
  }
}


// -----------------------------------------------------------------------------
// Check admissibility of each sub-pk
// -----------------------------------------------------------------------------
template<class PK_Base>
bool PK_MPCStrong<PK_Base>::IsAdmissible(Teuchos::RCP<const TreeVector> u)
{
  // First ensure each PK thinks we are admissible -- this will ensure
  // the residual can at least be evaluated.
  for (unsigned int i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    if (!sub_pks_[i]->IsAdmissible(pk_u)) {
      return false;
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
// Modify predictor from each sub pk.
// -----------------------------------------------------------------------------
template<class PK_Base>
bool PK_MPCStrong<PK_Base>::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                                         Teuchos::RCP<TreeVector> u)
{
  // loop over sub-PKs
  bool modified = false;
  for (unsigned int i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u0 = u0->SubVector(i);
    Teuchos::RCP<TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    modified |= sub_pks_[i]->ModifyPredictor(h, pk_u0, pk_u);
  }
  return modified;
}


// -----------------------------------------------------------------------------
// Modify correction from each sub pk.
// -----------------------------------------------------------------------------
template<class PK_Base>
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
    PK_MPCStrong<PK_Base>::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                                         Teuchos::RCP<const TreeVector> u,
                                         Teuchos::RCP<TreeVector> du)
{
  // loop over sub-PKs
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      modified = AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  for (unsigned int i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    Teuchos::RCP<const TreeVector> pk_res = res->SubVector(i);
    Teuchos::RCP<TreeVector> pk_du = du->SubVector(i);

    if (pk_u == Teuchos::null || pk_du == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    modified = std::max(modified, sub_pks_[i]->ModifyCorrection(h, pk_res, pk_u, pk_du));
  }
  return modified;
}

}  // namespace Amanzi

#endif
