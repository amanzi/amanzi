/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived StrongMPC class.  Is both a PK and a Model
Evalulator, providing needed methods for BDF time integration of the coupled
system.

Completely automated and generic to any sub PKs, this uses a block diagonal
preconditioner.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "errors.hh"
#include "bdf1_time_integrator.hh"

#include "strong_mpc.hh"

namespace Amanzi {

RegisteredPKFactory<StrongMPC> StrongMPC::reg_("strong MPC");

StrongMPC::StrongMPC(Teuchos::ParameterList &mpc_plist, const Teuchos::RCP<State>& S,
         const Teuchos::RCP<TreeVector>& solution) :
    MPC(mpc_plist, S, solution) {

  solution_ = solution;

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  Teuchos::ParameterList pks_list = mpc_plist.sublist("PKs");
  PKFactory pk_factory;

  if (mpc_plist.isParameter("PKs order")) {
    // ordered
    Teuchos::Array<std::string> pk_order = mpc_plist.get< Teuchos::Array<std::string> >("PKs order");
    int npks = pk_order.size();

    for (int i=0; i!=npks; ++i) {
      std::string name_i = pk_order[i];

      // set up the solution
      Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(name_i));
      solution_->PushBack(pk_soln);

      // get the list and set a flag so the PK knows it is strongly coupled
      Teuchos::ParameterList pk_list = pks_list.sublist(name_i);
      pk_list.set<bool>("strongly coupled PK", true);

      // instantiate the PK
      Teuchos::RCP<PK> pk = pk_factory.CreatePK(pk_list, S, pk_soln);
      pk->set_name(name_i);
      sub_pks_.push_back(pk);
    }

  } else {
    // no order, just loop over all sublists
    for (Teuchos::ParameterList::ConstIterator i = pks_list.begin();
         i != pks_list.end(); ++i) {

      const std::string &name_i  = pks_list.name(i);
      const Teuchos::ParameterEntry  &entry_i = pks_list.entry(i);
      if (entry_i.isList()) {

        // set up the solution
        Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(name_i));
        solution_->PushBack(pk_soln);

        // get the list and set a flag so the PK knows it is strongly coupled
        Teuchos::ParameterList pk_list = pks_list.sublist(name_i);
        pk_list.set<bool>("strongly coupled PK", true);

        // instantiate the PK
        Teuchos::RCP<PK> pk = pk_factory.CreatePK(pk_list, S, pk_soln);
        pk->set_name(name_i);
        sub_pks_.push_back(pk);
      }
    }
  }
};


void StrongMPC::initialize(const Teuchos::RCP<State>& S) {
  // Initialize all sub PKs.
  MPC::initialize(S);

  // Initialize the time integration scheme for the coupled problem.
  atol_ = mpc_plist_.get<double>("absolute error tolerance",1e-5);
  rtol_ = mpc_plist_.get<double>("relative error tolerance",1e-5);

  // Instantiate time stepper.
  Teuchos::RCP<Teuchos::ParameterList> bdf1_plist_p =
    Teuchos::rcp(new Teuchos::ParameterList(mpc_plist_.sublist("time integrator")));
  time_stepper_ = Teuchos::rcp(new BDF1TimeIntegrator(this, bdf1_plist_p, solution_));
  time_step_reduction_factor_ = bdf1_plist_p->get<double>("time step reduction factor");

  // Initialize time derivative.
  Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
  solution_dot->PutScalar(0.0);

  // Set initial state for timestepper.
  state_to_solution(S, solution_);
  time_stepper_->set_initial_state(S->time(), solution_, solution_dot);

  // Initialize my timestep size as the min of my sub PK timestep sizes.
  dt_ = MPC::get_dt();
};


// Advance the fully coupled system together.
// -- Advance from state S0 to state S1 at time S0.time + dt.
bool StrongMPC::advance(double dt) {
  state_to_solution(S_next_, solution_);
  n_iter_ = 0;

  // take a bdf timestep
  double h = dt;
  double dt_solver;
  try {
    dt_solver = time_stepper_->time_step(h, solution_);
  } catch (Exceptions::Amanzi_exception &error) {
    if (S_next_->GetMesh()->get_comm()->MyPID() == 0) {
      std::cout << "Timestepper called error: " << error.what() << std::endl;
    }
    if (error.what() == std::string("BDF time step failed") ||
        error.what() == std::string("Cut time step")) {
      // try cutting the timestep
      dt_ = h*time_step_reduction_factor_;
      return true;
    } else {
      throw error;
    }
  }

  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h,S_next_);

  // update the timestep size
  if (dt_solver < dt_ && dt_solver >= h) {
    // We took a smaller step than we recommended, and it worked fine (not
    // suprisingly).  Likely this was due to constraints from other PKs or
    // vis.  Do not reduce our recommendation.
  } else {
    dt_ = dt_solver;
  }

  return false;
};


// computes the non-linear functional g = g(t,u,udot)
void StrongMPC::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // loop over sub-PKs
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    // pull out the old solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_old = u_old->SubVector((*pk)->name());
    if (pk_u_old == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the new solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_new = u_new->SubVector((*pk)->name());
    if (pk_u_new == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the residual sub-vector
    Teuchos::RCP<TreeVector> pk_g = g->SubVector((*pk)->name());
    if (pk_g == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // fill the nonlinear function with each sub-PKs contribution
    (*pk)->fun(t_old, t_new, pk_u_old, pk_u_new, pk_g);
  }

  n_iter_++;
};


// Applies preconditioner to u and returns the result in Pu.
//

// In the Strong MPC case, we build a preconditioner from the preconditioners
// of the sub-PKs.  This PC does not include any coupling terms, so may not be
// any good, but it is automated.
void StrongMPC::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // loop over sub-PKs
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector((*pk)->name());
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the preconditioned u sub-vector
    Teuchos::RCP<TreeVector> pk_Pu = Pu->SubVector((*pk)->name());
    if (pk_Pu == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // Fill the preconditioned u as the block-diagonal product using each sub-PK.
    (*pk)->precon(pk_u, pk_Pu);
  }
};


// computes a norm on u-du and returns the result
//
// For a Strong MPC, the enorm is just the max of the sub PKs enorms.
double StrongMPC::enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du){
  double norm = 0.0;

  // loop over sub-PKs
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector((*pk)->name());
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the du sub-vector
    Teuchos::RCP<const TreeVector> pk_du = du->SubVector((*pk)->name());
    if (pk_du == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // norm is the max of the sub-PK norms
    double tmp_norm = (*pk)->enorm(pk_u, pk_du);
    norm = std::max(norm, tmp_norm);
  }
  return norm;
};


// updates the preconditioner
void StrongMPC::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  PK::solution_to_state(up, S_next_);

  // loop over sub-PKs
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    // pull out the up sub-vector
    Teuchos::RCP<const TreeVector> pk_up = up->SubVector((*pk)->name());
    if (pk_up == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // update precons of each of the sub-PKs
    (*pk)->update_precon(t, pk_up, h);
  };
};

} // namespace
