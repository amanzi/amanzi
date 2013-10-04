/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived StrongMPC class.  Is both a PK and a Model
Evalulator, providing needed methods for BDF time integration of the coupled
system.

Completely automated and generic to any sub PKs, this uses a block diagonal
preconditioner.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#ifndef PKS_MPC_STRONG_MPC_HH_
#define PKS_MPC_STRONG_MPC_HH_

#include <vector>

#include "mpc.hh"
#include "pk_bdf_base.hh"

namespace Amanzi {

// note this looks odd, but StrongMPC is both a MPC within a hierarchy of BDF
// PKs, but it also IS a BDF PK itself, in that it implements the BDF
// interface.
template <class PK_t>
class StrongMPC : public MPC<PK_t>,
                  public PKBDFBase {

public:
  StrongMPC(Teuchos::ParameterList& plist,
            const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPC<PK_t>(plist,soln),
      PKBDFBase(plist,soln) {}

  // Virtual destructor
  virtual ~StrongMPC() {}

  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // StrongMPC is a BDFFnBase
  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk fun().
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // -- enorm for the coupled system
  virtual double enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  // StrongMPC's preconditioner is, by default, just the block-diagonal
  // operator formed by placing the sub PK's preconditioners on the diagonal.
  // -- Apply preconditioner to u and returns the result in Pu.
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Update the preconditioner.
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void changed_solution();

  // -- Admissibility of the solution.
  virtual bool is_admissible(Teuchos::RCP<const TreeVector> u);

  // -- Modify the predictor.
  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> u);

  // -- Modify the correction.
  virtual bool modify_correction(double h, Teuchos::RCP<const TreeVector> res,
          Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

private:
  // factory registration
  static RegisteredPKFactory<StrongMPC> reg_;

};


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
template<class PK_t>
void StrongMPC<PK_t>::setup(const Teuchos::Ptr<State>& S) {
  MPC<PK_t>::setup(S);
  PKBDFBase::setup(S);

  // Set the initial timestep as the min of the sub-pk sizes.
  dt_ = 1.0e99;
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {
    dt_ = std::min<double>(dt_, (*pk)->get_dt());
  }

};


// -----------------------------------------------------------------------------
// Required unique initialize(), just calls both of its base class
// initialize() methods.
// -----------------------------------------------------------------------------
template<class PK_t>
void StrongMPC<PK_t>::initialize(const Teuchos::Ptr<State>& S) {
  // Just calls both subclass's initialize.  NOTE - order is important here --
  // MPC<PK_t> grabs the primary variables from each sub-PK and stuffs
  // them into the solution, which must be done prior to BDFBase initializing
  // the timestepper.

  // Initialize all sub PKs.
  MPC<PK_t>::initialize(S);

  // Initialize my timestepper.
  PKBDFBase::initialize(S);
};


// -----------------------------------------------------------------------------
// Compute the non-linear functional g = g(t,u,udot).
// -----------------------------------------------------------------------------
template<class PK_t>
void StrongMPC<PK_t>::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  solution_to_state(u_new, S_next_);

  // loop over sub-PKs
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {

    // pull out the old solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_old(Teuchos::null);
    if (u_old != Teuchos::null) {
      pk_u_old = u_old->SubVector((*pk)->name());
      if (pk_u_old == Teuchos::null) {
        Errors::Message message("MPC: vector structure does not match PK structure");
        Exceptions::amanzi_throw(message);
      }
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
};


// -----------------------------------------------------------------------------
// Applies preconditioner to u and returns the result in Pu.
// -----------------------------------------------------------------------------
template<class PK_t>
void StrongMPC<PK_t>::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // loop over sub-PKs
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {

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

//   std::cout<<*(((Pu->SubVector("flow"))->Data())->ViewComponent("cell", false));
//   cout<<"Exit from StrongMPC precon\n";
//   exit(0);
};


// -----------------------------------------------------------------------------
// Compute a norm on u-du and returns the result.
// For a Strong MPC, the enorm is just the max of the sub PKs enorms.
// -----------------------------------------------------------------------------
template<class PK_t>
double StrongMPC<PK_t>::enorm(Teuchos::RCP<const TreeVector> u,
                        Teuchos::RCP<const TreeVector> du){
  double norm = 0.0;

  // loop over sub-PKs
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {

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


// -----------------------------------------------------------------------------
// Update the preconditioner.
// -----------------------------------------------------------------------------
template<class PK_t>
void StrongMPC<PK_t>::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  PKDefaultBase::solution_to_state(up, S_next_);

  // loop over sub-PKs
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {

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


// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time integration
// scheme is changing the value of the solution in state.
// -----------------------------------------------------------------------------
template<class PK_t>
void StrongMPC<PK_t>::changed_solution() {
  // loop over sub-PKs
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {
    (*pk)->changed_solution();
  }
};

// -----------------------------------------------------------------------------
// Check admissibility of each sub-pk
// -----------------------------------------------------------------------------
template<class PK_t>
bool StrongMPC<PK_t>::is_admissible(Teuchos::RCP<const TreeVector> u) {
  // First ensure each PK thinks we are admissible -- this will ensure
  // the residual can at least be evaluated.
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {

    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector((*pk)->name());
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    if (!(*pk)->is_admissible(pk_u)) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "PK " << (*pk)->name() << " is not admissible." << std::endl;

      return false;
    }
  }

  // If that worked, check backtracking admissility.
  return PKBDFBase::is_admissible(u);
};


// -----------------------------------------------------------------------------
// Modify predictor from each sub pk.
// -----------------------------------------------------------------------------
template<class PK_t>
bool StrongMPC<PK_t>::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  // loop over sub-PKs
  bool modified = false;
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {

    // pull out the u sub-vector
    Teuchos::RCP<TreeVector> pk_u = u->SubVector((*pk)->name());
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    modified |= (*pk)->modify_predictor(h,pk_u);
  }
  return modified;
};


// -----------------------------------------------------------------------------
// Modify correction from each sub pk.
// -----------------------------------------------------------------------------
template<class PK_t>
bool StrongMPC<PK_t>::modify_correction(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  // loop over sub-PKs
  bool modified = false;
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {

    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector((*pk)->name());
    Teuchos::RCP<const TreeVector> pk_res = res->SubVector((*pk)->name());
    Teuchos::RCP<TreeVector> pk_du = du->SubVector((*pk)->name());

    if (pk_u == Teuchos::null || pk_du == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    modified |= (*pk)->modify_correction(h, pk_res, pk_u, pk_du);
  }
  return modified;
};


} // close namespace Amanzi

#endif
