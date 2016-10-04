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
<<<<<<< HEAD
  StrongMPC(Teuchos::Ptr<State> S,const Teuchos::RCP<Teuchos::ParameterList>& plist,
=======
  StrongMPC(const Teuchos::RCP<Teuchos::ParameterList>& plist,
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
            Teuchos::ParameterList& FElist,
            const Teuchos::RCP<TreeVector>& soln);

  // Virtual destructor
  virtual ~StrongMPC() {}

  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {
    PKBDFBase::commit_state(dt,S);
    MPC<PK_t>::commit_state(dt, S);
  }
  
  // StrongMPC is a BDFFnBase
  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk Functional().
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // -- enorm for the coupled system
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  // StrongMPC's preconditioner is, by default, just the block-diagonal
  // operator formed by placing the sub PK's preconditioners on the diagonal.
  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
<<<<<<< HEAD
  virtual void ChangedSolution(const Teuchos::Ptr<State>& S);
=======
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
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

protected:
  using MPC<PK_t>::sub_pks_;

private:
  // factory registration
  static RegisteredPKFactory<StrongMPC> reg_;

};

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template<class PK_t>
<<<<<<< HEAD
StrongMPC<PK_t>::StrongMPC(Teuchos::Ptr<State> S, const Teuchos::RCP<Teuchos::ParameterList>& plist,
                           Teuchos::ParameterList& FElist,
                           const Teuchos::RCP<TreeVector>& soln) :
  PKDefaultBase(S, plist, FElist, soln),
  MPC<PK_t>(S, plist, FElist, soln),
  PKBDFBase(S, plist, FElist, soln) {}
=======
StrongMPC<PK_t>::StrongMPC(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                           Teuchos::ParameterList& FElist,
                           const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln),
    MPC<PK_t>(plist, FElist, soln),
    PKBDFBase(plist, FElist, soln) {}
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
template<class PK_t>
void StrongMPC<PK_t>::setup(const Teuchos::Ptr<State>& S) {
  // tweak the sub-PK parameter lists
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(plist_, "PKs");
  for (Teuchos::ParameterList::ConstIterator param=pks_list->begin();
       param!=pks_list->end(); ++param) {
    std::string pname = param->first;
    if (pks_list->isSublist(pname)) {
      pks_list->sublist(pname).set("strongly coupled PK", true);
    } else {
      ASSERT(0);
    }
  }

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
void StrongMPC<PK_t>::Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  solution_to_state(*u_new, S_next_);

  // loop over sub-PKs
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
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
    sub_pks_[i]->Functional(t_old, t_new, pk_u_old, pk_u_new, pk_g);
  }
};


// -----------------------------------------------------------------------------
// Applies preconditioner to u and returns the result in Pu.
// -----------------------------------------------------------------------------
template<class PK_t>
int StrongMPC<PK_t>::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // loop over sub-PKs
  int ierr = 0;
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
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
    int icur_err = sub_pks_[i]->ApplyPreconditioner(pk_u, pk_Pu);
    ierr += icur_err;
  }

//   std::cout<<*(((Pu->SubVector("flow"))->Data())->ViewComponent("cell", false));
//   cout<<"Exit from StrongMPC precon\n";
//   exit(0);
  return ierr;
};


// -----------------------------------------------------------------------------
// Compute a norm on u-du and returns the result.
// For a Strong MPC, the enorm is just the max of the sub PKs enorms.
// -----------------------------------------------------------------------------
template<class PK_t>
double StrongMPC<PK_t>::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                        Teuchos::RCP<const TreeVector> du){
  double norm = 0.0;

  // loop over sub-PKs
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the du sub-vector
    Teuchos::RCP<const TreeVector> pk_du = du->SubVector(i);
    if (pk_du == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // norm is the max of the sub-PK norms
    double tmp_norm = sub_pks_[i]->ErrorNorm(pk_u, pk_du);
    norm = std::max(norm, tmp_norm);
  }
  return norm;
};


// -----------------------------------------------------------------------------
// Update the preconditioner.
// -----------------------------------------------------------------------------
template<class PK_t>
void StrongMPC<PK_t>::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  PKDefaultBase::solution_to_state(*up, S_next_);

  // loop over sub-PKs
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
    // pull out the up sub-vector
    Teuchos::RCP<const TreeVector> pk_up = up->SubVector(i);
    if (pk_up == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // update precons of each of the sub-PKs
    sub_pks_[i]->UpdatePreconditioner(t, pk_up, h);
  };
};


// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time integration
// scheme is changing the value of the solution in state.
// -----------------------------------------------------------------------------
template<class PK_t>
<<<<<<< HEAD
void StrongMPC<PK_t>::ChangedSolution(const Teuchos::Ptr<State>& S) {
  // loop over sub-PKs
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {
    (*pk)->ChangedSolution(S);
  }
};

  
// -----------------------------------------------------------------------------
// Calling this indicates that the time integration scheme is changing
// the value of the solution in state.
// -----------------------------------------------------------------------------
template<class PK_t>
=======
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
void StrongMPC<PK_t>::ChangedSolution() {
  // loop over sub-PKs
  for (typename MPC<PK_t>::SubPKList::iterator pk = MPC<PK_t>::sub_pks_.begin();
       pk != MPC<PK_t>::sub_pks_.end(); ++pk) {
    (*pk)->ChangedSolution();
  }
};

<<<<<<< HEAD

=======
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
// -----------------------------------------------------------------------------
// Check admissibility of each sub-pk
// -----------------------------------------------------------------------------
template<class PK_t>
bool StrongMPC<PK_t>::IsAdmissible(Teuchos::RCP<const TreeVector> u) {
  // First ensure each PK thinks we are admissible -- this will ensure
  // the residual can at least be evaluated.
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    if (!sub_pks_[i]->IsAdmissible(pk_u)) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "PK " << sub_pks_[i]->name() << " is not admissible." << std::endl;

      return false;
    }
  }

  // If that worked, check backtracking admissility.
  return PKBDFBase::IsAdmissible(u);
};


// -----------------------------------------------------------------------------
// Modify predictor from each sub pk.
// -----------------------------------------------------------------------------
template<class PK_t>
bool StrongMPC<PK_t>::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  // loop over sub-PKs
  bool modified = false;
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
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
};


// -----------------------------------------------------------------------------
// Modify correction from each sub pk.
// -----------------------------------------------------------------------------
template<class PK_t>
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
    StrongMPC<PK_t>::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                                      Teuchos::RCP<const TreeVector> u,
                                      Teuchos::RCP<TreeVector> du) {
  // loop over sub-PKs
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      modified = AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
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
};


} // close namespace Amanzi

#endif
