/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Multi process coupler for globally implicit (strong) coupling.
/*!

Globally implicit coupling solves all sub-PKs as a single system of equations.  This can be completely automated when all PKs are also `PK: BDF`_ PKs, using a block-diagonal preconditioner where each diagonal block is provided by its own sub-PK.

.. _strong-mpc-spec:
.. admonition:: strong-mpc-spec

    INCLUDES:

    - ``[mpc-spec]`` *Is a* MPC_.
    - ``[pk-bdf-default-spec]`` *Is a* `PK: BDF`_.

*/

#pragma once

#include <vector>

#include "MPC.hh"
#include "PK_BDF_Default.hh"

namespace Amanzi {

// note this looks odd, but MPCStrong is both a MPC within a hierarchy of BDF
// PKs, but it also IS a BDF PK itself, in that it implements the BDF
// interface and can be implicitly integrated.
//
// Does this need to be templated?  Can't SubPK_type = PK_BDF_Default? Or do we
// need it to sometimes be PK_PhysicalBDF_Default? --ETC
template <class SubPK_type>
class MPCStrong : public MPC<PK_BDF_Default, SubPK_type> {
 public:
  MPCStrong(const Comm_ptr_type& comm,
            Teuchos::ParameterList& pk_tree,
            const Teuchos::RCP<Teuchos::ParameterList>& global_list,
            const Teuchos::RCP<State>& S);

  virtual void modifyParameterList() override;
  virtual const std::string& getType() const override { return pk_type_; }

  virtual void markChangedSolutionPK(const Tag& tag) override;

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk FunctionalResidual().
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g) override;

  // -- enorm for the coupled system
  virtual double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;

  // MPCStrong's preconditioner is, by default, just the block-diagonal
  // operator formed by placing the sub PK's preconditioners on the diagonal.
  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // -- Admissibility of the solution.
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> u) override;

  // -- Modify the predictor.
  virtual bool
  ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override;

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h,
                   Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override;

 protected:
  using MPC<PK_BDF_Default, SubPK_type>::sub_pks_;
  using MPC<PK_BDF_Default, SubPK_type>::global_list_;
  using MPC<PK_BDF_Default, SubPK_type>::pk_tree_;
  using MPC<PK_BDF_Default, SubPK_type>::pks_list_;
  using MPC<PK_BDF_Default, SubPK_type>::plist_;
  using MPC<PK_BDF_Default, SubPK_type>::tag_next_;
  using MPC<PK_BDF_Default, SubPK_type>::S_;
  using MPC<PK_BDF_Default, SubPK_type>::vo_;

  static const std::string pk_type_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCStrong> reg_;
};

//
// Class implementation
//

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template <class SubPK_type>
MPCStrong<SubPK_type>::MPCStrong(const Comm_ptr_type& comm,
                           Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S)
  : MPC<PK_BDF_Default, SubPK_type>(comm, pk_tree, global_list, S)
{
  MPC<PK_BDF_Default, SubPK_type>::createSubPKs_();
}


template<class SubPK_type>
void
MPCStrong<SubPK_type>::modifyParameterList()
{
  // push on a parameter to indicate that sub-pks need not assemble their
  // operators, as we will do that here (or above here)
  auto pk_order = plist_->template get<Teuchos::Array<std::string>>("PKs order");
  for (const auto& pk_name : pk_order) {
    pks_list_->sublist(pk_name).set("strongly coupled PK", true);
  }

  MPC<PK_BDF_Default, SubPK_type>::modifyParameterList();
}


// -----------------------------------------------------------------------------
// Compute the non-linear functional g = g(t,u,udot).
// -----------------------------------------------------------------------------
template <class SubPK_type>
void
MPCStrong<SubPK_type>::FunctionalResidual(double t_old,
                                    double t_new,
                                    Teuchos::RCP<TreeVector> u_old,
                                    Teuchos::RCP<TreeVector> u_new,
                                    Teuchos::RCP<TreeVector> g)
{
  this->moveSolutionToState(*u_new, tag_next_);

  // loop over sub-PKs
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the old solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_old(Teuchos::null);
    if (u_old != Teuchos::null) {
      pk_u_old = u_old->getSubVector(i);
      if (pk_u_old == Teuchos::null) {
        Errors::Message message("MPC: vector structure does not match PK structure");
        Exceptions::amanzi_throw(message);
      }
    }

    // pull out the new solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_new = u_new->getSubVector(i);
    if (pk_u_new == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the residual sub-vector
    Teuchos::RCP<TreeVector> pk_g = g->getSubVector(i);
    if (pk_g == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // fill the nonlinear function with each sub-PKs contribution
    sub_pks_[i]->FunctionalResidual(t_old, t_new, pk_u_old, pk_u_new, pk_g);
  }
};


// -----------------------------------------------------------------------------
// Applies preconditioner to u and returns the result in Pu.
// -----------------------------------------------------------------------------
template <class SubPK_type>
int
MPCStrong<SubPK_type>::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  // loop over sub-PKs
  int ierr = 0;
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->getSubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the preconditioned u sub-vector
    Teuchos::RCP<TreeVector> pk_Pu = Pu->getSubVector(i);
    if (pk_Pu == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // Fill the preconditioned u as the block-diagonal product using each sub-PK.
    int icur_err = sub_pks_[i]->ApplyPreconditioner(pk_u, pk_Pu);
    ierr += icur_err;
  }
  return ierr;
};


// -----------------------------------------------------------------------------
// Compute a norm on u-du and returns the result.
// For a Strong MPC, the enorm is just the max of the sub PKs enorms.
// -----------------------------------------------------------------------------
template <class SubPK_type>
double
MPCStrong<SubPK_type>::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  double norm = 0.0;

  // loop over sub-PKs
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->getSubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the du sub-vector
    Teuchos::RCP<const TreeVector> pk_du = du->getSubVector(i);
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
template <class SubPK_type>
void
MPCStrong<SubPK_type>::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  this->moveSolutionToState(*up, tag_next_);

  // loop over sub-PKs
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the up sub-vector
    Teuchos::RCP<const TreeVector> pk_up = up->getSubVector(i);
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
template <class SubPK_type>
void
MPCStrong<SubPK_type>::markChangedSolutionPK(const Tag& tag)
{
  // loop over sub-PKs
  for (auto& pk : sub_pks_) pk->markChangedSolutionPK(tag);
};


// -----------------------------------------------------------------------------
// Check admissibility of each sub-pk
// -----------------------------------------------------------------------------
template <class SubPK_type>
bool
MPCStrong<SubPK_type>::IsAdmissible(Teuchos::RCP<const TreeVector> u)
{
  // First ensure each PK thinks we are admissible -- this will ensure
  // the residual can at least be evaluated.
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->getSubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    if (!sub_pks_[i]->IsAdmissible(pk_u)) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "PK " << sub_pks_[i]->getName() << " is not admissible." << std::endl;

      return false;
    }
  }

  // If that worked, check backtracking admissility.
  //  return PKBDFBase::IsAdmissible(u);
  return true;
};


// -----------------------------------------------------------------------------
// Modify predictor from each sub pk.
// -----------------------------------------------------------------------------
template <class SubPK_type>
bool
MPCStrong<SubPK_type>::ModifyPredictor(double h,
                                 Teuchos::RCP<const TreeVector> u0,
                                 Teuchos::RCP<TreeVector> u)
{
  // loop over sub-PKs
  bool modified = false;
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u0 = u0->getSubVector(i);
    Teuchos::RCP<TreeVector> pk_u = u->getSubVector(i);
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
template <class SubPK_type>
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MPCStrong<SubPK_type>::ModifyCorrection(double h,
                                  Teuchos::RCP<const TreeVector> res,
                                  Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> du)
{
  // loop over sub-PKs
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult modified =
    AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->getSubVector(i);
    Teuchos::RCP<const TreeVector> pk_res = res->getSubVector(i);
    Teuchos::RCP<TreeVector> pk_du = du->getSubVector(i);

    if (pk_u == Teuchos::null || pk_du == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    modified = std::max(modified, sub_pks_[i]->ModifyCorrection(h, pk_res, pk_u, pk_du));
  }
  return modified;
};


} // namespace Amanzi
