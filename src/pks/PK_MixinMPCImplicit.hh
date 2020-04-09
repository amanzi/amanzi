/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A mixin for an implicitly coupling MPC

/*!

This almost certainly needs to be combined with PK_MixinImplicit and
PK_MixinMPC.

E.g, the following should be instantiable as a globally implicit MPC:

PK_Implicit_Adaptor<PK_MixinImplicit<PK_MixinMPCImplicit<PK_Default,
PK_Implicit<TreeVector> > > >

*/

#ifndef AMANZI_PK_MIXIN_MPC_IMPLICIT_HH_
#define AMANZI_PK_MIXIN_MPC_IMPLICIT_HH_

#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_MixinMPC.hh"

namespace Amanzi {

template <class Base_t, class PK_Contained_t>
class PK_MixinMPCImplicit : public PK_MixinMPC<Base_t, PK_Contained_t> {
 public:
  using PK_MixinMPC<Base_t, PK_Contained_t>::template PK_MixinMPC;

  // IMPLEMENT ME! Is this different from MPC?
  void ConstructChildren() {}

  // PK methods
  // the BDFfnBase interface
  // computes the non-linear functional f = f(t,u,udot)
  void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                     Teuchos::RCP<TreeVector> u_new,
                     Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                          Teuchos::RCP<TreeVector> Pu);

  // computes a norm on u-du and returns the result
  double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<const TreeVector> du);

  // updates the preconditioner
  void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // check the admissibility of a solution
  // override final with the actual admissibility check
  bool IsAdmissible(Teuchos::RCP<const TreeVector> up);

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                       Teuchos::RCP<TreeVector> u);

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du);

  // update the continuation parameter
  void UpdateContinuationParameter(double lambda);

  // calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  void ChangedSolution();

 protected:
  using PK_MixinMPC<Base_t, PK_Contained_t>::sub_pks_;
};

// PK methods
// the BDFfnBase interface
// computes the non-linear functional f = f(t,u,udot)
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPCImplicit<Base_t, PK_Contained_t>::FunctionalResidual(
  double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
  Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f)
{
  int i = 0;
  for (auto& pk : sub_pks_) {
    pk->FunctionalResidual(
      t_old, t_new, u_old->SubVector(i), u_new->SubVector(i), f->SubVector(i));
    ++i;
  }
}

// applies preconditioner to u and returns the result in Pu
template <class Base_t, class PK_Contained_t>
int
PK_MixinMPCImplicit<Base_t, PK_Contained_t>::ApplyPreconditioner(
  Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  int ierr;
  int i = 0;
  for (auto& pk : sub_pks_) {
    ierr += pk->ApplyPreconditioner(u->SubVector(i), Pu->SubVector(i));
    ++i;
  }
  return ierr;
}

// computes a norm on u-du and returns the result
template <class Base_t, class PK_Contained_t>
double
PK_MixinMPCImplicit<Base_t, PK_Contained_t>::ErrorNorm(
  Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  double error;
  int i = 0;
  for (auto& pk : sub_pks_) {
    double err = pk->ErrorNorm(u->SubVector(i), du->SubVector(i));
    error = std::max(err, error);
    ++i;
  }
  return error;
}

// updates the preconditioner
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPCImplicit<Base_t, PK_Contained_t>::UpdatePreconditioner(
  double t, Teuchos::RCP<const TreeVector> up, double h)
{
  int i = 0;
  for (auto& pk : sub_pks_) {
    pk->UpdatePreconditioner(t, up->SubVector(i), h);
    ++i;
  }
}

// check the admissibility of a solution
// override final with the actual admissibility check
template <class Base_t, class PK_Contained_t>
bool
PK_MixinMPCImplicit<Base_t, PK_Contained_t>::IsAdmissible(
  Teuchos::RCP<const TreeVector> up)
{
  int i = 0;
  for (auto& pk : sub_pks_) {
    bool admissible = pk->IsAdmissible(up->SubVector(i));
    if (!admissible) return admissible;
    ++i;
  }
  return true;
}

// possibly modifies the predictor that is going to be used as a
// starting value for the nonlinear solve in the time integrator,
// the time integrator will pass the predictor that is computed
// using extrapolation and the time step that is used to compute
// this predictor this function returns true if the predictor was
// modified, false if not
template <class Base_t, class PK_Contained_t>
bool
PK_MixinMPCImplicit<Base_t, PK_Contained_t>::ModifyPredictor(
  double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u)
{
  bool modified = false;
  int i = 0;
  for (auto& pk : sub_pks_) {
    modified |= pk->ModifyPredictor(h, u0->SubVector(i), u->SubVector(i));
    ++i;
  }
  return modified;
}

// possibly modifies the correction, after the nonlinear solver (NKA)
// has computed it, will return true if it did change the correction,
// so that the nonlinear iteration can store the modified correction
// and pass it to NKA so that the NKA space can be updated
template <class Base_t, class PK_Contained_t>
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
PK_MixinMPCImplicit<Base_t, PK_Contained_t>::ModifyCorrection(
  double h, Teuchos::RCP<const TreeVector> res,
  Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du)
{
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult modified =
    AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;

  int i = 0;
  for (auto& pk : sub_pks_) {
    auto local_modified = pk->ModifyCorrection(
      h, res->SubVector(i), u->SubVector(i), du->SubVector(i));
    ++i;
    modified = std::max(modified, local_modified);
  }
  return modified;
}

// update the continuation parameter
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPCImplicit<Base_t, PK_Contained_t>::UpdateContinuationParameter(
  double lambda)
{
  for (auto& pk : sub_pks_) pk->UpdateContinuationParameter(lambda);
}

// calling this indicates that the time
// integration scheme is changing the value of the solution in
// state.
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPCImplicit<Base_t, PK_Contained_t>::ChangedSolution()
{
  for (auto& pk : sub_pks_) pk->ChangedSolution();
}

} // namespace Amanzi

#endif
