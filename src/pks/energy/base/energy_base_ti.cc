/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "boundary_function.hh"
#include "field_evaluator.hh"
#include "energy_base.hh"

namespace Amanzi {
namespace Energy {

#define DEBUG_FLAG 1
#define MORE_DEBUG_FLAG 0

// EnergyBase is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void EnergyBase::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  niter_++;

  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();
  double h = t_new - t_old;
  ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  Teuchos::RCP<CompositeVector> u = u_new->data();

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Residual calculation:" << std::endl;
    *out_ << "  T0: " << (*u)("cell",0) << " " << (*u)("face",3) << std::endl;
    *out_ << "  T1: " << (*u)("cell",99) << " " << (*u)("face",497) << std::endl;
  }
#endif

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_temperature_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  ApplyDiffusion_(S_next_.ptr(), res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after diffusion): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
    *out_ << "  res1 (after diffusion): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
  }
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after accumulation): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
    *out_ << "  res1 (after accumulation): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
  }
#endif

  // advection term, implicit
  AddAdvection_(S_next_.ptr(), res.ptr(), false);
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after advection): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
    *out_ << "  res1 (after advection): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
  }
#endif

  // source terms
  AddSources_(S_next_.ptr(), res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after sources): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
    *out_ << "  res1 (after sources): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
  }
#endif

  // Dump residual to state for visual debugging.
#if MORE_DEBUG_FLAG
  if (niter_ < 23) {
    std::stringstream namestream;
    namestream << domain_prefix_ << "energy_residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::stringstream solnstream;
    solnstream << domain_prefix_ << "energy_solution_" << niter_;
    *S_next_->GetFieldData(solnstream.str(),name_) = *u;
  }
#endif

};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void EnergyBase::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
#if DEBUG_FLAG
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  T0: " << (*u->data())("cell",0) << " " << (*u->data())("face",3) << std::endl;
    *out_ << "  T1: " << (*u->data())("cell",99) << " " << (*u->data())("face",497) << std::endl;
  }
#endif

  // apply the preconditioner
  preconditioner_->ApplyInverse(*u->data(), Pu->data().ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  PC*T0: " << (*Pu->data())("cell",0) << " " << (*Pu->data())("face",3) << std::endl;
    *out_ << "  PC*T1: " << (*Pu->data())("cell",99) << " " << (*Pu->data())("face",497) << std::endl;
  }
#endif
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void EnergyBase::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon update at t = " << t << std::endl;
  }

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);

  PKDefaultBase::solution_to_state(up, S_next_);

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // div K_e grad u
  S_next_->GetFieldEvaluator(conductivity_key_)
    ->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> conductivity =
    S_next_->GetFieldData(conductivity_key_);

  preconditioner_->CreateMFDstiffnessMatrices(conductivity.ptr());
  preconditioner_->CreateMFDrhsVectors();

  // update with accumulation terms
  // -- update the accumulation derivatives, de/dT
  S_next_->GetFieldEvaluator(energy_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);

  // -- get the accumulation deriv
  const Epetra_MultiVector& de_dT = *S_next_->GetFieldData(de_dT_key_)
      ->ViewComponent("cell",false);

  // -- get the matrices/rhs that need updating
  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();

  // -- update the diagonal
  int ncells = de_dT.MyLength();
  for (int c=0; c!=ncells; ++c) {
    Acc_cells[c] += de_dT[0][c] / h;
  }

  // Apply boundary conditions.
  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  // Assemble
  if (assemble_preconditioner_) {
    // -- assemble
    preconditioner_->AssembleGlobalMatrices();
    // -- form and prep the Schur complement for inversion
    preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
    preconditioner_->UpdatePreconditioner();
  }
};


double EnergyBase::enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) {

  Teuchos::RCP<const CompositeVector> res = du->data();
  const Epetra_MultiVector& res_c = *res->ViewComponent("cell",false);
  const Epetra_MultiVector& res_f = *res->ViewComponent("face",false);

  const Epetra_MultiVector& cv = *S_next_->GetFieldData(cell_vol_key_)
      ->ViewComponent("cell",false);
  const CompositeVector& temp = *u->data();
  double h = S_next_->time() - S_inter_->time();

  // Cell error is based upon error in energy conservation relative to
  // a characteristic energy
  double enorm_cell(0.);
  int ncells = res_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    double tmp = std::abs(h*res_c[0][c]) / (atol_+rtol_* (cv[0][c]*2.e6));
    enorm_cell = std::max<double>(enorm_cell, tmp);
  }

  // Face error is based upon temperature?  This is unclear!
  double enorm_face(0.);
  int nfaces = res_f.MyLength();
  for (int f=0; f!=nfaces; ++f) {
    double tmp = std::abs(res_f[0][f]) / (atol_+rtol_*273.15);
    enorm_face = std::max<double>(enorm_face, tmp);
  }

  // Write out Inf norms too.
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    double infnorm_c(0.), infnorm_f(0.);
    res_c.NormInf(&infnorm_c);
    res_f.NormInf(&infnorm_f);

#ifdef HAVE_MPI
    double buf_c(enorm_cell), buf_f(enorm_face);
    MPI_Allreduce(&buf_c, &enorm_cell, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&buf_f, &enorm_face, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    *out_ << "ENorm (cells) = " << enorm_cell << " (" << infnorm_c << ")  " << std::endl;
    *out_ << "ENorm (faces) = " << enorm_face << " (" << infnorm_f << ")  " << std::endl;
  }

  // Communicate and take the max.
  double enorm_val(std::max<double>(enorm_face, enorm_cell));
#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return enorm_val;
};


} // namespace Energy
} // namespace Amanzi
