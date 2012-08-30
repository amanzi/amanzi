/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "boost/math/special_functions/fpclassify.hpp"

#include "field_evaluator.hh"
#include "two_phase.hh"

namespace Amanzi {
namespace Energy {

#define DEBUG_FLAG 0

// TwoPhase is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void TwoPhase::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  Teuchos::RCP<CompositeVector> u = u_new->data();
#if DEBUG_FLAG
  std::cout << "Two-Phase Residual calculation:" << std::endl;
  std::cout << "  T0: " << (*u)("cell",0) << " " << (*u)("face",3) << std::endl;
  std::cout << "  T1: " << (*u)("cell",99) << " " << (*u)("face",497) << std::endl;
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
  ApplyDiffusion_(S_next_, res);
#if DEBUG_FLAG
  std::cout << "  res0 (after diffusion): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
  std::cout << "  res1 (after diffusion): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
#endif

  // accumulation term
  AddAccumulation_(res);
#if DEBUG_FLAG
  std::cout << "  res0 (after accumulation): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
  std::cout << "  res1 (after accumulation): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
#endif

  // advection term, explicit
  AddAdvection_(S_inter_, res, false);
#if DEBUG_FLAG
  std::cout << "  res0 (after advection): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
  std::cout << "  res1 (after advection): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
#endif
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void TwoPhase::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
#if DEBUG_FLAG
  std::cout << "Precon application:" << std::endl;
  std::cout << "  T0: " << (*u->data())("cell",0) << " " << (*u->data())("face",3) << std::endl;
  std::cout << "  T1: " << (*u->data())("cell",99) << " " << (*u->data())("face",497) << std::endl;
#endif

  preconditioner_->ApplyInverse(*u->data(), Pu->data());

#if DEBUG_FLAG
  std::cout << "  PC*T0: " << (*Pu->data())("cell",0) << " " << (*Pu->data())("face",3) << std::endl;
  std::cout << "  PC*T1: " << (*Pu->data())("cell",99) << " " << (*Pu->data())("face",497) << std::endl;
#endif
};


// -----------------------------------------------------------------------------
// Compute a norm on (u,du)
// -----------------------------------------------------------------------------
double TwoPhase::enorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) {
  Teuchos::RCP<const CompositeVector> temp = u->data();
  Teuchos::RCP<const CompositeVector> dtemp = du->data();


  double enorm_val_cell = 0.0;
  for (int c=0; c!=temp->size("cell",false); ++c) {
    if (boost::math::isnan<double>((*dtemp)("cell",c))) {
      std::cout << "Cutting time step due to NaN in correction." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }

    double tmp = abs((*dtemp)("cell",c)) / (atol_ + rtol_ * abs((*temp)("cell",c)));
    enorm_val_cell = std::max<double>(enorm_val_cell, tmp);
    //    printf("cell: %5i %14.7e %14.7e\n",lcv,(*(*dtemp_vec)(0))[lcv],tmp);
  }

  double enorm_val_face = 0.0;
  for (int f=0; f!=temp->size("face",false); ++f) {
    if (boost::math::isnan<double>((*dtemp)("face",f))) {
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }

    double tmp = abs((*dtemp)("face",f)) / (atol_ + rtol_ * abs((*temp)("face",f)));
    enorm_val_face = std::max<double>(enorm_val_face, tmp);
    //    printf("face: %5i %14.7e %14.7e\n",lcv,(*(*fdtemp_vec)(0))[lcv],tmp);
  }


  //  std::cout.precision(15);
  //  std::cout << "Temperature enorm (cell, face): " << std::scientific << enorm_val_cell
  //            << " / " << std::scientific << enorm_val_face << std::endl;

  double enorm_val = std::max<double>(enorm_val_cell, enorm_val_face);
#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return enorm_val;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void TwoPhase::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
#if DEBUG_FLAG
  std::cout << "Precon update at t = " << t << std::endl;
#endif

  // update state with the solution up.
  S_next_->set_time(t);
  PK::solution_to_state(up, S_next_);

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // div K_e grad u
  S_next_->GetFieldEvaluator("thermal_conductivity")
    ->HasFieldChanged(S_next_.ptr(), "energy_pk");
  Teuchos::RCP<const CompositeVector> thermal_conductivity =
    S_next_->GetFieldData("thermal_conductivity");

  preconditioner_->CreateMFDstiffnessMatrices(*thermal_conductivity);
  preconditioner_->CreateMFDrhsVectors();

  // update with accumulation terms
  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator("energy")
      ->HasFieldDerivativeChanged(S_next_.ptr(), "energy_pk", "temperature");

  // -- get the accumulation deriv
  Teuchos::RCP<const CompositeVector> de_dT =
      S_next_->GetFieldData("denergy_dtemperature");
  Teuchos::RCP<const CompositeVector> temp =
      S_next_->GetFieldData("temperature");

  // -- get the matrices/rhs that need updating
  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();

  // -- update the diagonal
  int ncells = temp->size("cell");
  for (int c=0; c!=ncells; ++c) {
    Acc_cells[c] += (*de_dT)("cell",c) / h;
    //    Fc_cells[c] += (*de_dT)("cell",c) / h * (*temp)("cell",c);
  }

  // -- assemble
  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  preconditioner_->AssembleGlobalMatrices();

  // -- form and prep the Schur complement for inversion
  preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
  preconditioner_->UpdatePreconditioner();

};

} // namespace Energy
} // namespace Amanzi
