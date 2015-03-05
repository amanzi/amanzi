/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "Epetra_Vector.h"
#include "advection_diffusion.hh"
#include "Op.hh"

namespace Amanzi {
namespace Energy {

// AdvectionDiffusion is a BDFFnBase
// computes the non-linear functional g = g(t,u,udot)
void AdvectionDiffusion::Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                 Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(*u_new, S_next_);

  bc_temperature_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  Teuchos::RCP<CompositeVector> u = u_new->Data();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Residual calculation:" << std::endl;
    *vo_->os() << "  u: " << (*u)("cell",0) << std::endl;
  }

  // get access to the solution
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  ApplyDiffusion_(S_next_, res);
  if (vo_->os_OK(Teuchos::VERB_HIGH)) 
    *vo_->os() << "  res (after diffusion): " << (*res)("cell",0) << std::endl;

  // accumulation term
  AddAccumulation_(res);
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  res (after accumulation): " << (*res)("cell",0) << std::endl;

  // advection term, explicit
  AddAdvection_(S_inter_, res, true);
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  res (after advection): " << (*res)("cell",0) << std::endl;
};

// applies preconditioner to u and returns the result in Pu
void AdvectionDiffusion::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Precon application:" << std::endl;
    *vo_->os() << "  u: " << (*u->Data())("cell",0) << std::endl;
  }
  // preconditioner for accumulation only:
  //  *Pu = *u;

  // MFD ML preconditioner
  preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Pu: " << (*Pu->Data())("cell",0)  << std::endl;
};


// updates the preconditioner
void AdvectionDiffusion::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(*up, S_next_);

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // div K_e grad u
  Teuchos::RCP<const CompositeVector> thermal_conductivity =
      S_next_->GetFieldData("thermal_conductivity");
  preconditioner_->Init();
  preconditioner_diff_->Setup(thermal_conductivity, Teuchos::null);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update with accumulation terms
  Teuchos::RCP<const CompositeVector> cell_volume =
    S_next_->GetFieldData("cell_volume");

  std::vector<double>& Acc_cells = preconditioner_acc_->local_matrices()->vals;
  int ncells = cell_volume->size("cell",false);
  for (int c=0; c!=ncells; ++c) {
    double factor = (*cell_volume)("cell",c);
    Acc_cells[c] += factor/h;
  }

  preconditioner_diff_->ApplyBCs(true);
  preconditioner_->AssembleMatrix();
};


} // namespace Energy
} // namespace Amanzi
