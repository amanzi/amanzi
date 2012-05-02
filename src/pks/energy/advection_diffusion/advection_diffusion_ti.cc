/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "Epetra_Vector.h"
#include "advection_diffusion.hh"

namespace Amanzi {
namespace Energy {

// AdvectionDiffusion is a BDFFnBase
// computes the non-linear functional g = g(t,u,udot)
void AdvectionDiffusion::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                 Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(u_old, S_inter_);
  solution_to_state(u_new, S_next_);

  bc_temperature_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  Teuchos::RCP<CompositeVector> u = u_new->data();
  std::cout << "Residual calculation:" << std::endl;
  std::cout << "  u: " << (*u)("cell",0,0) << " " << (*u)("face",0,0) << std::endl;

  // get access to the solution
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  ApplyDiffusion_(S_next_, res);
  std::cout << "  res (after diffusion): " << (*res)("cell",0,0) << " " << (*res)("face",0,0) << std::endl;

  // accumulation term
  AddAccumulation_(res);
  std::cout << "  res (after accumulation): " << (*res)("cell",0,0) << " " << (*res)("face",0,0) << std::endl;

  // advection term, explicit
  AddAdvection_(S_inter_, res, true);
  std::cout << "  res (after advection): " << (*res)("cell",0,0) << " " << (*res)("face",0,0) << std::endl;
};

// applies preconditioner to u and returns the result in Pu
void AdvectionDiffusion::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  std::cout << "Precon application:" << std::endl;
  std::cout << "  u: " << (*u->data())("cell",0,0) << " " << (*u->data())("face",0,0) << std::endl;
  // preconditioner for accumulation only:
  //  *Pu = *u;

  // MFD ML preconditioner
  preconditioner_->ApplyInverse(*u->data(), Pu->data());

  std::cout << "  Pu: " << (*Pu->data())("cell",0,0) << " " << (*Pu->data())("face",0,0) << std::endl;
};

// computes a norm on u-du and returns the result
double AdvectionDiffusion::enorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) {
  double enorm_val = 0.0;
  Teuchos::RCP<const Epetra_MultiVector> temp_vec = u->data()->ViewComponent("cell", false);
  Teuchos::RCP<const Epetra_MultiVector> dtemp_vec = du->data()->ViewComponent("cell", false);

  for (int lcv=0; lcv!=temp_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*dtemp_vec)(0))[lcv])/(atol_ + rtol_*abs((*(*temp_vec)(0))[lcv]));
    enorm_val = std::max<double>(enorm_val, tmp);
  }

  Teuchos::RCP<const Epetra_MultiVector> ftemp_vec = u->data()->ViewComponent("face", false);
  Teuchos::RCP<const Epetra_MultiVector> fdtemp_vec = du->data()->ViewComponent("face", false);

  for (int lcv=0; lcv!=ftemp_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*fdtemp_vec)(0))[lcv])/(atol_ + rtol_*abs((*(*ftemp_vec)(0))[lcv]));
    enorm_val = std::max<double>(enorm_val, tmp);
  }

#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  return enorm_val;
};

// updates the preconditioner
void AdvectionDiffusion::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  S_next_->set_time(t);
  PK::solution_to_state(up, S_next_); // not sure why this isn't getting found? --etc

  // div K_e grad u
  Teuchos::RCP<CompositeVector> thermal_conductivity =
    S_next_->GetFieldData("thermal_conductivity", "energy");

  for (int c=0; c != Ke_.size(); ++c) {
    Ke_[c](0,0) = (*thermal_conductivity)("cell", c);
  }

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  preconditioner_->CreateMFDstiffnessMatrices(Ke_);
  preconditioner_->CreateMFDrhsVectors();

  // update with accumulation terms
  Teuchos::RCP<const CompositeVector> temp1 =
    S_next_->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> temp0 =
    S_inter_->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> cell_volume =
    S_next_->GetFieldData("cell_volume");

  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();
  int ncells = S_next_->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    double factor = (*cell_volume)(c);
    Acc_cells[c] += factor/h;
    Fc_cells[c] += factor/h * (*temp0)("cell",0,c);
  }

  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  preconditioner_->AssembleGlobalMatrices();

  preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
  preconditioner_->UpdateMLPreconditioner();
};

} // namespace Energy
} // namespace Amanzi
