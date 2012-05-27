/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "Epetra_Vector.h"
#include "two_phase.hh"

namespace Amanzi {
namespace Energy {

// TwoPhase is a BDFFnBase
// computes the non-linear functional g = g(t,u,udot)
void TwoPhase::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  Teuchos::RCP<CompositeVector> u = u_new->data();
  std::cout << "Two-Phase Residual calculation:" << std::endl;
  std::cout << "  T0: " << (*u)("cell",0,0) << " " << (*u)("face",0,3) << std::endl;
  std::cout << "  T1: " << (*u)("cell",0,99) << " " << (*u)("face",0,497) << std::endl;

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(u_new, S_next_);
  UpdateSecondaryVariables_(S_next_);

  // update boundary conditions
  bc_temperature_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  ApplyDiffusion_(S_next_, res);
  std::cout << "  res0 (after diffusion): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after diffusion): " << (*res)("cell",0,99) << " " << (*res)("face",0,497) << std::endl;

  // accumulation term
  AddAccumulation_(res);
  std::cout << "  res0 (after accumulation): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after accumulation): " << (*res)("cell",0,99) << " " << (*res)("face",0,497) << std::endl;

  // advection term, explicit
  AddAdvection_(S_inter_, res, true);
  std::cout << "  res0 (after advection): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after advection): " << (*res)("cell",0,99) << " " << (*res)("face",0,497) << std::endl;
};

// applies preconditioner to u and returns the result in Pu
void TwoPhase::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  std::cout << "Precon application:" << std::endl;
  std::cout << "  T0: " << (*u->data())("cell",0,0) << " " << (*u->data())("face",0,3) << std::endl;
  std::cout << "  T1: " << (*u->data())("cell",0,99) << " " << (*u->data())("face",0,497) << std::endl;
  preconditioner_->ApplyInverse(*u->data(), Pu->data());
  std::cout << "  PC*T0: " << (*Pu->data())("cell",0,0) << " " << (*Pu->data())("face",0,3) << std::endl;
  std::cout << "  PC*T1: " << (*Pu->data())("cell",0,99) << " " << (*Pu->data())("face",0,497) << std::endl;
};

// computes a norm on u-du and returns the result
double TwoPhase::enorm(Teuchos::RCP<const TreeVector> u,
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
void TwoPhase::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  S_next_->set_time(t);
  PK::solution_to_state(up, S_next_); // not sure why this isn't getting found? --etc
  UpdateSecondaryVariables_(S_next_);
  UpdateThermalConductivity_(S_next_);

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
  Teuchos::RCP<const CompositeVector> temp =
    S_next_->GetFieldData("temperature");

  Teuchos::RCP<const CompositeVector> poro =
    S_next_->GetFieldData("porosity");

  Teuchos::RCP<const CompositeVector> dens_gas;
  if (iem_gas_->IsMolarBasis()) {
    dens_gas = S_next_->GetFieldData("molar_density_gas");
  } else {
    dens_gas = S_next_->GetFieldData("density_gas");
  }
  Teuchos::RCP<const CompositeVector> mol_frac_gas =
    S_next_->GetFieldData("mol_frac_gas");


  Teuchos::RCP<const CompositeVector> dens_liq;
  if (iem_liquid_->IsMolarBasis()) {
    dens_liq = S_next_->GetFieldData("molar_density_liquid");
  } else {
    dens_liq = S_next_->GetFieldData("density_liquid");
  }

  Teuchos::RCP<const CompositeVector> sat_liq =
    S_next_->GetFieldData("saturation_liquid");
  Teuchos::RCP<const CompositeVector> sat_gas =
    S_next_->GetFieldData("saturation_gas");

  Teuchos::RCP<const CompositeVector> cell_volume =
    S_next_->GetFieldData("cell_volume");

  Teuchos::RCP<const double> dens_rock =
    S_next_->GetScalarData("density_rock");

  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();
  int ncells = S_next_->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    // accumulation term is d/dt ( phi * (s_g*n_g*u_g + s_l*n_l*u_l) + (1-phi)*rho_r*u_r
    // note: s_g, s_l do not depend upon T
    // note: CURRENTLY IGNORING the density dependence in temperature in this
    //       preconditioner.  This is idiodic, but due to poor design on my
    //       part we don't have access to the equations of state here.  This
    //       must be fixed, but need to figure out a new paradigm for this
    //       sort of thing. -- etc
    // note: also assumes phi does not depend on temperature.
    double T = (*temp)("cell",0,c);
    double phi = (*poro)("cell",0,c);
    double du_liq_dT = iem_liquid_->DInternalEnergyDT(T);
    double du_gas_dT = iem_gas_->DInternalEnergyDT(T, (*mol_frac_gas)("cell",0,c));
    double du_rock_dT = iem_rock_->DInternalEnergyDT(T);

    double factor_gas = (*dens_gas)(c)*(*sat_gas)(c)*du_gas_dT;
    double factor_liq = (*dens_liq)(c)*(*sat_liq)(c)*du_liq_dT;
    double factor_rock = (*dens_rock)*du_rock_dT;

    double factor = (phi * (factor_gas + factor_liq) +
                     (1-phi) * factor_rock) * (*cell_volume)(c);

    Acc_cells[c] += factor/h;
    Fc_cells[c] += factor/h * T;
  }

  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  preconditioner_->AssembleGlobalMatrices();

  preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
  preconditioner_->UpdateMLPreconditioner();

};

} // namespace Energy
} // namespace Amanzi
