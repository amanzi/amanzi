/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "Epetra_Vector.h"
#include "air_water_rock.hh"

namespace Amanzi {
namespace Energy {

// AirWaterRock is a BDFFnBase
// computes the non-linear functional f = f(t,u,udot)
void AirWaterRock::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                 Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {
  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(u_old, S_inter_);
  solution_to_state(u_new, S_next_);
  UpdateSecondaryVariables_(S_inter_);
  UpdateSecondaryVariables_(S_next_);

  bc_temperature_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // get access to the solution
  Teuchos::RCP<CompositeVector> res = f->data();
  res->PutScalar(0.0);

  std::cout << "residual before: " << (*res)("cell",0,0) << " " << (*res)("face",0,0) << std::endl;

  // conduction term, implicit
  ApplyConduction_(S_next_, res);
  std::cout << "residual after diffusion: " << (*res)("cell",0,0) << " " << (*res)("face",0,0) << std::endl;

  // source term?

  // accumulation term
  AddAccumulation_(res);

  // advection term, explicit
  //advection_->set_bcs(bcs_advection_, bcs_advection_dofs_);
  //AddAdvection_(S_inter_, res, true);

  std::cout << "residual after: " << (*res)("cell",0,0) << " " << (*res)("face",0,0) << std::endl;
};

// applies preconditioner to u and returns the result in Pu
void AirWaterRock::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  std::cout << "before precon: " << (*u->data())("cell",0,0) << " " << (*u->data())("face",0,0) << std::endl;
  //  preconditioner_->ApplyInverse(*u->data(), Pu->data());

  Teuchos::RCP<const CompositeVector> udat = u->data();
  Teuchos::RCP<CompositeVector> pudat = Pu->data();

  Teuchos::RCP<const Epetra_Vector> Acc = preconditioner_->Acc();

  int ncells = S_next_->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    (*pudat)("cell",0,c) = (*udat)("cell",0,c)/(*Acc)[c];
  }

  std::cout << "after precon: " << (*Pu->data())("cell",0,0) << " " << (*Pu->data())("face",0,0) << std::endl;
};

// computes a norm on u-du and returns the result
double AirWaterRock::enorm(Teuchos::RCP<const TreeVector> u,
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
void AirWaterRock::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
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

  preconditioner_->CreateMFDstiffnessMatrices(Ke_, *thermal_conductivity);
  preconditioner_->CreateMFDrhsVectors();

  // update with accumulation terms
  Teuchos::RCP<const CompositeVector> temp =
    S_next_->GetFieldData("temperature");

  Teuchos::RCP<const CompositeVector> poro =
    S_next_->GetFieldData("porosity");

  Teuchos::RCP<const CompositeVector> dens_gas;
  if (internal_energy_gas_model_->IsMolarBasis()) {
    dens_gas = S_next_->GetFieldData("molar_density_gas");
  } else {
    dens_gas = S_next_->GetFieldData("density_gas");
  }
  Teuchos::RCP<const CompositeVector> mol_frac_gas =
    S_next_->GetFieldData("mol_frac_gas");


  Teuchos::RCP<const CompositeVector> dens_liq;
  if (internal_energy_liquid_model_->IsMolarBasis()) {
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
  int ncells = S_next_->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    double T = (*temp)("cell",0,c);
    double phi = (*poro)("cell",0,c);
    double du_liq_dT = internal_energy_liquid_model_->DInternalEnergyDT(T);
    double du_gas_dT = internal_energy_gas_model_->DInternalEnergyDT(T, (*mol_frac_gas)("cell",0,c));
    double du_rock_dT = internal_energy_rock_model_->DInternalEnergyDT(T);

    double factor_gas = (*dens_gas)(c)*(*sat_gas)(c)*du_gas_dT;
    double factor_liq = (*dens_liq)(c)*(*sat_liq)(c)*du_liq_dT;
    double factor_rock = (*dens_rock)*du_rock_dT;

    double factor = (phi * (factor_gas + factor_liq) +
                     (1-phi) * factor_rock) * (*cell_volume)(c);

    Acc_cells[c] += factor/h;
  }

  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  preconditioner_->AssembleGlobalMatrices();

  preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
  preconditioner_->UpdateMLPreconditioner();

};

} // namespace Energy
} // namespace Amanzi
