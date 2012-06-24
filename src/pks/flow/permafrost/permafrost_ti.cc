/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

// Permafrost is a BDFFnBase
// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void Permafrost::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  std::cout << "Precon update at t = " << t << std::endl;
  // update state with the solution up.
  S_next_->set_time(t);
  PK::solution_to_state(up, S_next_);
  UpdateSecondaryVariables_(S_next_);

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_);
  Teuchos::RCP<const CompositeVector> num_rel_perm =
    S_next_->GetFieldData("numerical_rel_perm");

  preconditioner_->CreateMFDstiffnessMatrices(*num_rel_perm);
  preconditioner_->CreateMFDrhsVectors();
  AddGravityFluxes_(S_next_, preconditioner_);

  // update with accumulation terms
  Teuchos::RCP<const CompositeVector> pres = S_next_->GetFieldData("pressure");
  Teuchos::RCP<const double> p_atm = S_next_->GetScalarData("atmospheric_pressure");
  Teuchos::RCP<const CompositeVector> temp = S_next_->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> poro = S_next_->GetFieldData("porosity");

  Teuchos::RCP<const CompositeVector> mol_frac_gas =
    S_next_->GetFieldData("mol_frac_gas");
  Teuchos::RCP<const CompositeVector> n_gas = S_next_->GetFieldData("molar_density_gas");
  Teuchos::RCP<const CompositeVector> sat_gas = S_next_->GetFieldData("saturation_gas");

  Teuchos::RCP<const CompositeVector> n_liq = S_next_->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> sat_liq = S_next_->GetFieldData("saturation_liquid");

  Teuchos::RCP<const CompositeVector> n_ice = S_next_->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> sat_ice = S_next_->GetFieldData("saturation_ice");

  Teuchos::RCP<const CompositeVector> cell_volume = S_next_->GetFieldData("cell_volume");

  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();

  Teuchos::RCP<CompositeVector> sat_star_gl = Teuchos::rcp(new CompositeVector(*sat_liq));
  UnfrozenSaturation_(S_next_, *pres, *p_atm, sat_star_gl);
  Teuchos::RCP<CompositeVector> dsat_star_gl = Teuchos::rcp(new CompositeVector(*sat_liq));
  DUnfrozenSaturationDp_(S_next_, *pres, *p_atm, dsat_star_gl);

  int ncells = pres->size("cell");
  for (int c=0; c!=ncells; ++c) {
    // many terms here...
    // accumulation term is d/dt ( phi * (omega_g*n_g*s_g + n_l*s_l) )
    // note: s_g = (1 - s_l)
    // note: mol_frac of vapor is not a function of pressure
    // note: assumes phi does not depend on pressure
    double p = (*pres)("cell",c);
    double T = (*temp)("cell",c);
    double phi = (*poro)("cell",c);

    if (c==0) std::cout << "    p =" << p << std::endl;
    if (c==0) std::cout << "    T =" << T << std::endl;
    if (c==0) std::cout << "    phi =" << phi << std::endl;
    if (c==99) std::cout << "    p =" << p << std::endl;
    if (c==99) std::cout << "    T =" << T << std::endl;
    if (c==99) std::cout << "    phi =" << phi << std::endl;


    //  omega_g * sat_g * d(n_g)/dp
    double result = (*mol_frac_gas)("cell",c) * (*sat_gas)("cell",c) * eos_gas_->DMolarDensityDp(T,p);
    if (c==0) std::cout << "    res0 (0) =" << result << std::endl;
    if (c==99) std::cout << "    res0 (99) =" << result << std::endl;

    // + sat_g * n_g * d(omega_g)/dp  (mol frac vapor not a function of pressure)

    // + sat_l * d(n_l)/dp
    result += (*sat_liq)("cell",c) * eos_liquid_->DMolarDensityDp(T,p);
    if (c==0) std::cout << "    res1 (0) =" << result << std::endl;
    if (c==99) std::cout << "    res1 (99) =" << result << std::endl;

    // + sat_i * d(n_i)/dp
    result += (*sat_ice)("cell",c) * eos_ice_->DMolarDensityDp(T,p);
    if (c==0) std::cout << "    res2 (0) =" << result << std::endl;
    if (c==99) std::cout << "    res2 (99) =" << result << std::endl;

    // + n_i  * d(sat_i)/dp + n_l * d(sat_l)/dp + omega_g * n_g * d(sat_g)/dp
    // not obvious... requires derivation from A,B, etc, and specific to dA/dp = 0
    // and B = 1/S*(p)
    result += (*sat_liq)("cell",c) * ( (*mol_frac_gas)("cell",c)*(*n_gas)("cell",c) *
            (1.0 - (*sat_gas)("cell",c))
            - (*n_liq)("cell",c)*(*sat_liq)("cell",c)
            - (*n_ice)("cell",c)*(*sat_ice)("cell",c)) *
                -(*dsat_star_gl)("cell",c) / (*sat_star_gl)("cell",c) /
                                                (*sat_star_gl)("cell",c);

    if (c==0) std::cout << "    res3 (0) =" << result << std::endl;
    if (c==99) std::cout << "    res3 (99) =" << result << std::endl;

    Acc_cells[c] += phi * result * (*cell_volume)("cell",c) / h;
    Fc_cells[c] += phi * result * (*cell_volume)("cell",c) / h * p;
  }

  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  preconditioner_->AssembleGlobalMatrices();
  preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);

  // Code to dump Schur complement to check condition number
  /*
  Teuchos::RCP<Epetra_FECrsMatrix> sc = preconditioner_->Schur();
  std::stringstream filename_s;
  filename_s << "schur_" << S_next_->cycle() << ".txt";
  //a  std::string filename = filename_s.str();
  EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
  std::cout << "updated precon " << S_next_->cycle() << std::endl;
  */

  preconditioner_->UpdateMLPreconditioner();

  // test_precon(t, up, h);
};

// Runs a very expensive FD test of the Jacobian and prints out an enorm
//  measure of the error.
void Permafrost::test_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::RCP<TreeVector> dp = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> f2 = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> df = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> uold = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> unew = Teuchos::rcp(new TreeVector(*up));

  double maxval = 0.0;

  int ncells = up->data()->size("cell");
  for (int c=0; c!=ncells; ++c) {
    *unew = (*up);
    fun(t-h, t, uold, unew, f1);

    dp->PutScalar(0.0);
    (*dp->data())("cell",c) = 0.0001;
    unew->Update(1.0, *dp, 1.0);
    fun(t-h, t, uold, unew, f2);

    preconditioner_->Apply(*dp->data(), df->data());
    df->Update(-1.0, *f2, 1.0, *f1, 1.0);
    maxval = std::max(maxval, enorm(f1, df));
  }
  std::cout << "Testing PC with FD.  Error: " << maxval << std::endl;
};

}  // namespace Flow
}  // namespace Amanzi



