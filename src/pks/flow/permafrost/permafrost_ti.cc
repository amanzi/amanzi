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

/* ******************************************************************
 * Calculate f(u, du/dt) = d(s u)/dt + A*u - g.
 ****************************************************************** */
void Permafrost::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  Teuchos::RCP<CompositeVector> u = u_new->data();
  std::cout << "Permafrost Residual calculation:" << std::endl;
  std::cout << "  p0: " << (*u)("cell",0,0) << " " << (*u)("face",0,3) << std::endl;
  std::cout << "  p1: " << (*u)("cell",0,99) << " " << (*u)("face",0,497) << std::endl;

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);
  UpdateSecondaryVariables_(S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  //bc_head_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_, res);
  std::cout << "  res0 (after diffusion): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after diffusion): " << (*res)("cell",0,99) << " " << (*res)("face",0,497) << std::endl;

  // accumulation term
  AddAccumulation_(res);
  std::cout << "  res0 (after accumulation): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after accumulation): " << (*res)("cell",0,99) << " " << (*res)("face",0,497) << std::endl;
};

/* ******************************************************************
* Apply preconditioner to u and return the result in Pu.
****************************************************************** */
void Permafrost::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  std::cout << "Precon application:" << std::endl;
  std::cout << "  p0: " << (*u->data())("cell",0,0) << " " << (*u->data())("face",0,3) << std::endl;
  std::cout << "  p1: " << (*u->data())("cell",0,99) << " " << (*u->data())("face",0,497) << std::endl;
  preconditioner_->ApplyInverse(*u->data(), Pu->data());
  std::cout << "  PC*p0: " << (*Pu->data())("cell",0,0) << " " << (*Pu->data())("face",0,3) << std::endl;
  std::cout << "  PC*p1: " << (*Pu->data())("cell",0,99) << " " << (*Pu->data())("face",0,497) << std::endl;
};


/* ******************************************************************
 * computes a norm on u-du and returns the result
 ****************************************************************** */
double Permafrost::enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) {
  double enorm_val = 0.0;
  Teuchos::RCP<const Epetra_MultiVector> pres_vec = u->data()->ViewComponent("cell", false);
  Teuchos::RCP<const Epetra_MultiVector> dpres_vec = du->data()->ViewComponent("cell", false);

  for (int lcv=0; lcv!=pres_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*dpres_vec)(0))[lcv])/(atol_ + rtol_*abs((*(*pres_vec)(0))[lcv]));
    if (std::isnan(tmp)) return 1.e99;
    enorm_val = std::max<double>(enorm_val, tmp);
  }

  Teuchos::RCP<const Epetra_MultiVector> fpres_vec = u->data()->ViewComponent("face", false);
  Teuchos::RCP<const Epetra_MultiVector> fdpres_vec = du->data()->ViewComponent("face", false);

  for (int lcv=0; lcv!=fpres_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*fdpres_vec)(0))[lcv])/(atol_ + rtol_*abs((*(*fpres_vec)(0))[lcv]));
    if (std::isnan(tmp)) return 1.e99;
    enorm_val = std::max<double>(enorm_val, tmp);
  }

#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  return enorm_val;
};

/* ******************************************************************
* Compute new preconditioner B(p, dT_prec).
****************************************************************** */
void Permafrost::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  std::cout << "Precon update at t = " << t << std::endl;
  // update state with the solution up.
  S_next_->set_time(t);
  PK::solution_to_state(up, S_next_);
  UpdateSecondaryVariables_(S_next_);

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  //bc_head_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_);
  Teuchos::RCP<const CompositeVector> rel_perm_faces =
    S_next_->GetFieldData("rel_perm_faces");

  preconditioner_->CreateMFDstiffnessMatrices(K_, rel_perm_faces);
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
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells = preconditioner_->Aff_cells();
  std::vector<Epetra_SerialDenseVector>& Acf_cells = preconditioner_->Acf_cells();
  std::vector<Epetra_SerialDenseVector>& Afc_cells = preconditioner_->Afc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();
  std::vector<Epetra_SerialDenseVector>& Ff_cells = preconditioner_->Ff_cells();

  Teuchos::RCP<CompositeVector> sat_star_gl = Teuchos::rcp(new CompositeVector(*sat_liq));
  UnfrozenSaturation_(S_next_, *pres, *p_atm, sat_star_gl);
  Teuchos::RCP<CompositeVector> dsat_star_gl = Teuchos::rcp(new CompositeVector(*sat_liq));
  DUnfrozenSaturationDp_(S_next_, *pres, *p_atm, dsat_star_gl);

  int ncells = S_next_->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    // many terms here...
    // accumulation term is d/dt ( phi * (omega_g*n_g*s_g + n_l*s_l) )
    // note: s_g = (1 - s_l)
    // note: mol_frac of vapor is not a function of pressure
    // note: assumes phi does not depend on pressure
    double p = (*pres)("cell",0,c);
    double T = (*temp)("cell",0,c);
    double phi = (*poro)("cell",0,c);

    if (c==0) std::cout << "    p =" << p << std::endl;
    if (c==0) std::cout << "    T =" << T << std::endl;
    if (c==0) std::cout << "    phi =" << phi << std::endl;
    if (c==99) std::cout << "    p =" << p << std::endl;
    if (c==99) std::cout << "    T =" << T << std::endl;
    if (c==99) std::cout << "    phi =" << phi << std::endl;


    //  omega_g * sat_g * d(n_g)/dp
    double result = (*mol_frac_gas)(c) * (*sat_gas)(c) * eos_gas_->DMolarDensityDp(T,p);
    if (c==0) std::cout << "    res0 (0) =" << result << std::endl;
    if (c==99) std::cout << "    res0 (99) =" << result << std::endl;

    // + sat_g * n_g * d(omega_g)/dp  (mol frac vapor not a function of pressure)

    // + sat_l * d(n_l)/dp
    result += (*sat_liq)(c) * eos_liquid_->DMolarDensityDp(T,p);
    if (c==0) std::cout << "    res1 (0) =" << result << std::endl;
    if (c==99) std::cout << "    res1 (99) =" << result << std::endl;

    // + sat_i * d(n_i)/dp
    result += (*sat_ice)(c) * eos_ice_->DMolarDensityDp(T,p);
    if (c==0) std::cout << "    res2 (0) =" << result << std::endl;
    if (c==99) std::cout << "    res2 (99) =" << result << std::endl;

    // + n_i  * d(sat_i)/dp + n_l * d(sat_l)/dp + omega_g * n_g * d(sat_g)/dp
    // not obvious... requires derivation from A,B, etc, and specific to dA/dp = 0
    // and B = 1/S*(p)
    result += (*sat_liq)(c) * ( (*mol_frac_gas)(c)*(*n_gas)(c)*(1.0 - (*sat_gas)(c))
                               - (*n_liq)(c)*(*sat_liq)(c) - (*n_ice)(c)*(*sat_ice)(c))
               * -(*dsat_star_gl)(c)/(*sat_star_gl)(c)/(*sat_star_gl)(c);

    if (c==0) std::cout << "    res3 (0) =" << result << std::endl;
    if (c==99) std::cout << "    res3 (99) =" << result << std::endl;

    Acc_cells[c] += phi * result * (*cell_volume)(c) / h;
    Fc_cells[c] += phi * result * (*cell_volume)(c) / h * p;
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

  int ncells = S_next_->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    *unew = (*up);
    fun(t-h, t, uold, unew, f1);

    dp->PutScalar(0.0);
    (*dp->data())("cell",0,c) = 0.0001;
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



