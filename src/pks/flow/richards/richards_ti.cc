/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 0

// Richards is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void Richards::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  Teuchos::RCP<CompositeVector> u = u_new->data();
#if DEBUG_FLAG
  std::cout << "Richards Residual calculation:" << std::endl;
  std::cout << "  p0: " << (*u)("cell",0,0) << " " << (*u)("face",0,3) << std::endl;
  std::cout << "  p1: " << (*u)("cell",0,99) << " " << (*u)("face",0,497) << std::endl;
#endif

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_, res);
#if DEBUG_FLAG
  std::cout << "  res0 (after diffusion): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after diffusion): " << (*res)("cell",0,99) << " " << (*res)("face",0,497) << std::endl;
#endif

  // accumulation term
  AddAccumulation_(res);
#if DEBUG_FLAG
  std::cout << "  res0 (after accumulation): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after accumulation): " << (*res)("cell",0,99) << " " << (*res)("face",0,497) << std::endl;
#endif
};

// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void Richards::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
#if DEBUG_FLAG
  std::cout << "Precon application:" << std::endl;
  std::cout << "  p0: " << (*u->data())("cell",0,0) << " " << (*u->data())("face",0,3) << std::endl;
  std::cout << "  p1: " << (*u->data())("cell",0,99) << " " << (*u->data())("face",0,497) << std::endl;
#endif
  preconditioner_->ApplyInverse(*u->data(), Pu->data());
#if DEBUG_FLAG
  std::cout << "  PC*p0: " << (*Pu->data())("cell",0,0) << " " << (*Pu->data())("face",0,3) << std::endl;
  std::cout << "  PC*p1: " << (*Pu->data())("cell",0,99) << " " << (*Pu->data())("face",0,497) << std::endl;
#endif
};


// -----------------------------------------------------------------------------
// Compute a norm on (u,du)
// -----------------------------------------------------------------------------
double Richards::enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) {
  Teuchos::RCP<const CompositeVector> pres = u->data();
  Teuchos::RCP<const CompositeVector> dpres = du->data();


  double enorm_val_cell = 0.0;
  for (int c=0; c!=pres->size("cell",false); ++c) {
    if (boost::math::isnan<double>((*dpres)("cell",c))) {
      std::cout << "Cutting time step due to NaN in correction." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }

    double tmp = abs((*dpres)("cell",c)) / (atol_ + rtol_ * abs((*pres)("cell",c)));
    enorm_val_cell = std::max<double>(enorm_val_cell, tmp);
    //    printf("cell: %5i %14.7e %14.7e\n",lcv,(*(*dpres_vec)(0))[lcv],tmp);
  }

  double enorm_val_face = 0.0;
  for (int f=0; f!=pres->size("face",false); ++f) {
    if (boost::math::isnan<double>((*dpres)("face",f))) {
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }

    double tmp = abs((*dpres)("face",f)) / (atol_ + rtol_ * abs((*pres)("face",f)));
    enorm_val_face = std::max<double>(enorm_val_face, tmp);
    //    printf("face: %5i %14.7e %14.7e\n",lcv,(*(*fdpres_vec)(0))[lcv],tmp);
  }


  //  std::cout.precision(15);
  //  std::cout << "enorm val (cell, face): " << std::scientific << enorm_val_cell
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
void Richards::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  std::cout << "Precon update at t = " << t << std::endl;
  // update state with the solution up.
  S_next_->set_time(t);
  PK::solution_to_state(up, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_);
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_next_->GetFieldData("numerical_rel_perm");

  preconditioner_->CreateMFDstiffnessMatrices(*rel_perm);
  preconditioner_->CreateMFDrhsVectors();
  AddGravityFluxes_(S_next_, preconditioner_);

  // update with accumulation terms
  // -- update the accumulation derivatives
  S_next_->GetFieldModel("water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), "richards_pk", "pressure");

  // -- get the accumulation deriv
  Teuchos::RCP<const CompositeVector> dwc_dp =
      S_next_->GetFieldData("dwater_content_dpressure");
  Teuchos::RCP<const CompositeVector> pres =
      S_next_->GetFieldData("pressure");

  // -- get the matrices/rhs that need updating
  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();

  int ncells = pres->size("cell");
  for (int c=0; c!=ncells; ++c) {
    Acc_cells[c] += (*dwc_dp)("cell",c) / h;
    Fc_cells[c] += (*dwc_dp)("cell",c) / h * (*pres)("cell",c);
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

  preconditioner_->UpdatePreconditioner();
};

}  // namespace Flow
}  // namespace Amanzi



