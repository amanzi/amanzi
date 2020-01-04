/*
  This is the multiphase flow component of the Amanzi code. 
  This routine implements the interface functions to the FnTimeIntegratorPK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
*/

// include stuff
#include "CompW_PK.hh"
#include "Epetra_Vector.h"
#include "Op.hh"

namespace Amanzi {
namespace Multiphase {

class Op;

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
This is basically the residual
****************************************************************** */
void CompW_PK::FunctionalResidual(double t_old, double t_new, 
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new, 
                                  Teuchos::RCP<TreeVector> f)
{
  const std::vector<int>& bc_model_p = op_bc_p_->bc_model();
  const std::vector<double>& bc_value_p = op_bc_p_->bc_value();

  const std::vector<int>& bc_model_s = op_bc_s_->bc_model();
  const std::vector<double>& bc_value_s = op_bc_s_->bc_value();

  double dTp(t_new - t_old);

  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_w = u_new->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w = u_new->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl = u_new->SubVector(2)->Data();

  // solution at previous time step
  Teuchos::RCP<const CompositeVector> pressure_w_old = u_old->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w_old = u_old->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl_old = u_old->SubVector(2)->Data(); 

  // Calculate total mobility needed to initialize diffusion operator
  rel_perm_w_->Compute(*saturation_w, bc_model_p, bc_value_p);
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);

  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *rel_perm_w_->Krel());
  rel_perm_w_->Krel()->Scale(rho_/mu_);

  // compute the upwinded coefficient for the molecular diffusion operator
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);
  cvs.AddComponent("dirichlet_faces", AmanziMesh::BOUNDARY_FACE, 1);
  Teuchos::RCP<CompositeVector> s_with_face = Teuchos::rcp(new CompositeVector(cvs));
  *s_with_face->ViewComponent("cell") = *saturation_w->ViewComponent("cell");
  DeriveFaceValuesFromCellValues(*s_with_face->ViewComponent("cell"), *s_with_face->ViewComponent("face"),
                                 bc_model_s, bc_value_s);

  upwind_vw_->Scale(-1.0);
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *s_with_face);
  s_with_face->Scale(-phi_);
  upwind_vw_->Scale(-1.0);

  // compute the water component density for gravity term using the density of hydrogen in liquid
  // from the previous time step.
  // rho_w_ = rho_w_std + rho_h_l
  rho_w_->PutScalar(rho_);
  rho_w_->Update(1.0, *rhl_old, 1.0);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_matrix_->SetDensity(rho_w_);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_matrix_->ApplyBCs(true, true, true);

  op2_matrix_->global_operator()->Init();
  op2_matrix_->Setup(D1ptr, s_with_face, Teuchos::null);
  op2_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op2_matrix_->ApplyBCs(true, true, true);

  // Update source term
  // Teuchos::RCP<CompositeVector> rhs = ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->rhs();
  // if (src_sink_ != NULL) AddSourceTerms(*rhs);
  
  Teuchos::RCP<CompositeVector> f1 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f2 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->global_operator()->ComputeNegativeResidual(*pressure_w, *f1);
  op2_matrix_->global_operator()->ComputeNegativeResidual(*rhl, *f2);
  f->Data()->Update(1.0, *f1, 1.0, *f2, 0.0);
  
  // Add time derivative
  Epetra_MultiVector& f_cell = *f->Data()->ViewComponent("cell", true);

  const Epetra_MultiVector& S_new_c = *u_new->SubVector(1)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& S_old_c = *u_old->SubVector(1)->Data()->ViewComponent("cell");

  // Add accumulation term
  double s0, s1, volume;
  for (int c = 0; c < f_cell.MyLength(); c++) {
    s1 = S_new_c[0][c];
    s0 = S_old_c[0][c];

    double factor = phi_ * mesh_->cell_volume(c) / dTp;
    double tmp_acc_term = rho_ * (s1 - s0) * factor;
    f_cell[0][c] += rho_ * (s1 - s0) * factor;
  }
  f->Scale(dTp);
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int CompW_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                  Teuchos::RCP<TreeVector> Pu) {
  return 0;
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void CompW_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  const std::vector<int>& bc_model_p = op_bc_p_->bc_model();
  const std::vector<double>& bc_value_p = op_bc_p_->bc_value();

  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_w = u->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w = u->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl = u->SubVector(2)->Data();

  // Calculate relative perm needed to initialize diffusion operator
  rel_perm_w_->Compute(*saturation_w, bc_model_p, bc_value_p);
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);

  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *rel_perm_w_->Krel());
  rel_perm_w_->Krel()->Scale(dTp*rho_/mu_);

  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *rel_perm_w_->dKdS());
  rel_perm_w_->dKdS()->Scale(dTp*rho_ / mu_); 

  // the flux is only used for finding the upwind cells
  Teuchos::RCP<CompositeVector> tmp_flux_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_wet")));
  tmp_flux_->PutScalar(0.0);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->Setup(Kptr, Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateFlux(pressure_w.ptr(), tmp_flux_.ptr());

  // rel_perm_w_->dKdS() is used as the coefficients of the matrix
  // Epetra_MultiVector& krel_f = *rel_perm_w_->Krel()->ViewComponent("face");
  Epetra_MultiVector& fface = *tmp_flux_->ViewComponent("face", true);
  Epetra_MultiVector& dKdS_f = *rel_perm_w_->dKdS()->ViewComponent("face", true);
  AMANZI_ASSERT(dKdS_f.MyLength() == fface.MyLength());

  for (int f = 0; f < dKdS_f.MyLength(); f++) {
    dKdS_f[0][f] *= fface[0][f]; 
    //if (fabs(krel_f[0][f]) < 1e-12) {
    //  dKdS_f[0][f] = 0.0;
    //} else {
    //  dKdS_f[0][f] *= fface[0][f]/krel_f[0][f]; 
    //}
  }

  // A_11 block wrt Pw
  op1_preconditioner_->global_operator()->Init();
  op1_preconditioner_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_preconditioner_->ApplyBCs(true, true, true);
  // op1_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  // op1_preconditioner_->global_operator()->AssembleMatrix();

  // A_12 block wrt Sw
  //tmp_flux_->Scale(-1.0);
  op2_preconditioner_->global_operator()->Init();
  op2_preconditioner_->Setup(*tmp_flux_);
  op2_preconditioner_->UpdateMatrices(rel_perm_w_->dKdS().ptr());
  op2_preconditioner_->ApplyBCs(true, true, true);
  op2_preconditioner_->global_operator()->Rescale(-1.0);

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);
  Teuchos::RCP<CompositeVector> s_with_face = Teuchos::rcp(new CompositeVector(cvs));
  /*
  *s_with_face->ViewComponent("cell") = *saturation_w->ViewComponent("cell");
  //DeriveFaceValuesFromCellValues(*s_with_face->ViewComponent("cell"), *s_with_face->ViewComponent("face"),
  //                               bc_model_, bc_value_s_);
  CoefUpwindFn1 func3 = &MPCoeff::ValuePrimaryVar;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_s_,
                     *s_with_face, *s_with_face, func3);
  */
  s_with_face->PutScalar(-dTp*phi_);

  tmp_flux_->PutScalar(0.0);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  op2_matrix_->global_operator()->Init();
  op2_matrix_->Setup(D1ptr, s_with_face, Teuchos::null);
  op2_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op2_matrix_->UpdateFlux(rhl.ptr(), tmp_flux_.ptr());

  op_prec_sat_->Setup(*tmp_flux_);
  tmp_flux_->Scale(-1.0);
  op_prec_sat_->UpdateMatrices(tmp_flux_.ptr());
  op_prec_sat_->ApplyBCs(true, true, true);
  op_prec_sat_->global_operator()->Rescale(-1.0);

  // create operator for accumulation term
  //op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_prec_sat_->global_operator()));

  // Create cvs for accumulation term 
  /*
  CompositeVectorSpace cvs1; 
  cvs1.SetMesh(mesh_);
  cvs1.SetGhosted(true);
  cvs1.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs1.SetOwned(false);
  */
  
  CompositeVector porosity(*S_->GetFieldData("pressure_w"));
  porosity.PutScalar(phi_);
  porosity.Scale(rho_);
  if (dTp > 0.0) {
    op_acc_->AddAccumulationDelta(*saturation_w, porosity, porosity, 1.0, "cell");
  } 
  // op_prec_sat_->global_operator()->SymbolicAssembleMatrix();
  // op_prec_sat_->global_operator()->AssembleMatrix();

  // A_13 block wrt rhl
  // s_with_face->Scale(-1.0);
  op3_preconditioner_->global_operator()->Init();
  op3_preconditioner_->Setup(D1ptr, s_with_face, Teuchos::null);
  op3_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op3_preconditioner_->ApplyBCs(true, true, true);
  // op3_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  // op3_preconditioner_->global_operator()->AssembleMatrix();
}

}  // namespace Multiphase
}  // namespace Amanzi
