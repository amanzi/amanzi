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
#include "CompH_PK.hh"
#include "Epetra_Vector.h"
#include "Op.hh"

namespace Amanzi {
namespace Multiphase {

class Op;

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
This is basically the residual
****************************************************************** */
void CompH_PK::FunctionalResidual(double t_old, double t_new, 
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new, 
                                  Teuchos::RCP<TreeVector> f)
{
  const std::vector<int>& bc_model_p = op_bc_p_->bc_model();
  const std::vector<double>& bc_value_p = op_bc_p_->bc_value();

  const std::vector<int>& bc_model_p_n = op_bc_p_n_->bc_model();
  const std::vector<double>& bc_value_p_n = op_bc_p_n_->bc_value();

  double dTp(t_new - t_old);
  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_w = u_new->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w = u_new->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl = u_new->SubVector(2)->Data();

  // solution at previous time step
  Teuchos::RCP<const CompositeVector> pressure_w_old = u_old->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w_old = u_old->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl_old = u_old->SubVector(2)->Data(); 

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);  
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);

  // old nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n_old = Teuchos::rcp(new CompositeVector(*pressure_w_old));
  capillary_pressure_->Compute(*saturation_w_old);
  pressure_n_old->Update(1.0, *capillary_pressure_->Pc(), 1.0);

  // compute density of gas rho_g = rho_h_g = P_g * Cg
  Teuchos::RCP<CompositeVector> rho_g = Teuchos::rcp(new CompositeVector(*pressure_n));
  rho_g->Scale(Cg_);

  // compute old molar fractions
  Teuchos::RCP<CompositeVector> rho_g_old = Teuchos::rcp(new CompositeVector(*pressure_n_old));
  rho_g_old->Scale(Cg_);

  // Calculate total mobility needed to initialize diffusion operator
  coef_w_->Compute(*rhl, *saturation_w, bc_model_p, bc_value_p);
  coef_n_->Compute(*pressure_w, *saturation_w, bc_model_p_n, bc_value_p_n);
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);
  upwind_vn_ = S_->GetFieldData("velocity_nonwet", passwd_);

  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *coef_w_->Coeff());
  coef_w_->Coeff()->Scale(1.0 / mu1_);
  upwind_n_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *coef_n_->Coeff());
  coef_n_->Coeff()->Scale(1.0 / mu2_);

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
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *s_with_face);
  s_with_face->Scale(phi_);

  // compute the water component density for gravity term using the density of hydrogen in liquid
  // from the previous time step.
  // rho_w_ = rho_w_std + rho_h_l
  rho_w_->PutScalar(rho_);
  rho_w_->Update(1.0, *rhl_old, 1.0);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->Setup(Kptr, coef_w_->Coeff(), Teuchos::null);
  op1_matrix_->SetDensity(rho_w_);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->ApplyBCs(true, true, true);
  //((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->SymbolicAssembleMatrix();
  //((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->AssembleMatrix();

  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->Setup(Kptr, coef_n_->Coeff(), Teuchos::null);
  op2_matrix_->SetDensity(rho_g_old);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->ApplyBCs(true, true, true);
  //((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->SymbolicAssembleMatrix();
  //((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->AssembleMatrix(); 

  op3_matrix_->global_operator()->Init();
  op3_matrix_->Setup(D1ptr, s_with_face, Teuchos::null);
  op3_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op3_matrix_->ApplyBCs(true, true, true);
  //op3_matrix_->SymbolicAssembleMatrix();
  //op3_matrix_->AssembleMatrix(); 

  Teuchos::RCP<CompositeVector> f1 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f2 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f3 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->global_operator()->ComputeNegativeResidual(*pressure_w, *f1);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->ComputeNegativeResidual(*pressure_n, *f2);
  op3_matrix_->global_operator()->ComputeNegativeResidual(*rhl, *f3);
  f->Data()->Update(1.0, *f1, 1.0, *f2, 0.0);
  f->Data()->Update(1.0, *f3, 1.0);
  
  // Add time derivative
  Epetra_MultiVector& f_cell = *f->Data()->ViewComponent("cell", true);

  const Epetra_MultiVector& S_new_c = *u_new->SubVector(1)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& S_old_c = *u_old->SubVector(1)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& rhl_new_c = *rhl->ViewComponent("cell");
  const Epetra_MultiVector& rhl_old_c = *rhl_old->ViewComponent("cell");
  const Epetra_MultiVector& rhg_new_c = *rho_g->ViewComponent("cell");
  const Epetra_MultiVector& rhg_old_c = *rho_g_old->ViewComponent("cell");

  // Add accumulation term
  double s0, s1, rhl1, rhl0, rhg1, rhg0, volume;
  for (int c = 0; c < f_cell.MyLength(); c++) {
    s1 = S_new_c[0][c];
    s0 = S_old_c[0][c];
    rhl1 = rhl_new_c[0][c];
    rhl0 = rhl_old_c[0][c];
    rhg1 = rhg_new_c[0][c];
    rhg0 = rhg_old_c[0][c];

    double factor = phi_ * mesh_->cell_volume(c) / dTp;
    double tmp_acc_term = (rhg1*(1.0-s1) + rhl1*s1 - rhg0*(1.0-s0) - rhl0*s0) * factor;
    f_cell[0][c] += tmp_acc_term;
  }
  f->Scale(dTp);
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int CompH_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                  Teuchos::RCP<TreeVector> Pu) {
  return 0;
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void CompH_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  const std::vector<int>& bc_model_p = op_bc_p_->bc_model();
  const std::vector<double>& bc_value_p = op_bc_p_->bc_value();

  const std::vector<int>& bc_model_p_n = op_bc_p_n_->bc_model();
  const std::vector<double>& bc_value_p_n = op_bc_p_n_->bc_value();

  Teuchos::RCP<const CompositeVector> pressure_w = u->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w = u->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl = u->SubVector(2)->Data();

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);

  // compute density of gas rho_g = rho_h_g = P_g * Cg
  Teuchos::RCP<CompositeVector> rho_g = Teuchos::rcp(new CompositeVector(*pressure_n));
  rho_g->Scale(Cg_);

  // Calculate total mobility needed to initialize diffusion operator
  coef_w_->Compute(*rhl, *saturation_w, bc_model_p, bc_value_p);
  coef_n_->Compute(*pressure_w, *saturation_w, bc_model_p_n, bc_value_p_n);
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);
  upwind_vn_ = S_->GetFieldData("velocity_nonwet", passwd_);

  // compute flux for upwind direction of advection operator for nonwetting phase
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->Setup(Kptr, Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateFlux(pressure_n.ptr(), upwind_vn_.ptr());

  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *coef_w_->Coeff());
  // coef_w_->Coeff()->Scale(dTp / mu1_);
  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_p, *coef_n_->Coeff());
  // coef_n_->Coeff()->Scale(dTp / mu2_);
  Teuchos::RCP<CompositeVector> total_coef = Teuchos::rcp(new CompositeVector(*coef_w_->Coeff()));
  total_coef->Update(dTp/mu2_, *coef_n_->Coeff(), dTp/mu1_);

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);
  cvs.AddComponent("dirichlet_faces", AmanziMesh::BOUNDARY_FACE, 1);
  Teuchos::RCP<CompositeVector> s_with_face = Teuchos::rcp(new CompositeVector(cvs));
  *s_with_face->ViewComponent("cell") = *saturation_w->ViewComponent("cell");
  //DeriveFaceValuesFromCellValues(*s_with_face->ViewComponent("cell"), *s_with_face->ViewComponent("face"),
  //                               bc_model_, bc_value_s_);
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *s_with_face);
  s_with_face->Scale(dTp*phi_);

  // A_21 block wrt Pw
  op1_preconditioner_->global_operator()->Init();
  op1_preconditioner_->Setup(Kptr, total_coef, Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_preconditioner_->ApplyBCs(true, true, true);
  // op1_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  // op1_preconditioner_->global_operator()->AssembleMatrix();

  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_p, *coef_n_->Krel());
  Teuchos::RCP<CompositeVector> lambda_n = Teuchos::rcp(new CompositeVector(*coef_n_->Krel()));
  lambda_n->Scale(Cg_*dTp / mu2_);

  Teuchos::RCP<CompositeVector> tmp_flux_ = Teuchos::rcp(new CompositeVector(*upwind_vn_));
  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->Setup(Kptr, lambda_n, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateFlux(pressure_n.ptr(), tmp_flux_.ptr());

  // Teuchos::ParameterList olist_adv = comp_list_->sublist("operators").sublist("advection operator");
  // op_prec_pres_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, op1_preconditioner_->global_operator()));
  op_prec_pres_->Setup(*tmp_flux_);
  // tmp_flux_->Scale(-1.0);
  op_prec_pres_->UpdateMatrices(tmp_flux_.ptr());
  op_prec_pres_->ApplyBCs(true, true, true);

  // create operator for accumulation term
  //op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_prec_pres_->global_operator()));

  CompositeVector acc_term(*saturation_w);
  acc_term.Scale(-1.0);
  acc_term.Shift(1.0);
  acc_term.Scale(phi_*Cg_);
  if (dTp > 0.0) {
    op_acc_->AddAccumulationDelta(*saturation_w, acc_term, acc_term, 1.0, "cell");
  } 
  // op_prec_pres_->global_operator()->SymbolicAssembleMatrix();
  // op_prec_pres_->global_operator()->AssembleMatrix();

  // A_22 block wrt Sw
  Teuchos::RCP<CompositeVector> adv_coef = Teuchos::rcp(new CompositeVector(*tmp_flux_));
  adv_coef->PutScalar(0.0);
  Teuchos::RCP<CompositeVector> tmp_coef = Teuchos::rcp(new CompositeVector(*coef_w_->Krel()));
  tmp_coef->PutScalar(0.0);

  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p, *coef_w_->RhoDerivKrel());
  Teuchos::RCP<CompositeVector> rhoDerivLambdaW = Teuchos::rcp(new CompositeVector(*coef_w_->RhoDerivKrel()));
  rhoDerivLambdaW->Scale(dTp/mu1_);
  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_p, *coef_n_->RhoDerivKrel());
  Teuchos::RCP<CompositeVector> rhoDerivLambdaN = Teuchos::rcp(new CompositeVector(*coef_n_->RhoDerivKrel()));
  rhoDerivLambdaN->Scale(dTp/mu2_);
  // coef_n_->RhoDerivKrel()->Scale(dTp / mu2_);

  // Compute RHL * K * lambda_w' * grad Pw
  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->Setup(Kptr, rhoDerivLambdaW, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  //((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->ApplyBCs(true, true);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateFlux(pressure_w.ptr(), tmp_flux_.ptr());
  adv_coef->Update(1.0, *tmp_flux_, 1.0);

  // Compute RHG * K * lambda_n' * grad Pn
  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->Setup(Kptr, rhoDerivLambdaN, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  //((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->ApplyBCs(true, true);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateFlux(pressure_n.ptr(), tmp_flux_.ptr());
  adv_coef->Update(1.0, *tmp_flux_, 1.0);

  // Compute Cg * Pc' * K * lambda_n * grad Pn
  // upwind Pc'
  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_p, *coef_n_->DPc()); 
  tmp_coef->Update(Cg_*dTp/mu2_, *coef_n_->Krel(), 0.0);
  tmp_coef->Multiply(1.0, *tmp_coef, *coef_n_->DPc(), 0.0);

  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->Setup(Kptr, tmp_coef, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  //((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->ApplyBCs(true, true);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateFlux(pressure_n.ptr(), tmp_flux_.ptr());
  adv_coef->Update(0.0, *tmp_flux_, 1.0);

  // Compute RHG * K * lambda_n * grad Pc'
  tmp_flux_->PutScalar(0.0);
  Teuchos::RCP<CompositeVector> rhoLambdaN = Teuchos::rcp(new CompositeVector(*coef_n_->Coeff()));
  rhoLambdaN->Scale(dTp/mu2_);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->Setup(Kptr, rhoLambdaN, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  //((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->ApplyBCs(true, true);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateFlux(coef_n_->DPc().ptr(), tmp_flux_.ptr());
  adv_coef->Update(0.0, *tmp_flux_, 1.0);

  // Compute -phi * D * grad RHL
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  tmp_flux_->PutScalar(0.0);
  op3_matrix_->global_operator()->Init();
  op3_matrix_->Setup(D1ptr, Teuchos::null, Teuchos::null);
  op3_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  //op3_matrix_->ApplyBCs(true, true);
  op3_matrix_->UpdateFlux(rhl.ptr(), tmp_flux_.ptr());
  tmp_flux_->Scale(dTp*phi_);
  adv_coef->Update(1.0, *tmp_flux_, 1.0);

  // diffusion operator wrt Sw
  // RHG * K * lambda_n * Pc'
  tmp_coef->Update(1.0, *coef_n_->Coeff(), 0.0);
  tmp_coef->Multiply(1.0, *tmp_coef, *coef_n_->DPc(), 0.0);
  tmp_coef->Scale(dTp/mu2_);
  op3_preconditioner_->global_operator()->Init();
  op3_preconditioner_->Setup(Kptr, tmp_coef, Teuchos::null);
  op3_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op3_preconditioner_->ApplyBCs(true, true, true);
  //op3_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op3_preconditioner_->global_operator()->AssembleMatrix();
  op3_preconditioner_->global_operator()->Rescale(-1.0);

  //op_prec_sat_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, op3_preconditioner_->global_operator()));
  op_prec_sat_->Setup(*upwind_vn_);
  op_prec_sat_->UpdateMatrices(adv_coef.ptr());
  op_prec_sat_->ApplyBCs(true, true, true);
  op_prec_sat_->global_operator()->Rescale(-1.0);
  //op_prec_sat_->global_operator()->SymbolicAssembleMatrix();
  //op_prec_sat_->global_operator()->AssembleMatrix();

  //op2_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_prec_sat_->global_operator()));
  acc_term.Update(phi_, *rhl, -phi_, *rho_g, 0.0);
  if (dTp > 0.0) {
    op2_acc_->AddAccumulationDelta(*saturation_w, acc_term, acc_term, 1.0, "cell");
  } 
  //op_prec_sat_->global_operator()->SymbolicAssembleMatrix();
  //op_prec_sat_->global_operator()->AssembleMatrix();

  // A_23 block wrt rhl
  //s_with_face->Scale(-1.0);
  op2_preconditioner_->global_operator()->Init();
  op2_preconditioner_->Setup(D1ptr, s_with_face, Teuchos::null);
  op2_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op2_preconditioner_->ApplyBCs(true, true, true);

  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->Setup(Kptr, coef_w_->Krel(), Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateFlux(pressure_w.ptr(), tmp_flux_.ptr());

  //op_prec_rho_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, op2_preconditioner_->global_operator()));
  op_prec_rho_->Setup(*tmp_flux_);
  op_prec_rho_->UpdateMatrices(tmp_flux_.ptr());
  op_prec_rho_->ApplyBCs(true, true, true);

  //op1_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_prec_rho_->global_operator()));
  acc_term.Update(phi_, *saturation_w, 0.0);
  if (dTp > 0.0) {
    op1_acc_->AddAccumulationDelta(*saturation_w, acc_term, acc_term, 1.0, "cell");
  } 
  // op_prec_rho_->global_operator()->SymbolicAssembleMatrix();
  // op_prec_rho_->global_operator()->AssembleMatrix();
}


void CompH_PK::NumericalJacobian(double t_old, double t_new, 
                                 Teuchos::RCP<const TreeVector> u, double eps)
{
  const std::vector<int>& bc_model_p = op_bc_p_->bc_model();

  for (op_iter op_it = ops_.begin(); op_it != ops_.end(); op_it++) {
    (*op_it)->global_operator()->Init();
  }

  Teuchos::RCP<TreeVector> u_copy = Teuchos::rcp(new TreeVector(*u));

  // Create data structures to use functional
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false); 

  Teuchos::RCP<CompositeVector> deriv = Teuchos::rcp(new CompositeVector(cvs));
  deriv->PutScalar(0.0);

  Teuchos::RCP<TreeVector> f_ref = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<CompositeVector> f_ref_cv = Teuchos::rcp(new CompositeVector(cvs));
  f_ref->SetData(f_ref_cv);

  FunctionalResidual(t_old, t_new, u_copy, u_copy, f_ref);

  Teuchos::RCP<TreeVector> u_diff = Teuchos::rcp(new TreeVector(*u));
  Teuchos::RCP<TreeVector> f_diff = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<CompositeVector> f_diff_cv = Teuchos::rcp(new CompositeVector(cvs));
  f_diff->SetData(f_diff_cv);

  for (int i = 0; i < ops_.size(); i++) {
  local_op_iter local_op_it = ops_[i]->global_operator()->begin();
  Epetra_MultiVector& var_diff_c = *u_diff->SubVector(i)->Data()->ViewComponent("cell");

  for (int c = 0; c < var_diff_c.MyLength(); c++) {
    double s_save = var_diff_c[0][c];
    double h;
    h = std::pow(eps, 0.5)*var_diff_c[0][c];
    if (std::abs(h) < 1e-18) h = eps;
    var_diff_c[0][c] += h;
    FunctionalResidual(t_old, t_new, u_copy, u_diff, f_diff);
    var_diff_c[0][c] = s_save;
    deriv->Update(1.0/h, *f_diff->Data(), -1.0/h, *f_ref->Data(), 0.0);
    Epetra_MultiVector& deriv_c = *deriv->ViewComponent("cell");

    AmanziMesh::Entity_ID_List faces;
    mesh_->cell_get_faces(c, &faces);
    int nfaces_none = 0;
    for (int f_it = 0; f_it < faces.size(); ++f_it) {
      int f_id = faces[f_it];
      if (bc_model_p[f_id] != Operators::OPERATOR_BC_NEUMANN) nfaces_none++;
    }
    for (int f_it = 0; f_it < faces.size(); ++f_it) {
      int f_id = faces[f_it];
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f_id, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      // Epetra_MultiVector& deriv_c = *deriv->ViewComponent("cell");
      WhetStone::DenseMatrix Aface(ncells, ncells);
      Aface = 0.0;

      if (bc_model_p[f_id] != Operators::OPERATOR_BC_NEUMANN) {
      for (int i = 0; i != ncells; ++i) {
        if (cells[i] == c) {
          Aface(i, i) = deriv_c[0][cells[i]]/(double)nfaces_none;
          for (int j = i + 1; j != ncells; ++j) {
            Aface(j, i) = deriv_c[0][cells[j]];
          }
        } else {
          for (int j = i + 1; j != ncells; ++j) {
            Aface(i, j) = deriv_c[0][cells[i]];
          }
        }
      }
      WhetStone::DenseMatrix& tmp_matrix = (*local_op_it)->matrices[f_id];
      Aface += tmp_matrix;
      tmp_matrix = Aface;
      }
    }
  }

  ops_[i]->ApplyBCs(true, true, true);
  // ops_[i]->global_operator()->SymbolicAssembleMatrix();
  // ops_[i]->global_operator()->AssembleMatrix();
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
