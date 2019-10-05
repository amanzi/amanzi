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
#include "CompH.hh"
#include "Epetra_Vector.h"
#include "Op.hh"

namespace Amanzi {
namespace Multiphase {

class Op;

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
This is basically the residual
****************************************************************** */
void CompH_PK::Functional(double t_old, double t_new, 
              Teuchos::RCP<TreeVector> u_old,
              Teuchos::RCP<TreeVector> u_new, 
              Teuchos::RCP<TreeVector> f)
{
  //std::cout << "Phase1: Functional\n";
  double dTp(t_new - t_old);
  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_w = u_new->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w = u_new->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl = u_new->SubVector(2)->Data();
  //std::cout << "saturation_n: " << *saturation_n->ViewComponent("cell") << "\n";
  //std::cout << "fugacity: " << *fuga_comp->ViewComponent("cell") << "\n";

  // solution at previous time step
  Teuchos::RCP<const CompositeVector> pressure_w_old = u_old->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w_old = u_old->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl_old = u_old->SubVector(2)->Data(); 

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);  
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);
  //std::cout << "pressure_n: " << *pressure_n->ViewComponent("cell") << "\n";

  // old nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n_old = Teuchos::rcp(new CompositeVector(*pressure_w_old));
  capillary_pressure_->Compute(*saturation_w_old);
  pressure_n_old->Update(1.0, *capillary_pressure_->Pc(), 1.0);

  // compute density of gas rho_g = rho_h_g = P_g * Cg
  Teuchos::RCP<CompositeVector> rho_g = Teuchos::rcp(new CompositeVector(*pressure_n));
  rho_g->Scale(Cg_);
  //std::cout << "rho_g: " << *xw->ViewComponent("cell") << "\n";

  // compute old molar fractions
  Teuchos::RCP<CompositeVector> rho_g_old = Teuchos::rcp(new CompositeVector(*pressure_n_old));
  rho_g_old->Scale(Cg_);
  //std::cout << "rho_g_old: " << *xw->ViewComponent("cell") << "\n";

  // Calculate total mobility needed to initialize diffusion operator
  coef_w_->Compute(*rhl, *saturation_w);
  coef_n_->Compute(*pressure_w, *saturation_w);
  //std::cout << "krel_w before upwind " << *coef_w_->Krel()->ViewComponent("cell") << "\n";
  //std::cout << "krel_n before upwind " << *coef_n_->Krel()->ViewComponent("cell") << "\n";
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);
  upwind_vn_ = S_->GetFieldData("velocity_nonwet", passwd_);
  //std::cout << "upwind velocity: " << *upwind_velocity_->ViewComponent("face") << "\n";
  //std::cout << "x * krel_w before upwind " << *coef_w_->Krel()->ViewComponent("cell") << "\n";
  //std::cout << "x * krel_n before upwind " << *coef_n_->Krel()->ViewComponent("cell") << "\n";

  CoefUpwindFn2 func1 = &MPCoeff::ValueRhoKrel;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_rhl_, bc_value_s_, 
                     *coef_w_->Coeff(), *coef_w_->Coeff(), func1);
  coef_w_->Coeff()->Scale(1.0/mu1_);
  upwind_n_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_p_, bc_value_s_, 
                     *coef_n_->Coeff(), *coef_n_->Coeff(), func1);
  coef_n_->Coeff()->Scale(1.0/mu2_);
  //std::cout << "x * krel_w after upwind " << *coef_w_->Krel()->ViewComponent("face") << "\n";
  //std::cout << "x * krel_n after upwind " << *coef_n_->Krel()->ViewComponent("face") << "\n";
  //std::cout << "coef_w_ after upwind " << *coef_w_->Krel()->ViewComponent("cell") << "\n";
  //std::cout << "coef_n_ after upwind " << *coef_n_->Krel()->ViewComponent("cell") << "\n";

  // compute the upwinded coefficient for the molecular diffusion operator
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);
  Teuchos::RCP<CompositeVector> s_with_face = Teuchos::rcp(new CompositeVector(cvs));
  *s_with_face->ViewComponent("cell") = *saturation_w->ViewComponent("cell");
  //DeriveFaceValuesFromCellValues(*s_with_face->ViewComponent("cell"), *s_with_face->ViewComponent("face"),
  //                               bc_model_, bc_value_s_);
  CoefUpwindFn1 func3 = &MPCoeff::ValuePrimaryVar;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_s_,
                     *s_with_face, *s_with_face, func3);
  s_with_face->Scale(phi_);

  // compute the water component density for gravity term using the density of hydrogen in liquid
  // from the previous time step.
  // rho_w_ = rho_w_std + rho_h_l
  rho_w_->PutScalar(rho_);
  rho_w_->Update(1.0, *rhl_old, 1.0);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->Setup(Kptr, coef_w_->Coeff(), Teuchos::null);
  op1_matrix_->SetDensity(rho_w_);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->ApplyBCs(true, true);
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->SymbolicAssembleMatrix();
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->AssembleMatrix();

  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->Setup(Kptr, coef_n_->Coeff(), Teuchos::null);
  op2_matrix_->SetDensity(rho_g_old);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->ApplyBCs(true, true);
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->SymbolicAssembleMatrix();
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->AssembleMatrix(); 
  //std::cout << "op2_rhs: \n";
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->rhs()->Print(std::cout);
  //std::cout << "((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_): " << *((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->A() << "\n";

  op3_matrix_->global_operator()->Init();
  op3_matrix_->Setup(D1ptr, s_with_face, Teuchos::null);
  op3_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op3_matrix_->ApplyBCs(true, true);
  //op3_matrix_->SymbolicAssembleMatrix();
  //op3_matrix_->AssembleMatrix(); 

  Teuchos::RCP<CompositeVector> f1 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f2 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f3 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->global_operator()->ComputeNegativeResidual(*pressure_w, *f1);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->ComputeNegativeResidual(*pressure_n, *f2);
  op3_matrix_->global_operator()->ComputeNegativeResidual(*rhl, *f3);
  //std::cout << "f1: " << *f1->ViewComponent("cell") << "\n";
  //std::cout << "f2: " << *f2->ViewComponent("cell") << "\n";
  //std::cout << "f3: " << *f3->ViewComponent("cell") << "\n";
  //std::cout << "f4: " << *f4->ViewComponent("cell") << "\n";
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
    //std::cout << "accumulation term cell " << c << ": " << tmp_acc_term << "\n";
    f_cell[0][c] += tmp_acc_term;
  }
  f->Scale(dTp);
  //std::cout << "Functional matrix A: " << *op_matrix_->A_ << "\n";
  //std::cout << "Functional RHS: " << *op_matrix_->rhs_->ViewComponent("cell") << "\n";
  //std::cout << "Functional u_new: " << *u_new->Data()->ViewComponent("cell") << "\n";
  //std::cout << "Functional u_old: " << *u_old->Data()->ViewComponent("cell") << "\n";
  //std::cout << "Functional Residual Component " << comp_name_ << ": " << *f->Data()->ViewComponent("cell") << "\n";
  //std::cout << "Component " << comp_id_ << ": Done Functional.\n";
}

/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int CompH_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                    Teuchos::RCP<TreeVector> Pu) {}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void CompH_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  Teuchos::RCP<const CompositeVector> pressure_w = u->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w = u->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl = u->SubVector(2)->Data();
  //std::cout << "pressure_w: " << *pressure_w->ViewComponent("cell") << "\n";
  //std::cout << "saturation_w: " << *saturation_w->ViewComponent("cell") << "\n";
  //std::cout << "rhl: " << *rhl->ViewComponent("cell") << "\n";

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);  
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);
  //std::cout << "pressure_n: " << *pressure_n->ViewComponent("cell") << "\n";

  // compute density of gas rho_g = rho_h_g = P_g * Cg
  Teuchos::RCP<CompositeVector> rho_g = Teuchos::rcp(new CompositeVector(*pressure_n));
  rho_g->Scale(Cg_);
  //std::cout << "rho_g: " << *rho_g->ViewComponent("cell") << "\n";

  // Calculate total mobility needed to initialize diffusion operator
  coef_w_->Compute(*rhl, *saturation_w);
  coef_n_->Compute(*pressure_w, *saturation_w);
  //std::cout << "krel_w before upwind " << *coef_w_->Krel()->ViewComponent("cell") << "\n";
  //std::cout << "krel_n before upwind " << *coef_n_->Krel()->ViewComponent("cell") << "\n";
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);
  upwind_vn_ = S_->GetFieldData("velocity_nonwet", passwd_);

  // compute flux for upwind direction of advection operator for nonwetting phase
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->Setup(Kptr, Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateFlux(*pressure_n, *upwind_vn_);

  //std::cout << "upwind velocity wet: " << *upwind_vw_->ViewComponent("face") << "\n";
  //std::cout << "upwind velocity nonwet: " << *upwind_vn_->ViewComponent("face") << "\n";

  CoefUpwindFn2 func1 = &MPCoeff::ValueRhoKrel;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_rhl_, bc_value_s_, 
                     *coef_w_->Coeff(), *coef_w_->Coeff(), func1);
  //coef_w_->Coeff()->Scale(dTp/mu1_);
  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_, bc_value_p_, bc_value_s_, 
                     *coef_n_->Coeff(), *coef_n_->Coeff(), func1);
  //coef_n_->Coeff()->Scale(dTp/mu2_);
  Teuchos::RCP<CompositeVector> total_coef = Teuchos::rcp(new CompositeVector(*coef_w_->Coeff()));
  total_coef->Update(dTp/mu2_, *coef_n_->Coeff(), dTp/mu1_);
  //std::cout << "rhl * krel_w cell" << *coef_w_->Coeff()->ViewComponent("cell") << "\n";
  //std::cout << "rhg * krel_n cell" << *coef_n_->Coeff()->ViewComponent("cell") << "\n";
  //std::cout << "coef_w_ after upwind " << *coef_w_->Krel()->ViewComponent("cell") << "\n";
  //std::cout << "coef_n_ after upwind " << *coef_n_->Krel()->ViewComponent("cell") << "\n";

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);
  Teuchos::RCP<CompositeVector> s_with_face = Teuchos::rcp(new CompositeVector(cvs));
  *s_with_face->ViewComponent("cell") = *saturation_w->ViewComponent("cell");
  //DeriveFaceValuesFromCellValues(*s_with_face->ViewComponent("cell"), *s_with_face->ViewComponent("face"),
  //                               bc_model_, bc_value_s_);
  CoefUpwindFn1 func3 = &MPCoeff::ValuePrimaryVar;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_s_,
                     *s_with_face, *s_with_face, func3);
  s_with_face->Scale(dTp*phi_);

  // A_21 block wrt Pw
  op1_preconditioner_->global_operator()->Init();
  op1_preconditioner_->Setup(Kptr, total_coef, Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_preconditioner_->ApplyBCs(true, true);
  //op1_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op1_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "A_21 coef_w_ + coef_n_: " << *op1_preconditioner_->global_operator()->A() << "\n";

  //std::cout << "krel_n before upwind " << *coef_n_->Krel()->ViewComponent("cell") << "\n";
  CoefUpwindFn1 func2 = &MPCoeff::ValueKrel;
  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_, bc_value_s_,
                     *coef_n_->Krel(), *coef_n_->Krel(), func2);
  Teuchos::RCP<CompositeVector> lambda_n = Teuchos::rcp(new CompositeVector(*coef_n_->Krel()));
  lambda_n->Scale(Cg_*dTp/mu2_);

  Teuchos::RCP<CompositeVector> tmp_flux_ = Teuchos::rcp(new CompositeVector(*upwind_vn_));
  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->Setup(Kptr, lambda_n, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateFlux(*pressure_n, *tmp_flux_);

  //Teuchos::ParameterList olist_adv = comp_list_->sublist("operators").sublist("advection operator");
  //op_prec_pres_ = Teuchos::rcp(new Operators::OperatorAdvection(olist_adv, op1_preconditioner_->global_operator()));
  op_prec_pres_->Setup(*tmp_flux_);
  //tmp_flux_->Scale(-1.0);
  op_prec_pres_->UpdateMatrices(*tmp_flux_);
  op_prec_pres_->ApplyBCs(op_bc_p_n_, true);

  // create operator for accumulation term
  //op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_prec_pres_->global_operator()));

  CompositeVector acc_term(*saturation_w);
  acc_term.Scale(-1.0);
  acc_term.Shift(1.0);
  acc_term.Scale(phi_*Cg_);
  if (dTp > 0.0) {
    op_acc_->AddAccumulationTerm(*saturation_w, acc_term, 1.0, "cell");
  } 
  //op_prec_pres_->global_operator()->SymbolicAssembleMatrix();
  //op_prec_pres_->global_operator()->AssembleMatrix();
  //std::cout << "A_21: " << *op_prec_pres_->global_operator()->A() << "\n";

  // A_22 block wrt Sw
  Teuchos::RCP<CompositeVector> adv_coef = Teuchos::rcp(new CompositeVector(*tmp_flux_));
  adv_coef->PutScalar(0.0);
  Teuchos::RCP<CompositeVector> tmp_coef = Teuchos::rcp(new CompositeVector(*coef_w_->Krel()));
  tmp_coef->PutScalar(0.0);

  CoefUpwindFn2 func4 = &MPCoeff::ValueRhoDerivKrel; 
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_rhl_, bc_value_s_, 
                     *coef_w_->RhoDerivKrel(), *coef_w_->RhoDerivKrel(), func4);
  //std::cout << "rhl * deriv krel_w cell" << *coef_w_->RhoDerivKrel()->ViewComponent("cell") << "\n";
  Teuchos::RCP<CompositeVector> rhoDerivLambdaW = Teuchos::rcp(new CompositeVector(*coef_w_->RhoDerivKrel()));
  rhoDerivLambdaW->Scale(dTp/mu1_);
  //std::cout << "rhl * deriv krel_w cell after upwind" << *coef_w_->RhoDerivKrel()->ViewComponent("face") << "\n";
  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_, bc_value_p_, bc_value_s_, 
                     *coef_n_->RhoDerivKrel(), *coef_n_->RhoDerivKrel(), func4);
  //std::cout << "rhg * deriv krel_n cell" << *coef_n_->RhoDerivKrel()->ViewComponent("cell") << "\n";
  Teuchos::RCP<CompositeVector> rhoDerivLambdaN = Teuchos::rcp(new CompositeVector(*coef_n_->RhoDerivKrel()));
  rhoDerivLambdaN->Scale(dTp/mu2_);
  //coef_n_->RhoDerivKrel()->Scale(dTp/mu2_);

  // Compute RHL * K * lambda_w' * grad Pw
  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->Setup(Kptr, rhoDerivLambdaW, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->ApplyBCs(true, true);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->UpdateFlux(*pressure_w, *tmp_flux_);
  //std::cout << "tmp_flux_: " << *tmp_flux_->ViewComponent("face") << "\n";
  adv_coef->Update(1.0, *tmp_flux_, 1.0);
  //std::cout << "adv_coef 1: " << *adv_coef->ViewComponent("face") << "\n";

  // Compute RHG * K * lambda_n' * grad Pn
  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->Setup(Kptr, rhoDerivLambdaN, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->ApplyBCs(true, true);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateFlux(*pressure_n, *tmp_flux_);
  //std::cout << "tmp_flux_: " << *tmp_flux_->ViewComponent("face") << "\n";
  adv_coef->Update(1.0, *tmp_flux_, 1.0);
  //std::cout << "adv_coef 2: " << *adv_coef->ViewComponent("face") << "\n";

  // Compute Cg * Pc' * K * lambda_n * grad Pn
  // upwind Pc'
  func2 = &MPCoeff::DerivativePc;
  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_, bc_value_s_,
                     *coef_n_->DPc(), *coef_n_->DPc(), func2); 
  //std::cout << "Pc': \n"; 
  //coef_n_->DPc()->Print(std::cout);
  tmp_coef->Update(Cg_*dTp/mu2_, *coef_n_->Krel(), 0.0);
  tmp_coef->Multiply(1.0, *tmp_coef, *coef_n_->DPc(), 0.0);
  //tmp_coef->Print(std::cout);

  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->Setup(Kptr, tmp_coef, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->ApplyBCs(true, true);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateFlux(*pressure_n, *tmp_flux_);
  //std::cout << "tmp_flux_: " << *tmp_flux_->ViewComponent("face") << "\n";
  adv_coef->Update(0.0, *tmp_flux_, 1.0);
  //std::cout << "adv_coef 3: " << *adv_coef->ViewComponent("face") << "\n";

  // Compute RHG * K * lambda_n * grad Pc'
  tmp_flux_->PutScalar(0.0);
  Teuchos::RCP<CompositeVector> rhoLambdaN = Teuchos::rcp(new CompositeVector(*coef_n_->Coeff()));
  rhoLambdaN->Scale(dTp/mu2_);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->Setup(Kptr, rhoLambdaN, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->ApplyBCs(true, true);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op2_matrix_)->UpdateFlux(*coef_n_->DPc(), *tmp_flux_);
  //std::cout << "tmp_flux_: " << *tmp_flux_->ViewComponent("face") << "\n";
  adv_coef->Update(0.0, *tmp_flux_, 1.0);
  //std::cout << "adv_coef 4: " << *adv_coef->ViewComponent("face") << "\n";

  // Compute -phi * D * grad RHL
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  tmp_flux_->PutScalar(0.0);
  op3_matrix_->global_operator()->Init();
  op3_matrix_->Setup(D1ptr, Teuchos::null, Teuchos::null);
  op3_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  //op3_matrix_->ApplyBCs(true, true);
  op3_matrix_->UpdateFlux(*rhl, *tmp_flux_);
  tmp_flux_->Scale(dTp*phi_);
  //std::cout << "tmp_flux_: " << *tmp_flux_->ViewComponent("face") << "\n";
  adv_coef->Update(1.0, *tmp_flux_, 1.0);
  //std::cout << "adv_coef 5: " << *adv_coef->ViewComponent("face") << "\n";

  // diffusion operator wrt Sw
  // RHG * K * lambda_n * Pc'
  tmp_coef->Update(1.0, *coef_n_->Coeff(), 0.0);
  tmp_coef->Multiply(1.0, *tmp_coef, *coef_n_->DPc(), 0.0);
  //tmp_coef->Print(std::cout);
  tmp_coef->Scale(dTp/mu2_);
  op3_preconditioner_->global_operator()->Init();
  op3_preconditioner_->Setup(Kptr, tmp_coef, Teuchos::null);
  op3_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op3_preconditioner_->ApplyBCs(true, true);
  //op3_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op3_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "A_22 excluded the advection term: " << *op3_preconditioner_->global_operator()->A() << "\n";
  //EpetraExt::RowMatrixToMatlabFile("A_22_diffusion.txt", *op3_preconditioner_->global_operator()->A());
  op3_preconditioner_->global_operator()->Rescale(-1.0);

  //op_prec_sat_ = Teuchos::rcp(new Operators::OperatorAdvection(olist_adv, op3_preconditioner_->global_operator()));
  op_prec_sat_->Setup(*upwind_vn_);
  op_prec_sat_->UpdateMatrices(*adv_coef);
  op_prec_sat_->ApplyBCs(op_bc_s_, true);
  op_prec_sat_->global_operator()->Rescale(-1.0);
  //op_prec_sat_->global_operator()->SymbolicAssembleMatrix();
  //op_prec_sat_->global_operator()->AssembleMatrix();
  //EpetraExt::RowMatrixToMatlabFile("A_22_before_acc_debug.txt", *op_prec_sat_->global_operator()->A());
  //std::cout << "A_22 excluded the accumulation term: " << *op_prec_sat_->global_operator()->A() << "\n";

  //op2_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_prec_sat_->global_operator()));
  acc_term.Update(phi_, *rhl, -phi_, *rho_g, 0.0);
  //acc_term.Print(std::cout);
  if (dTp > 0.0) {
    op2_acc_->AddAccumulationTerm(*saturation_w, acc_term, 1.0, "cell");
  } 
  //op_prec_sat_->global_operator()->SymbolicAssembleMatrix();
  //op_prec_sat_->global_operator()->AssembleMatrix();
  //EpetraExt::RowMatrixToMatlabFile("A_22_debug.txt", *op_prec_sat_->global_operator()->A());
  //std::cout << "A_22: " << *op_prec_sat_->global_operator()->A() << "\n";

  // A_23 block wrt rhl
  //s_with_face->Scale(-1.0);
  op2_preconditioner_->global_operator()->Init();
  op2_preconditioner_->Setup(D1ptr, s_with_face, Teuchos::null);
  op2_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op2_preconditioner_->ApplyBCs(true, true);

  tmp_flux_->PutScalar(0.0);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->Setup(Kptr, coef_w_->Krel(), Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->UpdateFlux(*pressure_w, *tmp_flux_);

  //op_prec_rho_ = Teuchos::rcp(new Operators::OperatorAdvection(olist_adv, op2_preconditioner_->global_operator()));
  op_prec_rho_->Setup(*tmp_flux_);
  op_prec_rho_->UpdateMatrices(*tmp_flux_);
  op_prec_rho_->ApplyBCs(op_bc_p_, true);

  //op1_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_prec_rho_->global_operator()));
  acc_term.Update(phi_, *saturation_w, 0.0);
  if (dTp > 0.0) {
    op1_acc_->AddAccumulationTerm(*saturation_w, acc_term, 1.0, "cell");
  } 
  //op_prec_rho_->global_operator()->SymbolicAssembleMatrix();
  //op_prec_rho_->global_operator()->AssembleMatrix();
  //std::cout << "A_23: " << *op_prec_rho_->global_operator()->A() << "\n";

}


void CompH_PK::NumericalJacobian(double t_old, double t_new, 
                                Teuchos::RCP<const TreeVector> u, double eps)
{
  for (op_iter op_it = ops_.begin(); op_it != ops_.end(); op_it++) {
    (*op_it)->global_operator()->Init();
  }
  /*
  for (int i = 1; i < ops_.size()-1; i++) {
    ops_[i]->global_operator()->Init();
  }
  */

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

  Functional(t_old, t_new, u_copy, u_copy, f_ref);

  Teuchos::RCP<TreeVector> u_diff = Teuchos::rcp(new TreeVector(*u));
  Teuchos::RCP<TreeVector> f_diff = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<CompositeVector> f_diff_cv = Teuchos::rcp(new CompositeVector(cvs));
  f_diff->SetData(f_diff_cv);

  for (int i = 0; i < ops_.size(); i++)
  {
  local_op_iter local_op_it = ops_[i]->global_operator()->OpBegin();
  Epetra_MultiVector& var_diff_c = *u_diff->SubVector(i)->Data()->ViewComponent("cell");
  for (int c = 0; c < var_diff_c.MyLength(); c++)
  {
    double s_save = var_diff_c[0][c];
    double h;
    h = std::pow(eps, 0.5)*var_diff_c[0][c];
    if (std::abs(h) < 1e-18) h = eps;
    var_diff_c[0][c] += h;
    Functional(t_old, t_new, u_copy, u_diff, f_diff);
    var_diff_c[0][c] = s_save;
    deriv->Update(1.0/h, *f_diff->Data(), -1.0/h, *f_ref->Data(), 0.0);
    Epetra_MultiVector& deriv_c = *deriv->ViewComponent("cell");
    //std::cout << "numerical derivative: " << *deriv->ViewComponent("cell") << "\n";

    AmanziMesh::Entity_ID_List faces;
    mesh_->cell_get_faces(c, &faces);
    int nfaces_none = 0;
    for (int f_it = 0; f_it < faces.size(); ++f_it) {
      int f_id = faces[f_it];
      if (bc_model_[f_id] != Operators::OPERATOR_BC_NEUMANN) nfaces_none++;
    }
    for (int f_it = 0; f_it < faces.size(); ++f_it) {
      int f_id = faces[f_it];
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f_id, AmanziMesh::USED, &cells);
      int ncells = cells.size();
      //std::cout << "Face: " << f_id << "; bc type: " << bc_model_[f_id] << "\n";

      //Epetra_MultiVector& deriv_c = *deriv->ViewComponent("cell");
      //std::cout << "numerical deriv: " << deriv_c << "\n";
      WhetStone::DenseMatrix Aface(ncells, ncells);
      Aface = 0.0;

      if (bc_model_[f_id] != Operators::OPERATOR_BC_NEUMANN)
      {
      for (int i = 0; i != ncells; ++i) {
        //std::cout << "adjacent cells: " << cells[i] << "\n";
        //std::cout << "deriv value diagonal: " << deriv_c[0][cells[i]] << "\n";
        if (cells[i] == c) 
        {
          Aface(i, i) = deriv_c[0][cells[i]]/(double)nfaces_none;
          for (int j = i + 1; j != ncells; ++j) {
            //std::cout << "deriv value off diag: " << deriv_c[0][cells[j]] << "\n";
            Aface(j, i) = deriv_c[0][cells[j]];
          }
        } else {
          for (int j = i + 1; j != ncells; ++j) {
            Aface(i, j) = deriv_c[0][cells[i]];
          }
        }
      }
      WhetStone::DenseMatrix& tmp_matrix = (*local_op_it)->matrices[f_id];
      //std::cout << "current local matrix: \n";
      //std::cout << tmp_matrix(0, 0) << " " << tmp_matrix(0, 1) << "\n";
      //std::cout << tmp_matrix(1, 0) << " " << tmp_matrix(1, 1) << "\n";
      Aface += tmp_matrix;
      //std::cout << "Aface matrix: \n";
      //std::cout << Aface(0, 0) << " " << Aface(0, 1) << "\n";
      //std::cout << Aface(1, 0) << " " << Aface(1, 1) << "\n";
      tmp_matrix = Aface;
      //std::cout << "local matrix after update: \n";
      //std::cout << tmp_matrix(0, 0) << " " << tmp_matrix(0, 1) << "\n";
      //std::cout << tmp_matrix(1, 0) << " " << tmp_matrix(1, 1) << "\n";
      }
    }
  }

  ops_[i]->ApplyBCs(true, true);
  //ops_[i]->global_operator()->SymbolicAssembleMatrix();
  //ops_[i]->global_operator()->AssembleMatrix();

  //std::cout << "Numerical jacobian: " << *ops_[i]->global_operator()->A() << "\n";
  }
} // End NumericalJacobian

} // End namespace Amanzi
} // End namespace Flow