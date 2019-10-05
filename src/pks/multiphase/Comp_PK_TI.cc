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
#include "Comp_PK.hh"
#include "Epetra_Vector.h"
#include "Op.hh"

namespace Amanzi {
namespace Multiphase {

class Op;

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
This is basically the residual
****************************************************************** */
void Comp_PK::Functional(double t_old, double t_new, 
              Teuchos::RCP<TreeVector> u_old,
              Teuchos::RCP<TreeVector> u_new, 
              Teuchos::RCP<TreeVector> f)
{
  //std::cout << "Phase1: Functional\n";
  double dTp(t_new - t_old);
  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_w = u_new->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_n = u_new->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> fuga_comp = u_new->SubVector(1+comp_id_)->Data();
  //std::cout << "saturation_n: " << *saturation_n->ViewComponent("cell") << "\n";
  //std::cout << "fugacity: " << *fuga_comp->ViewComponent("cell") << "\n";

  // solution at previous time step
  Teuchos::RCP<const CompositeVector> pressure_w_old = u_old->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_n_old = u_old->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> fuga_comp_old = u_old->SubVector(1+comp_id_)->Data(); 

  // compute saturation non-wetting phase and pressure wetting phase  
  Teuchos::RCP<CompositeVector> saturation_w = Teuchos::rcp(new CompositeVector(*saturation_n));
  saturation_w->Scale(-1.0);
  saturation_w->Shift(1.0);

  // old wet saturation
  Teuchos::RCP<CompositeVector> saturation_w_old = Teuchos::rcp(new CompositeVector(*saturation_n_old));
  saturation_w_old->Scale(-1.0);
  saturation_w_old->Shift(1.0);

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);  
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);
  //std::cout << "pressure_n: " << *pressure_n->ViewComponent("cell") << "\n";

  // old nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n_old = Teuchos::rcp(new CompositeVector(*pressure_w_old));
  capillary_pressure_old_->Compute(*saturation_w_old);
  pressure_n_old->Update(1.0, *capillary_pressure_old_->Pc(), 1.0);

  // compute molar fractions
  Teuchos::RCP<CompositeVector> xw = Teuchos::rcp(new CompositeVector(*fuga_comp));
  Teuchos::RCP<CompositeVector> xn = Teuchos::rcp(new CompositeVector(*fuga_comp));
  xn->ReciprocalMultiply(1.0, *pressure_n, *fuga_comp, 0.0);
  xw->Scale(1.0/comp_coeff_);
  //std::cout << "xw: " << *xw->ViewComponent("cell") << "\n";
  //std::cout << "xn: " << *xn->ViewComponent("cell") << "\n";

  // compute old molar fractions
  Teuchos::RCP<CompositeVector> xw_old = Teuchos::rcp(new CompositeVector(*fuga_comp_old));
  Teuchos::RCP<CompositeVector> xn_old = Teuchos::rcp(new CompositeVector(*fuga_comp_old));
  xn_old->ReciprocalMultiply(1.0, *pressure_n_old, *fuga_comp_old, 0.0);
  xw_old->Scale(1.0/comp_coeff_);

  /*
  std::cout << "phase_pressure: " << *pressure_w->ViewComponent("cell") << "\n";
  std::cout << "phase_saturation: " << *saturation_n->ViewComponent("cell") << "\n";
  std::cout << "x11: " <<  fuga_comp ->ViewComponent("cell") << "\n";
  std::cout << "pressure old: " << *u_old->SubVector(0)->Data()->ViewComponent("cell") << "\n";
  std::cout << "saturation old: " << *u_old->SubVector(1)->Data()->ViewComponent("cell") << "\n";
  std::cout << "x11 old: " << *u_old->SubVector(2)->Data()->ViewComponent("cell") << "\n";
  */

  // Create data structures to initialize operators
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);
  Teuchos::RCP<CompositeVector> Krel = Teuchos::rcp(new CompositeVector(cvs));
  Teuchos::RCP<CompositeVector> xw_wface = Teuchos::rcp(new CompositeVector(cvs));
  Teuchos::RCP<CompositeVector> xn_wface = Teuchos::rcp(new CompositeVector(cvs));
  Krel->PutScalarMasterAndGhosted(1.0);
  xw_wface->PutScalarMasterAndGhosted(0.0);
  xn_wface->PutScalarMasterAndGhosted(0.0);
  Epetra_MultiVector& xw_wf_c = *xw_wface->ViewComponent("cell");
  Epetra_MultiVector& xn_wf_c = *xn_wface->ViewComponent("cell");
  const Epetra_MultiVector& xw_c = *xw->ViewComponent("cell");
  const Epetra_MultiVector& xn_c = *xn->ViewComponent("cell");
  for (int c = 0; c < xw_c.MyLength(); c++) {
    xw_wf_c[0][c] = xw_c[0][c];
    xn_wf_c[0][c] = xn_c[0][c];
  }

  // Calculate total mobility needed to initialize diffusion operator
  rel_perm_w_->Compute(*saturation_w);
  rel_perm_n_->Compute(*saturation_w);
  //std::cout << "krel_w before upwind " << *rel_perm_w_->Krel()->ViewComponent("cell") << "\n";
  //std::cout << "krel_n before upwind " << *rel_perm_n_->Krel()->ViewComponent("cell") << "\n";
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);
  upwind_vn_ = S_->GetFieldData("velocity_nonwet", passwd_);
  //std::cout << "upwind velocity: " << *upwind_velocity_->ViewComponent("face") << "\n";
  rel_perm_w_->Krel()->Multiply(1.0, *rel_perm_w_->Krel(), *xw_wface, 0.0);
  rel_perm_n_->Krel()->Multiply(1.0, *rel_perm_n_->Krel(), *xn_wface, 0.0);
  //std::cout << "x * krel_w before upwind " << *rel_perm_w_->Krel()->ViewComponent("cell") << "\n";
  //std::cout << "x * krel_n before upwind " << *rel_perm_n_->Krel()->ViewComponent("cell") << "\n";

  FracKrelUpwindFn func1 = &RelativePermeability::Value;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_s_, bc_value_xw_, bc_value_s_, 
                     *rel_perm_w_->Krel(), *rel_perm_w_->Krel(), func1);
  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_s_, bc_value_xn_, bc_value_s_, 
                     *rel_perm_n_->Krel(), *rel_perm_n_->Krel(), func1);
  //std::cout << "x * krel_w after upwind " << *rel_perm_w_->Krel()->ViewComponent("face") << "\n";
  //std::cout << "x * krel_n after upwind " << *rel_perm_n_->Krel()->ViewComponent("face") << "\n";
  //std::cout << "krel after upwind " << *rel_perm_w_->Krel()->ViewComponent("face") << "\n";

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D2ptr = Teuchos::rcpFromRef(D2_);
  op1_matrix_->global_operator()->Init();
  op1_matrix_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_matrix_->ApplyBCs(true, true);
  //op1_matrix_->SymbolicAssembleMatrix();
  //op1_matrix_->AssembleMatrix();

  op2_matrix_->global_operator()->Init();
  op2_matrix_->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
  op2_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op2_matrix_->ApplyBCs(true, true);
  //op2_matrix_->SymbolicAssembleMatrix();
  //op2_matrix_->AssembleMatrix(); 

  op3_matrix_->global_operator()->Init();
  op3_matrix_->Setup(D1ptr, Teuchos::null, Teuchos::null);
  op3_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op3_matrix_->ApplyBCs(true, true);
  //op3_matrix_->SymbolicAssembleMatrix();
  //op3_matrix_->AssembleMatrix(); 

  op4_matrix_->global_operator()->Init();
  op4_matrix_->Setup(D2ptr, Teuchos::null, Teuchos::null);
  op4_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op4_matrix_->ApplyBCs(true, true);
  //op4_matrix_->SymbolicAssembleMatrix();
  //op4_matrix_->AssembleMatrix(); 

  // Update source term
  //Teuchos::RCP<CompositeVector> rhs = op1_matrix_->rhs();
  //if (src_sink_ != NULL) AddSourceTerms(*rhs);
  
  Teuchos::RCP<CompositeVector> f1 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f2 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f3 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f4 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  op1_matrix_->global_operator()->ComputeNegativeResidual(*u_new->SubVector(0)->Data(), *f1);
  op2_matrix_->global_operator()->ComputeNegativeResidual(*pressure_n, *f2);
  op3_matrix_->global_operator()->ComputeNegativeResidual(*xw, *f3);
  op4_matrix_->global_operator()->ComputeNegativeResidual(*xn, *f4);
  //std::cout << "f1: " << *f1->ViewComponent("cell") << "\n";
  //std::cout << "f2: " << *f2->ViewComponent("cell") << "\n";
  //std::cout << "f3: " << *f3->ViewComponent("cell") << "\n";
  //std::cout << "f4: " << *f4->ViewComponent("cell") << "\n";
  f->Data()->Update(1.0, *f1, 1.0, *f2, 0.0);
  f->Data()->Update(1.0, *f3, 1.0, *f4, 1.0);
  
  // Add time derivative
  Epetra_MultiVector& f_cell = *f->Data()->ViewComponent("cell", true);

  const Epetra_MultiVector& S_new_c = *u_new->SubVector(1)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& S_old_c = *u_old->SubVector(1)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& X1_new_c = *xw->ViewComponent("cell");
  const Epetra_MultiVector& X1_old_c = *xw_old->ViewComponent("cell");
  const Epetra_MultiVector& X2_new_c = *xn->ViewComponent("cell");
  const Epetra_MultiVector& X2_old_c = *xn_old->ViewComponent("cell");

  // Add accumulation term
  double s0, s1, mfw_0, mfw, mfn_0, mfn, volume;
  for (int c = 0; c < f_cell.MyLength(); c++) {
    s1 = S_new_c[0][c];
    s0 = S_old_c[0][c];
    mfw = X1_new_c[0][c];
    mfw_0 = X1_old_c[0][c];
    mfn = X2_new_c[0][c];
    mfn_0 = X2_old_c[0][c];

    double factor = phi_ * mesh_->cell_volume(c) / dTp;
    double tmp_acc_term = (mfw*(1.0-s1) + mfn*s1 - mfw_0*(1.0-s0) - mfn_0*s0) * factor;
    //std::cout << "accumulation term cell " << c << ": " << tmp_acc_term << "\n";
    f_cell[0][c] += (mfw*(1.0-s1) + mfn*s1 - mfw_0*(1.0-s0) - mfn_0*s0) * factor;
  }
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
int Comp_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                    Teuchos::RCP<TreeVector> Pu) {}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void Comp_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  Teuchos::RCP<const CompositeVector> pressure_w = u->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_n = u->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> fuga_comp = u->SubVector(1+comp_id_)->Data();

  // compute saturation non-wetting phase and pressure wetting phase  
  Teuchos::RCP<CompositeVector> saturation_w = Teuchos::rcp(new CompositeVector(*saturation_n));
  saturation_w->Scale(-1.0);
  saturation_w->Shift(1.0);

  // new nonwet pressure
  Teuchos::RCP<CompositeVector> pressure_n = Teuchos::rcp(new CompositeVector(*pressure_w));
  capillary_pressure_->Compute(*saturation_w);  
  pressure_n->Update(1.0, *capillary_pressure_->Pc(), 1.0);
  //std::cout << "pressure_n: " << *pressure_n->ViewComponent("cell") << "\n";

  // compute molar fractions
  Teuchos::RCP<CompositeVector> xw = Teuchos::rcp(new CompositeVector(*fuga_comp));
  Teuchos::RCP<CompositeVector> xn = Teuchos::rcp(new CompositeVector(*fuga_comp));
  xn->ReciprocalMultiply(1.0, *pressure_n, *fuga_comp, 0.0);
  xw->Scale(1.0/comp_coeff_);
  //std::cout << "xw: " << *xw->ViewComponent("cell") << "\n";
  //std::cout << "xn: " << *xn->ViewComponent("cell") << "\n";

  // Create data structures to initialize operators
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);
  Teuchos::RCP<CompositeVector> xw_wface = Teuchos::rcp(new CompositeVector(cvs));
  Teuchos::RCP<CompositeVector> xn_wface = Teuchos::rcp(new CompositeVector(cvs));
  xw_wface->PutScalarMasterAndGhosted(0.0);
  xn_wface->PutScalarMasterAndGhosted(0.0);
  Epetra_MultiVector& xw_wf_c = *xw_wface->ViewComponent("cell");
  Epetra_MultiVector& xn_wf_c = *xn_wface->ViewComponent("cell");
  const Epetra_MultiVector& xw_c = *xw->ViewComponent("cell");
  const Epetra_MultiVector& xn_c = *xn->ViewComponent("cell");
  for (int c = 0; c < xw_c.MyLength(); c++) {
    xw_wf_c[0][c] = xw_c[0][c];
    xn_wf_c[0][c] = xn_c[0][c];
  }

  // Calculate total mobility needed to initialize diffusion operator
  rel_perm_w_->Compute(*saturation_w);
  rel_perm_n_->Compute(*saturation_w);
  Teuchos::RCP<CompositeVector> krel_w_copy = Teuchos::rcp(new CompositeVector(*rel_perm_w_->Krel()));
  Teuchos::RCP<CompositeVector> krel_n_copy = Teuchos::rcp(new CompositeVector(*rel_perm_n_->Krel()));
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);
  upwind_vn_ = S_->GetFieldData("velocity_nonwet", passwd_);
  //std::cout << "upwind velocity: " << *upwind_velocity_->ViewComponent("face") << "\n";
  rel_perm_w_->Krel()->Multiply(1.0, *rel_perm_w_->Krel(), *xw_wface, 0.0);
  rel_perm_n_->Krel()->Multiply(1.0, *rel_perm_n_->Krel(), *xn_wface, 0.0);
  //std::cout << "krel_w before upwind " << *rel_perm_w_->Krel()->ViewComponent("cell") << "\n";
  //std::cout << "krel_n before upwind " << *rel_perm_n_->Krel()->ViewComponent("cell") << "\n";

  FracKrelUpwindFn func1 = &RelativePermeability::Value;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_s_, bc_value_xw_, bc_value_s_, 
                     *rel_perm_w_->Krel(), *rel_perm_w_->Krel(), func1);

  upwind_n_->Compute(*upwind_vn_, *upwind_vn_, bc_model_s_, bc_value_xn_, bc_value_s_, 
                     *rel_perm_n_->Krel(), *rel_perm_n_->Krel(), func1);

  FracKrelUpwindFn func2 = &RelativePermeability::Derivative;
  rel_perm_w_->dKdS()->Scale(-1.0);
  rel_perm_n_->dKdS()->Scale(-1.0);
  upwind_n_->Compute(*upwind_vw_, *upwind_vw_, bc_model_s_, bc_value_xw_, bc_value_s_, 
                     *rel_perm_w_->dKdS(), *rel_perm_w_->dKdS(), func2);

  upwind_w_->Compute(*upwind_vn_, *upwind_vn_, bc_model_s_, bc_value_xn_, bc_value_s_, 
                     *rel_perm_n_->dKdS(), *rel_perm_n_->dKdS(), func2);

  Teuchos::RCP<CompositeVector> dKdS = Teuchos::rcp(new CompositeVector(cvs));
  Teuchos::RCP<CompositeVector> Krel_m = Teuchos::rcp(new CompositeVector(cvs));
  Teuchos::RCP<CompositeVector> Krel = Teuchos::rcp(new CompositeVector(cvs));
  Teuchos::RCP<CompositeVector> porosity = Teuchos::rcp(new CompositeVector(cvs));
  Krel->PutScalarMasterAndGhosted(0.0);
  Krel_m->PutScalarMasterAndGhosted(1.0);
  dKdS->PutScalarMasterAndGhosted(0.0);
  //porosity->PutScalarMasterAndGhosted(phi_);
  porosity->PutScalarMasterAndGhosted(1.0);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D2ptr = Teuchos::rcpFromRef(D2_);
  op1_preconditioner_->global_operator()->Init();
  op1_preconditioner_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_preconditioner_->ApplyBCs(true, true);
  //std::cout << "Done op1 \n";

  op1_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  op1_preconditioner_->global_operator()->AssembleMatrix();
  std::cout << "Matrix prec 1: " << *op1_preconditioner_->global_operator()->A() << "\n";

  Teuchos::ParameterList& op_list = comp_list_->sublist("operators").sublist("diffusion operator").sublist("preconditioner");
  op_pres_prec_ = Teuchos::rcp(new Operators::OperatorDiffusionFV(op_list, op1_preconditioner_->global_operator()));
  op_pres_prec_->SetBCs(op_bc_p_, op_bc_p_);
  op_pres_prec_->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
  op_pres_prec_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_pres_prec_->ApplyBCs(true, true);
  //std::cout << "Done op_pres \n";

  Teuchos::RCP<CompositeVector> accumulation_pres = Teuchos::rcp(new CompositeVector(*saturation_n));
  accumulation_pres->Multiply(1.0, *xn, *accumulation_pres, 0.0);
  accumulation_pres->ReciprocalMultiply(-1.0, *pressure_n, *accumulation_pres, 0.0);
  accumulation_pres->Scale(phi_);
  op_pres_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_pres_prec_->global_operator()));
  if (dTp > 0.0) {
    op_pres_acc_->AddAccumulationTerm(*saturation_n, *accumulation_pres, dTp, "cell");
  }

  op_pres_prec_->global_operator()->SymbolicAssembleMatrix();
  op_pres_prec_->global_operator()->AssembleMatrix();
  std::cout << "Matrix pres prec: " << *op_pres_prec_->global_operator()->A() << "\n";

  //Teuchos::RCP<CompositeVector> tmp_flux1_ = Teuchos::rcp(new CompositeVector(*upwind_vw_));
  //Teuchos::RCP<CompositeVector> tmp_flux2_ = Teuchos::rcp(new CompositeVector(*upwind_vn_));
  tmp_flux1_->PutScalar(0.0);
  tmp_flux2_->PutScalar(0.0);
  op1_matrix_->global_operator()->Init();
  op1_matrix_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op1_matrix_->UpdateFlux(*pressure_w, *tmp_flux1_);
  op1_matrix_->UpdateFlux(*pressure_n, *tmp_flux2_);
  //std::cout << "tmp_flux1_: " << *tmp_flux1_->ViewComponent("face");
  //std::cout << "tmp_flux2_: " << *tmp_flux2_->ViewComponent("face");

  Epetra_MultiVector& fface1 = *tmp_flux1_->ViewComponent("face");
  Epetra_MultiVector& dKdS_f_w = *rel_perm_w_->dKdS()->ViewComponent("face");
  ASSERT(dKdS_f_w.MyLength() == fface1.MyLength());
  for (int f = 0; f < dKdS_f_w.MyLength(); f++) {
    dKdS_f_w[0][f] *= fface1[0][f]; 
  }
  Epetra_MultiVector& fface2 = *tmp_flux2_->ViewComponent("face");
  Epetra_MultiVector& dKdS_f_n = *rel_perm_n_->dKdS()->ViewComponent("face");
  ASSERT(dKdS_f_n.MyLength() == fface2.MyLength());
  for (int f = 0; f < dKdS_f_n.MyLength(); f++) {
    dKdS_f_n[0][f] *= fface2[0][f]; 
  }
  //std::cout << "dKdS_w * tmp_flux1_: " << *rel_perm_w_->dKdS()->ViewComponent("face");
  //std::cout << "dKdS_n * tmp_flux2_: " << *rel_perm_n_->dKdS()->ViewComponent("face");

  //std::cout << "Operator 2\n";
  op2_preconditioner_->global_operator()->Init();
  op2_preconditioner_->Setup(*tmp_flux1_);
  op2_preconditioner_->UpdateMatrices(*rel_perm_w_->dKdS()); 
  op2_preconditioner_->ApplyBCs(op_bc_p_, true);
  op2_preconditioner_->global_operator()->Rescale(-1.0);
  //std::cout << "Done op2 \n";

  Teuchos::ParameterList adv_list = comp_list_->sublist("operators").sublist("advection operator");
  op_sat_prec_ = Teuchos::rcp(new Operators::OperatorAdvection(adv_list, op2_preconditioner_->global_operator()));
  op_sat_prec_->Setup(*tmp_flux2_);
  op_sat_prec_->UpdateMatrices(*rel_perm_n_->dKdS());
  op_sat_prec_->ApplyBCs(op_bc_p_, true);
  //std::cout << "Done op_sat \n";
  op_sat_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_sat_prec_->global_operator()));

  CompositeVector accumulation_factor(*xn);
  accumulation_factor.Update(-1.0, *xw, 1.0);
  accumulation_factor.Scale(phi_);
  if (dTp > 0.0) {
    op_sat_acc_->AddAccumulationTerm(*saturation_n, accumulation_factor, dTp, "cell");
  }

  Teuchos::RCP<CompositeVector> mf_coef1 = Teuchos::rcp(new CompositeVector(cvs));
  mf_coef1->PutScalarMasterAndGhosted(1.0);
  mf_coef1->Scale(1.0/comp_coeff_);
  //std::cout << "mf_coef1: " << *mf_coef1->ViewComponent("face") << "\n";
  op3_preconditioner_->global_operator()->Init();
  op3_preconditioner_->Setup(D1ptr, mf_coef1, Teuchos::null);
  op3_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op3_preconditioner_->ApplyBCs(true, true);

  //op3_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op3_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "Matrix prec 3: " << *op3_preconditioner_->global_operator()->A() << "\n";

  Teuchos::RCP<CompositeVector> mf_coef2 = Teuchos::rcp(new CompositeVector(cvs));
  Teuchos::RCP<CompositeVector> pressure_n_wf = Teuchos::rcp(new CompositeVector(cvs));
  DeriveFaceValuesFromCellValues(*pressure_n->ViewComponent("cell"), *pressure_n_wf->ViewComponent("face"),
                                 bc_model_pc_, bc_value_pc_);
  //std::cout << "pressure_n_wf: " << *pressure_n_wf->ViewComponent("face") << "\n";
  mf_coef2->PutScalarMasterAndGhosted(1.0);
  mf_coef2->ReciprocalMultiply(1.0, *pressure_n_wf, *mf_coef2, 0.0);
  //std::cout << "mf_coef2: " << *mf_coef2->ViewComponent("face") << "\n";
  op4_preconditioner_ = Teuchos::rcp(new Operators::OperatorDiffusionFV(op_list, op3_preconditioner_->global_operator()));
  op4_preconditioner_->SetBCs(op_bc_xn_, op_bc_xn_);
  op4_preconditioner_->Setup(D2ptr, mf_coef2, Teuchos::null);
  op4_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op4_preconditioner_->ApplyBCs(true, true);

  //op4_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op4_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "Matrix prec 4: " << *op4_preconditioner_->global_operator()->A() << "\n";

  tmp_flux1_->PutScalar(0.0);
  tmp_flux2_->PutScalar(0.0);
  op1_matrix_->global_operator()->Init();
  op1_matrix_->Setup(Kptr, krel_w_copy, Teuchos::null);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_matrix_->UpdateFlux(*pressure_w, *tmp_flux1_);
  tmp_flux1_->Scale(1.0/comp_coeff_);

  op2_matrix_->global_operator()->Init();
  op2_matrix_->Setup(Kptr, krel_n_copy, Teuchos::null);
  op2_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op2_matrix_->UpdateFlux(*pressure_n, *tmp_flux2_);
  tmp_flux2_->ReciprocalMultiply(1.0, *pressure_n_wf, *tmp_flux2_, 0.0);

  op5_preconditioner_ = Teuchos::rcp(new Operators::OperatorAdvection(adv_list, op4_preconditioner_->global_operator()));
  op5_preconditioner_->Setup(*tmp_flux1_);
  op5_preconditioner_->UpdateMatrices(*tmp_flux1_);
  op5_preconditioner_->ApplyBCs(op_bc_p_, true);

  //op5_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op5_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "Matrix prec 5: " << *op5_preconditioner_->global_operator()->A() << "\n";

  op6_preconditioner_ = Teuchos::rcp(new Operators::OperatorAdvection(adv_list, op5_preconditioner_->global_operator()));
  op6_preconditioner_->Setup(*tmp_flux2_);
  op6_preconditioner_->UpdateMatrices(*tmp_flux2_);
  op6_preconditioner_->ApplyBCs(op_bc_p_, true);

  //op6_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op6_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "Matrix prec 6: " << *op6_preconditioner_->global_operator()->A() << "\n";

  op_fug_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op6_preconditioner_->global_operator()));
  CompositeVector accumulation_fuga(*saturation_n);
  accumulation_fuga.ReciprocalMultiply(1.0, *pressure_n, accumulation_fuga, 0.0);
  accumulation_fuga.Update(1.0/comp_coeff_, *saturation_w, 1.0);
  accumulation_fuga.Scale(phi_);
  if (dTp > 0.0) {
    op_fug_acc_->AddAccumulationTerm(*saturation_n, accumulation_fuga, dTp, "cell");
  }

  //op6_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op6_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "Matrix prec 4: " << *op6_preconditioner_->global_operator()->A() << "\n";
  /*
  op_pres_prec_->global_operator()->SymbolicAssembleMatrix();
  op_pres_prec_->global_operator()->AssembleMatrix();
  op_sat_prec_->global_operator()->SymbolicAssembleMatrix();
  op_sat_prec_->global_operator()->AssembleMatrix();
  */
}


void Comp_PK::NumericalJacobian(double t_old, double t_new, 
                                Teuchos::RCP<const TreeVector> u, double eps)
{
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
      if (bc_model_p_[f_id] != Operators::OPERATOR_BC_NEUMANN) nfaces_none++;
    }
    for (int f_it = 0; f_it < faces.size(); ++f_it) {
      int f_id = faces[f_it];
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f_id, AmanziMesh::USED, &cells);
      int ncells = cells.size();
      //std::cout << "Face: " << f_id << "; bc type: " << bc_model_p_[f_id] << "\n";

      //Epetra_MultiVector& deriv_c = *deriv->ViewComponent("cell");
      //std::cout << "numerical deriv: " << deriv_c << "\n";
      WhetStone::DenseMatrix Aface(ncells, ncells);
      Aface = 0.0;

      if (bc_model_p_[f_id] != Operators::OPERATOR_BC_NEUMANN)
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