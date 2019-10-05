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
#include "CompW.hh"
#include "Epetra_Vector.h"
#include "Op.hh"

namespace Amanzi {
namespace Multiphase {

class Op;

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
This is basically the residual
****************************************************************** */
void CompW_PK::Functional(double t_old, double t_new, 
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

  // solution at previous time step
  Teuchos::RCP<const CompositeVector> pressure_w_old = u_old->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w_old = u_old->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl_old = u_old->SubVector(2)->Data(); 

  //std::cout << "phase_pressure: " << *pressure_w->ViewComponent("cell") << "\n";
  //std::cout << "phase_saturation: " << *saturation_w->ViewComponent("cell") << "\n";
  /*
  std::cout << "rhl: " <<  rhl ->ViewComponent("cell") << "\n";
  std::cout << "pressure old: " << *u_old->SubVector(0)->Data()->ViewComponent("cell") << "\n";
  std::cout << "saturation old: " << *u_old->SubVector(1)->Data()->ViewComponent("cell") << "\n";
  std::cout << "rhl old: " << *u_old->SubVector(2)->Data()->ViewComponent("cell") << "\n";
  */

  // Calculate total mobility needed to initialize diffusion operator
  rel_perm_w_->Compute(*saturation_w);
  //std::cout << "krel_w before upwind " << *rel_perm_w_->Krel()->ViewComponent("cell") << "\n";
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);
  //std::cout << "upwind velocity: " << *upwind_vw_->ViewComponent("face") << "\n";

  CoefUpwindFn1 func1 = &MPCoeff::ValueKrel;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_s_, 
                     *rel_perm_w_->Krel(), *rel_perm_w_->Krel(), func1);
  rel_perm_w_->Krel()->Scale(rho_/mu_);
  //std::cout << "krel after upwind " << *rel_perm_w_->Krel()->ViewComponent("face") << "\n";

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
  CoefUpwindFn1 func2 = &MPCoeff::ValuePrimaryVar;
  upwind_vw_->Scale(-1.0);
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_s_,
                     *s_with_face, *s_with_face, func2);
  s_with_face->Scale(-phi_);
  upwind_vw_->Scale(-1.0);
  //std::cout << "s_with_face: " << *s_with_face->ViewComponent("cell") << "\n";

  // compute the water component density for gravity term using the density of hydrogen in liquid
  // from the previous time step.
  // rho_w_ = rho_w_std + rho_h_l
  rho_w_->PutScalar(rho_);
  rho_w_->Update(1.0, *rhl_old, 1.0);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_matrix_->SetDensity(rho_w_);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->ApplyBCs(true, true);
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->SymbolicAssembleMatrix();
  //((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->AssembleMatrix();

  op2_matrix_->global_operator()->Init();
  op2_matrix_->Setup(D1ptr, s_with_face, Teuchos::null);
  op2_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op2_matrix_->ApplyBCs(true, true);
  //op2_matrix_->SymbolicAssembleMatrix();
  //op2_matrix_->AssembleMatrix(); 

  // Update source term
  //Teuchos::RCP<CompositeVector> rhs = ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->rhs();
  //if (src_sink_ != NULL) AddSourceTerms(*rhs);
  
  Teuchos::RCP<CompositeVector> f1 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f2 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->global_operator()->ComputeNegativeResidual(*pressure_w, *f1);
  op2_matrix_->global_operator()->ComputeNegativeResidual(*rhl, *f2);
  //std::cout << "f1: " << *f1->ViewComponent("cell") << "\n";
  //std::cout << "f2: " << *f2->ViewComponent("cell") << "\n";
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
    //std::cout << "accumulation term cell " << c << ": " << tmp_acc_term << "\n";
    f_cell[0][c] += rho_ * (s1 - s0) * factor;
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
int CompW_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                    Teuchos::RCP<TreeVector> Pu) {}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void CompW_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_w = u->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_w = u->SubVector(1)->Data();
  Teuchos::RCP<const CompositeVector> rhl = u->SubVector(2)->Data();

  // Calculate relative perm needed to initialize diffusion operator
  rel_perm_w_->Compute(*saturation_w);
  //std::cout << "krel_w before upwind " << *rel_perm_w_->Krel()->ViewComponent("cell") << "\n";
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);
  //std::cout << "upwind velocity: " << *upwind_vw_->ViewComponent("face") << "\n";

  CoefUpwindFn1 func1 = &MPCoeff::ValueKrel;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_s_, 
                     *rel_perm_w_->Krel(), *rel_perm_w_->Krel(), func1);
  rel_perm_w_->Krel()->Scale(dTp*rho_/mu_);
  //std::cout << "Krel face: " << *rel_perm_w_->Krel()->ViewComponent("face") << "\n";

  CoefUpwindFn1 func2 = &MPCoeff::DerivativeKrel;
  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, bc_value_s_, *rel_perm_w_->dKdS(), *rel_perm_w_->dKdS(), func2);
  //std::cout << "dKdS face: " << *rel_perm_w_->dKdS()->ViewComponent("face") << "\n";
  rel_perm_w_->dKdS()->Scale(dTp*rho_/mu_); 
  //rel_perm_w_->dKdS()->Scale(0.0);

  // the flux is only used for finding the upwind cells
  Teuchos::RCP<CompositeVector> tmp_flux_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_wet")));
  tmp_flux_->PutScalar(0.0);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->Setup(Kptr, Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::OperatorDiffusion>)op1_matrix_)->UpdateFlux(*pressure_w, *tmp_flux_);

  // rel_perm_w_->dKdS() is used as the coefficients of the matrix
  //Epetra_MultiVector& krel_f = *rel_perm_w_->Krel()->ViewComponent("face");
  Epetra_MultiVector& fface = *tmp_flux_->ViewComponent("face", true);
  Epetra_MultiVector& dKdS_f = *rel_perm_w_->dKdS()->ViewComponent("face", true);
  ASSERT(dKdS_f.MyLength() == fface.MyLength());
  for (int f = 0; f < dKdS_f.MyLength(); f++) {
    dKdS_f[0][f] *= fface[0][f]; 
    //std::cout << "dKdS: " << dKdS_f[0][f] << "; fface: " << fface[0][f] << "\n";
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
  op1_preconditioner_->ApplyBCs(true, true);
  //op1_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op1_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "Analytic Jacobian A_11: " << *op1_preconditioner_->global_operator()->A() << "\n";

  // A_12 block wrt Sw
  //tmp_flux_->Scale(-1.0);
  op2_preconditioner_->global_operator()->Init();
  op2_preconditioner_->Setup(*tmp_flux_);
  op2_preconditioner_->UpdateMatrices(*rel_perm_w_->dKdS());
  op2_preconditioner_->ApplyBCs(op_bc_s_, true);
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
  op2_matrix_->UpdateFlux(*rhl, *tmp_flux_);
  //std::cout << "tmp_flux_ for rhl: " << *tmp_flux_->ViewComponent("face") << "\n";

  //Teuchos::ParameterList olist_adv = comp_list_->sublist("operators").sublist("advection operator");
  //op_prec_sat_ = Teuchos::rcp(new Operators::OperatorAdvection(olist_adv, op2_preconditioner_->global_operator()));
  //tmp_flux_->Scale(-1.0);
  op_prec_sat_->Setup(*tmp_flux_);
  tmp_flux_->Scale(-1.0);
  op_prec_sat_->UpdateMatrices(*tmp_flux_);
  op_prec_sat_->ApplyBCs(op_bc_s_, true);
  op_prec_sat_->global_operator()->Rescale(-1.0);

  // create operator for accumulation term
  //op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_prec_sat_->global_operator()));

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
    op_acc_->AddAccumulationTerm(*saturation_w, porosity, 1.0, "cell");
  } 
  //op_prec_sat_->global_operator()->SymbolicAssembleMatrix();
  //op_prec_sat_->global_operator()->AssembleMatrix();
  //std::cout << "Analytic Jacobian A_12: " << *op_prec_sat_->global_operator()->A() << "\n";

  // A_13 block wrt rhl
  //s_with_face->Scale(-1.0);
  op3_preconditioner_->global_operator()->Init();
  op3_preconditioner_->Setup(D1ptr, s_with_face, Teuchos::null);
  op3_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op3_preconditioner_->ApplyBCs(true, true);
  //op3_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op3_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "Analytic Jacobian A_13: " << *op3_preconditioner_->global_operator()->A() << "\n";
}


void CompW_PK::NumericalJacobian(double t_old, double t_new, 
                                Teuchos::RCP<const TreeVector> u, double eps)
{
  for (op_iter op_it = ops_.begin(); op_it != ops_.end(); op_it++) {
    (*op_it)->global_operator()->Init();
  }

  //std::cout << "solution: \n";
  //u->Print(std::cout);
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
  //std::cout << "f_ref: " << *f_ref->Data()->ViewComponent("cell") << "\n";

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
    //std::cout << "h: " << h << "\n";
    var_diff_c[0][c] += h;
    Functional(t_old, t_new, u_copy, u_diff, f_diff);
    //std::cout << "f_diff: " << *f_diff->Data()->ViewComponent("cell") << "\n";
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