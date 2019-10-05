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
#include "TotalPhase_PK.hh"
#include "Epetra_Vector.h"
#include "Op.hh"

namespace Amanzi {
namespace Multiphase {

class Op;

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
This is basically the residual
****************************************************************** */
void TotalPhase_PK::Functional(double t_old, double t_new, 
              Teuchos::RCP<TreeVector> u_old,
              Teuchos::RCP<TreeVector> u_new, 
              Teuchos::RCP<TreeVector> f)
{
  double dTp(t_new - t_old);
  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_phase1 = u_new->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_phase2 = u_new->SubVector(1)->Data();

  // compute saturation phase1 from saturation phase2
  Teuchos::RCP<CompositeVector> saturation_phase1 = Teuchos::rcp(new CompositeVector(*saturation_phase2));
  saturation_phase1->PutScalarMasterAndGhosted(1.0);
  saturation_phase1->Update(-1.0, *saturation_phase2, 1.0);

  /* Compute pressure of phase2 */
  //pressure_phase2_->Update(1.0, *pressure_phase1, 0.0); // set P2 = P1
  if (include_capillary_) {
    capillary_pressure_->Compute(*saturation_phase1);
    //pressure_phase2_->Update(1.0, *capillary_pressure_->Pc(), 1.0); // P2 = P1 + Pc    
  }

  // Calculate relative perm needed to initialize diffusion operator
  rel_perm_w_->Compute(*saturation_phase1);
  upwind_velocity_ = S_->GetFieldData("phase1_velocity", passwd_);
  RelativePermeabilityUpwindFn func1 = &RelativePermeability::Value;
  upwind_vw_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_, bc_value_s_, *rel_perm_w_->Krel(), *rel_perm_w_->Krel(), func1);
  rel_perm_w_->Krel()->Scale(1.0/mu1_);

  rel_perm_n_->Compute(*saturation_phase1);
  upwind_velocity_ = S_->GetFieldData("phase1_velocity", passwd_);
  upwind_vn_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_, bc_value_s_, *rel_perm_n_->Krel(), *rel_perm_n_->Krel(), func1);
  rel_perm_n_->Krel()->Scale(1.0/mu2_);
  //std::cout << "rel_perm_n_: \n";
  //rel_perm_n_->Krel()->Print(std::cout);
  tot_mobility_->Update(1.0, *rel_perm_w_->Krel(), 1.0, *rel_perm_n_->Krel(), 0.0);
  //std::cout << "total mobility: \n";
  //tot_mobility_->Print(std::cout);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  op1_matrix_->global_operator()->Init();
  op1_matrix_->Setup(Kptr, tot_mobility_, Teuchos::null);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null); 
  //std::cout << "op1_matrix_ ApplyBCs\n";
  op1_matrix_->ApplyBCs(true, true);
  //std::cout << "end op1_matrix_ ApplyBCs\n";
  
  if (include_capillary_) {
    op2_matrix_->global_operator()->Init();
    op2_matrix_->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
    op2_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null); 
    op2_matrix_->ApplyBCs(true, true);  // Update source term
  }

  //Teuchos::RCP<CompositeVector> rhs = op1_matrix_->rhs();
  //if (src_sink_ != NULL) AddSourceTerms(*rhs);
  
  //std::cout << "pressure phase2: \n";
  //pressure_phase2_->Print(std::cout);
  //op1_matrix_->global_operator()->SymbolicAssembleMatrix();
  //op1_matrix_->global_operator()->AssembleMatrix();
  //std::cout << "op1_matrix: \n";
  //std::cout << *op1_matrix_->global_operator()->A();
  //std::cout << "rhs: \n";
  //op1_matrix_->global_operator()->rhs()->Print(std::cout);
  //op1_matrix_->global_operator()->A()->Apply(*pressure_phase2_->ViewComponent("cell"), *f->Data()->ViewComponent("cell"));
  //f->Data()->Update(-1.0, *op1_matrix_->global_operator()->rhs(), 1.0);
  Teuchos::RCP<CompositeVector> f1 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  Teuchos::RCP<CompositeVector> f2 = Teuchos::rcp(new CompositeVector(f->Data()->Map()));
  f1->PutScalar(0.0);
  f2->PutScalar(0.0);
  op1_matrix_->global_operator()->ComputeNegativeResidual(*pressure_phase1, *f1);
  //std::cout << "TotalPhase_PK Functional f1: \n";
  //f1->Print(std::cout);
  if (include_capillary_) op2_matrix_->global_operator()->ComputeNegativeResidual(*capillary_pressure_->Pc(), *f2);
  f->Data()->Update(1.0, *f1, 1.0, *f2, 0.0);


  //std::cout << "TotalPhase_PK residual before adding accumulation term: " << *f->Data()->ViewComponent("cell") << "\n";
  
  // Add time derivative
  /*
  Epetra_MultiVector& f_cell = *f->Data()->ViewComponent("cell", true);

  const Epetra_MultiVector& S1_cell = *u_new->SubVector(1)->Data()->ViewComponent("cell", true);
  const Epetra_MultiVector& S2_cell = *u_old->SubVector(1)->Data()->ViewComponent("cell", true);
  const Epetra_MultiVector& por_cell = *S_->GetFieldData("porosity")->ViewComponent("cell", true);

  // Add accumulation term
  double s1, s2, phi, volume;
  for (int c = 0; c < f_cell.MyLength(); c++) {
    s1 = S1_cell[0][c];
    s2 = S2_cell[0][c];
    phi = por_cell[0][c];

    double factor = rho_ * phi * mesh_->cell_volume(c) / dTp;
    f_cell[0][c] += (s1 - s2) * factor;
  }
  */
  //std::cout << "TotalPhase_PK residual: \n";
  //f->Data()->Print(std::cout);
}

/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int TotalPhase_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                    Teuchos::RCP<TreeVector> Pu) 
{
  return 0;
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void TotalPhase_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_phase1 = u->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_phase2 = u->SubVector(1)->Data();
  //std::cout << "Saturation 2: " << endl;
  //saturation_phase2->Print(std::cout);

  // compute saturation phase1
  Teuchos::RCP<CompositeVector> phase1_sat = Teuchos::rcp(new CompositeVector(*saturation_phase2));
  phase1_sat->Scale(-1.0);
  phase1_sat->Shift(1.0);

  /*
  // Create data structures to initialize operators
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  //cvs.AddComponent("face", AmanziMesh::FACE, 1);

  Teuchos::RCP<CompositeVector> Krel = Teuchos::rcp(new CompositeVector(cvs));
  Krel->PutScalarMasterAndGhosted(0.0);
  */

  // Calculate total mobility needed to initialize diffusion operator
  rel_perm_w_->Compute(*phase1_sat);
  upwind_velocity_ = S_->GetFieldData("phase1_velocity", passwd_);
  RelativePermeabilityUpwindFn func1 = &RelativePermeability::Value;
  upwind_vw_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_, bc_value_s_, *rel_perm_w_->Krel(), *rel_perm_w_->Krel(), func1);
  rel_perm_w_->Krel()->Scale(1.0/mu1_); 

  rel_perm_n_->Compute(*phase1_sat);
  upwind_velocity_ = S_->GetFieldData("phase2_velocity", passwd_);
  upwind_vn_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_, bc_value_s_, *rel_perm_n_->Krel(), *rel_perm_n_->Krel(), func1);
  rel_perm_n_->Krel()->Scale(1.0/mu2_);

  tot_mobility_->Update(1.0, *rel_perm_w_->Krel(), 1.0, *rel_perm_n_->Krel(), 0.0);

  RelativePermeabilityUpwindFn func2 = &RelativePermeability::Derivative;
  upwind_vn_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_, bc_value_s_, *rel_perm_n_->dKdS(), *rel_perm_n_->dKdS(), func2);
  rel_perm_n_->dKdS()->Scale(-1.0/mu2_); // Must scale derivative by -1

  // Compute capillary pressure and its derivative
  if (include_capillary_) {
    capillary_pressure_->Compute(*phase1_sat);
    capillary_pressure_->dPc_dS()->Scale(-1.0); // Must scale derivative by -1
    //std::cout << "dPc_dS cell: " << *capillary_pressure_->dPc_dS()->ViewComponent("cell") << "\n";
    CapillaryPressureUpwindFn func3 = &CapillaryPressure::Derivative;
    upwind_pc_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_, bc_value_s_, 
                        *capillary_pressure_->dPc_dS(), *capillary_pressure_->dPc_dS(), func3);
  }

  // compute pressure phase2 from pressure phase1
  pressure_phase2_->Update(1.0, *pressure_phase1, 0.0); // P2 = P1
  if (include_capillary_) pressure_phase2_->Update(1.0, *capillary_pressure_->Pc(), 1.0); // P2 = P1 + Pc
  //std::cout << "Pressure 2: " << endl;
  //pressure_phase2_->Print(std::cout);

  tmp_flux_->PutScalar(0.0);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  op_matrix_copy1_->global_operator()->Init();
  op_matrix_copy1_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op_matrix_copy1_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_matrix_copy1_->UpdateFlux(*pressure_phase1, *tmp_flux_);

  // compute K * dKdS * grad(P1)
  Epetra_MultiVector& fface = *tmp_flux_->ViewComponent("face", true);
  Epetra_MultiVector& dKW_dS_f = *rel_perm_w_->dKdS()->ViewComponent("face", true);
  //std::cout << "dKdS: " << dKW_dS_f << endl;
  //std::cout << "grad P2: " << fface << endl;
  ASSERT(dKW_dS_f.MyLength() == fface.MyLength());
  for (int f = 0; f < dKW_dS_f.MyLength(); f++) {
    dKW_dS_f[0][f] *= fface[0][f]; 
  }
  //std::cout << "adv_coef: " << dKW_dS_f << endl;

  // Init operators
  // diffusion operator wrt P1
  op1_preconditioner_->global_operator()->Init();
  op1_preconditioner_->Setup(Kptr, tot_mobility_, Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_preconditioner_->ApplyBCs(true, true);

  // Advection operator wrt S2
  op2_preconditioner_->global_operator()->Init();
  op2_preconditioner_->Setup(*tmp_flux_);
  op2_preconditioner_->UpdateMatrices(*rel_perm_w_->dKdS());
  //op2_preconditioner_->ApplyBCs(op_bc_s_, true);
  op2_preconditioner_->global_operator()->Rescale(-1.0);

  tmp_flux_->PutScalar(0.0);
  op_matrix_copy1_->global_operator()->Init();
  op_matrix_copy1_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op_matrix_copy1_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_matrix_copy1_->UpdateFlux(*pressure_phase2_, *tmp_flux_);

  // compute K * dKdS * grad(P2)
  fface = *tmp_flux_->ViewComponent("face", true);
  Epetra_MultiVector& dKdS_f = *rel_perm_n_->dKdS()->ViewComponent("face", true);
  //std::cout << "dKdS: " << dKdS_f << endl;
  //std::cout << "grad P2: " << fface << endl;
  ASSERT(dKdS_f.MyLength() == fface.MyLength());
  for (int f = 0; f < dKdS_f.MyLength(); f++) {
    dKdS_f[0][f] *= fface[0][f]; 
  }

  //op3_preconditioner_->global_operator()->Init();
  op3_preconditioner_->Setup(*tmp_flux_);
  op3_preconditioner_->UpdateMatrices(*rel_perm_n_->dKdS());
  op3_preconditioner_->ApplyBCs(op_bc_s_, true);
  //op2_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op2_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "advection matrix alone: \n" << *op2_preconditioner_->global_operator()->A();
  double norm_adv = 0.0;
  rel_perm_n_->dKdS()->NormInf(&norm_adv);

  // norms of diffusion coefficient and velocity
  double pec_inf, pec_1, pec_2;
  CompositeVector peclet(*rel_perm_n_->dKdS());
  Epetra_MultiVector& peclet_f = *peclet.ViewComponent("face", true);
  for (int f = 0; f < peclet_f.MyLength(); f++) {
    peclet_f[0][f] *= mesh_->face_area(f);
  }

  if (include_capillary_)
  {
    tmp_flux_->PutScalar(0.0);
    //op_matrix_->SetBCs(op_bc_pc_prime_, op_bc_pc_prime_);
    Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
    op_matrix_copy_->global_operator()->Init();
    op_matrix_copy_->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
    op_matrix_copy_->UpdateMatrices(Teuchos::null, Teuchos::null);
    op_matrix_copy_->UpdateFlux(*capillary_pressure_->dPc_dS(), *tmp_flux_);
    //std::cout << "K * Krel * grad(Pc') \n";
    //tmp_flux_->Print(std::cout);
     
    // additional advection operator for K * lambda_n * grad(Pc')
    //Teuchos::ParameterList olist_adv = mp_list_.sublist("operators").sublist("advection operator");
    //op_sum1_ = Teuchos::rcp(new Operators::OperatorAdvection(olist_adv, op2_preconditioner_->global_operator()));
    op_sum1_->Setup(*tmp_flux_);
    op_sum1_->UpdateMatrices(*tmp_flux_);
    op_sum1_->ApplyBCs(op_bc_s_, true);
    //op_sum1_->global_operator()->Rescale(-1.0);
    //op_sum1_->global_operator()->SymbolicAssembleMatrix();
    //op_sum1_->global_operator()->AssembleMatrix();
    //std::cout << "advection matrix + lambda_n * grad(Pc'): \n" << *op_sum1_->global_operator()->A();

    // Compute the coefficient for the diffusion operator wrt s2
    //std::cout << "Phase2_PK_TI: Krel: " << *rel_perm_n_->Krel()->ViewComponent("face") << "\n";
    //std::cout << "Phase2_PK_TI: dPc_dS: " << *capillary_pressure_->dPc_dS()->ViewComponent("face") << "\n";
    rel_perm_n_->Krel()->Multiply(1.0, *rel_perm_n_->Krel(), *capillary_pressure_->dPc_dS(), 0.0);
    //std::cout << "Phase2_PK_TI: Krel * dPc_dS: " << *rel_perm_n_->Krel()->ViewComponent("face") << "\n";

    // compute the grid peclet number
    //peclet.ReciprocalMultiply(1.0, *rel_perm_n_->Krel(), peclet, 0.0);
    
    Epetra_MultiVector& diff_f = *rel_perm_n_->Krel()->ViewComponent("face", true);
    Epetra_MultiVector& adv_f = *peclet.ViewComponent("face", true);
    Epetra_MultiVector& pec_f = *peclet.ViewComponent("face", true);
    for (int f = 0; f < diff_f.MyLength(); f++) {
      //std::cout << "diff_f: " << diff_f[0][f] << "; adv_f: " << adv_f[0][f] << endl;
      if (std::abs(diff_f[0][f]) > 1e-12) pec_f[0][f] = adv_f[0][f] / diff_f[0][f];
      else pec_f[0][f] = 0.0;
      //std::cout << "peclet: " << pec_f[0][f] << endl;
    }
    
    //std::cout << "Peclet number: " << endl;
    //peclet.Print(std::cout);
    //pec_f.NormInf(&pec_inf);
    //pec_f.Norm2(&pec_2);
    //pec_f.Norm1(&pec_1);
    //std::cout << "pec_inf: " << pec_inf << "; pec_1: " << pec_1 << "; pec_2: " << pec_2 << endl;

    // diffusion operator wrt S2, corresponding to the term K * lambda_n * Pc' * grad
    //Teuchos::ParameterList& op_list = mp_list_.sublist("operators").sublist("diffusion operator").sublist("preconditioner1");
    //op_sum_ = Teuchos::rcp(new Operators::OperatorDiffusionFV(op_list, op_sum1_->global_operator()));
    //op_sum_->SetGravity(gravity_);
    op_sum_->SetBCs(op_bc_s_, op_bc_s_);
    op_sum_->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
    op_sum_->UpdateMatrices(Teuchos::null, Teuchos::null);
    op_sum_->ApplyBCs(true, true);
    //op_sum_->global_operator()->SymbolicAssembleMatrix();
    //op_sum_->global_operator()->AssembleMatrix();
    //std::cout << "advection matrix + lambda_n * grad(Pc') + diffsion wrt s2: \n" << *op_sum_->global_operator()->A();
  }

  // create accumulation term
  // op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op2_preconditioner_->global_operator())); 

  /*
  CompositeVector porosity(*S_->GetFieldData("porosity"));
  porosity.Scale(rho_);
  //porosity.PutScalar(phi_);
  //porosity.Scale(rho_);
  if (dTp > 0.0) {
    op_acc_->AddAccumulationTerm(*saturation_phase2, porosity, dTp, "cell");
  } 
  */

  //op1_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op1_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "TotalPhase_PK, prec matrix1: " << *op1_preconditioner_->global_operator()->A() << "\n";
  //std::cout << "TotalPhase_PK, prec matrix2: " << *op2_preconditioner_->global_operator()->A() << "\n";
  //op_sum_->global_operator()->SymbolicAssembleMatrix();
  //op_sum_->global_operator()->AssembleMatrix();
  //std::cout << "TotalPhase_PK, prec sum matrix: " << *op_sum_->global_operator()->A() << "\n";
} // End UpdatePreconditioner


void TotalPhase_PK::NumericalJacobian(double t_old, double t_new, 
                                  Teuchos::RCP<const TreeVector> u, double eps)
{
  //std::cout << "saturation phase1: " << *u->SubVector(1)->Data()->ViewComponent("cell") << "\n";
  //std::cout << "pressure phase1: " << *u->SubVector(0)->Data()->ViewComponent("cell") << "\n";
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

  for (int ii = 0; ii < ops_.size(); ii++)
  {
    local_op_iter local_op_it = ops_[ii]->global_operator()->OpBegin();
    Epetra_MultiVector& var_diff_c = *u_diff->SubVector(ii)->Data()->ViewComponent("cell", true);
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
      Epetra_MultiVector& deriv_c = *deriv->ViewComponent("cell", true);
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
        //std::cout << "Face: " << f_id << "; ncells: " << ncells << "\n";

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
        ASSERT(Aface.NumRows() == tmp_matrix.NumRows() && Aface.NumCols() == tmp_matrix.NumCols());
        tmp_matrix += Aface;
        }
      }
    }

    ops_[ii]->ApplyBCs(true, true);
    //ops_[i]->global_operator()->SymbolicAssembleMatrix();
    //ops_[i]->global_operator()->AssembleMatrix();
    //std::cout << "Numerical jacobian" << i << ": " << *ops_[i]->global_operator()->A() << "\n";
  }
}

} // End namespace Amanzi
} // End namespace Flow
