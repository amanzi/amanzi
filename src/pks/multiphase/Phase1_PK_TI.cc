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
#include "Phase1_PK.hh"
#include "Epetra_Vector.h"

namespace Amanzi {
namespace Multiphase {


/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
This is basically the residual
****************************************************************** */
void Phase1_PK::Functional(double t_old, double t_new, 
              Teuchos::RCP<TreeVector> u_old,
              Teuchos::RCP<TreeVector> u_new, 
              Teuchos::RCP<TreeVector> f)
{
  double dTp(t_new - t_old);
  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_phase1 = u_new->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_phase2 = u_new->SubVector(1)->Data();
  //std::cout << "saturation phase2: " << *saturation_phase2->ViewComponent("cell") << "\n";
  //std::cout << "pressure phase1: " << *pressure_phase1->ViewComponent("cell") << "\n";

  // Compute saturation phase1
  *saturation_phase1_ = *saturation_phase2;
  saturation_phase1_->Scale(-1.0);
  saturation_phase1_->Shift(1.0);

  // Calculate relative perm needed to initialize diffusion operator
  rel_perm_w_->Compute(*saturation_phase1_);
  //rel_perm_w_->Krel()->Scale(rho_/mu_);
  //rel_perm_w_->dKdS()->Scale(-rho_/mu_);
  upwind_velocity_ = S_->GetFieldData("phase1_velocity", passwd_);
  RelativePermeabilityUpwindFn func1 = &RelativePermeability::Value;
  upwind_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_s_, bc_value_s_, *rel_perm_w_->Krel(), *rel_perm_w_->Krel(), func1);
  rel_perm_w_->Krel()->Scale(rho_/mu_);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  op1_matrix_->global_operator()->Init();
  op1_matrix_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_matrix_->ApplyBCs(true, true); // arguments have no effect for Finite Volume
  //op1_matrix_->global_operator()->SymbolicAssembleMatrix();
  //op1_matrix_->global_operator()->AssembleMatrix();
  //std::cout << "Phase1_PK::Functional matrix: " << *op1_matrix_->global_operator()->A() << "\n";

  // Update source term
  Teuchos::RCP<CompositeVector> rhs = op1_matrix_->global_operator()->rhs();
  if (src_sink_ != NULL) AddSourceTerms(*rhs);
  
  //op1_matrix_->global_operator()->SymbolicAssembleMatrix();
  //op1_matrix_->global_operator()->AssembleMatrix();
  //std::cout << "op_matrix phase1: \n";
  //std::cout << *op1_matrix_->global_operator()->A();
  op1_matrix_->global_operator()->ComputeNegativeResidual(*u_new->SubVector(0)->Data(), *f->Data());
  //std::cout << "Phase1 Functional: \n";
  //f->Print(std::cout);
  
  // Add time derivative
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
    f_cell[0][c] += -(s1 - s2) * factor;
  }
  //std::cout << "Phase1_PK residual: " << f_cell << "\n";
}

/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Phase1_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                    Teuchos::RCP<TreeVector> Pu) 
{
  return 0;
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void Phase1_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  // Get the new pressure and saturation from the solution tree vector
  Teuchos::RCP<const CompositeVector> pressure_phase1 = u->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> saturation_phase2 = u->SubVector(1)->Data();
  //std::cout << "saturation phase2: " << *saturation_phase2->ViewComponent("cell") << "\n";
  //std::cout << "pressure phase1: " << *pressure_phase1->ViewComponent("cell") << "\n";

  // compute saturation phase1
  *saturation_phase1_ = *saturation_phase2; 
  saturation_phase1_->Scale(-1.0);
  saturation_phase1_->Shift(1.0);

  // Calculate relative perm needed to initialize diffusion operator
  rel_perm_w_->Compute(*saturation_phase1_);
  upwind_velocity_ = S_->GetFieldData("phase1_velocity", passwd_);
  //std::cout << "upwind_velocity_: " << *upwind_velocity_->ViewComponent("face") << "\n";
  RelativePermeabilityUpwindFn func1 = &RelativePermeability::Value;
  upwind_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_s_, bc_value_s_, *rel_perm_w_->Krel(), *rel_perm_w_->Krel(), func1);
  rel_perm_w_->Krel()->Scale(rho_/mu_); 
  //std::cout << "Krel: " << *rel_perm_w_->Krel()->ViewComponent("face") << "\n";

  RelativePermeabilityUpwindFn func2 = &RelativePermeability::Derivative;
  upwind_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_s_, bc_value_s_, *rel_perm_w_->dKdS(), *rel_perm_w_->dKdS(), func2);
  //std::cout << "dKdS face: " << *rel_perm_w_->dKdS()->ViewComponent("face") << "\n";
  rel_perm_w_->dKdS()->Scale(-rho_/mu_); // negative sign since s2 is the primary variable, but dKdS is computed wrt s1.

  // the flux is only used for finding the upwind cells
  tmp_flux_->PutScalar(0.0);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  op1_matrix_->global_operator()->Init();
  op1_matrix_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_matrix_->UpdateFlux(*pressure_phase1, *tmp_flux_);

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

  op1_preconditioner_->global_operator()->Init();
  op1_preconditioner_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_preconditioner_->ApplyBCs(true, true);

  op2_preconditioner_->global_operator()->Init();
  op2_preconditioner_->Setup(*tmp_flux_);
  op2_preconditioner_->UpdateMatrices(*rel_perm_w_->dKdS());
  op2_preconditioner_->ApplyBCs(op_bc_s_, true);
  op2_preconditioner_->global_operator()->Rescale(-1.0);

  // create operator for accumulation term
  //op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op2_preconditioner_->global_operator()));

  // Create cvs for accumulation term 
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  
  CompositeVector porosity(*S_->GetFieldData("porosity"));
  //porosity.PutScalar(-phi_); 
  porosity.Scale(-rho_);
  if (dTp > 0.0) {
    op_acc_->AddAccumulationTerm(*saturation_phase2, porosity, dTp, "cell");
  } 
  
  //op1_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op1_preconditioner_->global_operator()->AssembleMatrix();
  //op2_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  //op2_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "Phase1_PK, prec matrix1: " << *op1_preconditioner_->global_operator()->A() << "\n";
  //std::cout << "Phase1_PK, prec matrix2: " << *op2_preconditioner_->global_operator()->A() << "\n"; 
  
}


void Phase1_PK::NumericalJacobian(double t_old, double t_new, 
                                  Teuchos::RCP<const TreeVector> u, double eps)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  //std::cout << "saturation phase2: " << *u->SubVector(1)->Data()->ViewComponent("cell") << "\n";
  //std::cout << "pressure phase1: " << *u->SubVector(0)->Data()->ViewComponent("cell") << "\n";
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  for (int ii = 0; ii < ops_.size(); ii++) {
    ops_[ii]->global_operator()->Init();
    ops_[ii]->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
    ops_[ii]->UpdateMatrices(Teuchos::null, Teuchos::null);
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
      //std::cout << "MyPID: " << MyPID << "; nfaces_owned: " << nfaces_owned << "; number of local matrices: " << ops_[ii]->local_matrices()->matrices.size() << "\n";
      int nfaces_none = 0;
      for (int f_it = 0; f_it < faces.size(); ++f_it) {
        int f_id = faces[f_it];
        //std::cout << "MyPID: " << MyPID << "; cell: " << c << "; faces for this cell: " << f_id << "\n"; 
        if (bc_model_p_[f_id] != Operators::OPERATOR_BC_NEUMANN) nfaces_none++;
      }
      for (int f_it = 0; f_it < faces.size(); ++f_it) {
        int f_id = faces[f_it];
        //if (f_id >= nfaces_owned) continue;
        AmanziMesh::Entity_ID_List cells;
        mesh_->face_get_cells(f_id, AmanziMesh::Parallel_type::ALL, &cells);
        int ncells = cells.size();
        //std::cout << "MyPID: " << MyPID << "; Face: " << f_id << "; ncells: " << ncells << "\n";

        WhetStone::DenseMatrix Aface(ncells, ncells);
        Aface = 0.0;

        //if (bc_model_p_[f_id] != Operators::OPERATOR_BC_NEUMANN)
        //{
        for (int i = 0; i != ncells; ++i) {
          //std::cout << "MyPID: " << MyPID << "; adjacent cells: " << cells[i] << "\n";
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
        WhetStone::DenseMatrix& tmp_matrix = ops_[ii]->local_matrices()->matrices[f_id];
        //std::cout << "MyPID: " << MyPID << "; face: "<< f_id << "; size of Aface: " << Aface.NumRows() << "; size of local matrix: " << tmp_matrix.NumRows() << "\n";
        //std::cout << "MyPID: " << MyPID << "; face: "<< f_id << "; size of Aface: " << Aface.NumCols() << "; size of local matrix: " << tmp_matrix.NumCols() << "\n";
        ASSERT(Aface.NumRows() == tmp_matrix.NumRows() && Aface.NumCols() == tmp_matrix.NumCols());
        tmp_matrix += Aface;
        //}
      }
    }

    ops_[ii]->ApplyBCs(true, true);
    //ops_[ii]->global_operator()->SymbolicAssembleMatrix();
    //ops_[ii]->global_operator()->AssembleMatrix();
    //std::cout << "Numerical jacobian " << ii << ": " << *ops_[ii]->global_operator()->A() << "\n";
  }
}

} // End namespace Amanzi
} // End namespace Flow
