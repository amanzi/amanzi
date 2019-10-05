/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)

  The routine implements interface to the BDF1 time integrator.  
*/

#include <algorithm>
#include <string>
#include <vector>

#include "LinearOperatorFactory.hh"
#include "Saturation_PK.hh"
#include "OperatorDefs.hh"
#include "Op.hh"

namespace Amanzi {
namespace Multiphase {

class Op;

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void Saturation_PK::Functional(double Told, double Tnew, 
                             Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
                             Teuchos::RCP<TreeVector> f)
{ 
  double Tp(Tnew), dTp(Tnew - Told);
  
  // Compute the fractional flow function with new solution
  Teuchos::RCP<const CompositeVector> water_saturation = u_new->Data();
  frac_flow_->Compute(*water_saturation);
  fractional_flow_ = frac_flow_->Frac_Flow();

  // Now upwind fractional flow
  FractionalFlowUpwindFn func1 = &FractionalFlow::Value; 
  upwind_->Compute(*darcy_flux_, *darcy_flux_, bc_model, bc_value, *fractional_flow_, *fractional_flow_, func1);

  // create the diffusion coefficient for capillary pressure
  Teuchos::RCP<CompositeVector> diff_coef = Teuchos::rcp(new CompositeVector(*rel_perm_n_->Krel()));
  diff_coef->PutScalar(0.0);
  if (include_capillary_) 
  {
    capillary_pressure_->Compute(*water_saturation);
    //std::cout << "capillary_pressure_: " << *capillary_pressure_->Pc()->ViewComponent("cell") << "\n";
    rel_perm_n_->Compute(*water_saturation);
    RelativePermeabilityUpwindFn func1 = &RelativePermeability::Value;
    upwind_n_->Compute(*darcy_flux_, *darcy_flux_, bc_model, bc_value, *rel_perm_n_->Krel(), *rel_perm_n_->Krel(), func1);
    rel_perm_n_->Krel()->Scale(1.0/mu_[1]);
    diff_coef->Multiply(1.0, *rel_perm_n_->Krel(), *fractional_flow_, 0.0);
    //diff_coef->Multiply(1.0, *diff_coef, *capillary_pressure_->dPc_dS(), 0.0);
    //std::cout << "dPc_dS: " << *capillary_pressure_->dPc_dS()->ViewComponent("cell") << "\n";
    //std::cout << "diff_coef: " << *diff_coef->ViewComponent("face") << "\n";
  }

  // Multiply fractional flow by the velocity/flux computed in Pressure_PK
  fractional_flow_->Multiply(1.0, *fractional_flow_, *darcy_flux_, 0.0);
  //std::cout << "fractional_flow_: " << *fractional_flow_->ViewComponent("face") << "\n";
  /*
  Epetra_MultiVector& darcy_flux_f = *fractional_flow_->ViewComponent("face");
  Epetra_MultiVector& fflow_f = *fractional_flow_->ViewComponent("face");
  for (int nf = 0; nf < fflow_f.MyLength(); nf++) {
    fflow_f[0][nf] *= darcy_flux_f[0][nf];
  }
  */

  // Initialize advection operator (need to do this again when compute functional and update preconditioner)
  // the assumption is that darcy_flux_ and fractional_flow_ have been changed, so we need to update the operators.
  op_matrix_->global_operator()->Init();
  op_matrix_->Setup(*darcy_flux_);
  op_matrix_->UpdateMatrices(*fractional_flow_);
  op_matrix_->ApplyBCs(op_bc_, true);
  op_matrix_->global_operator()->SymbolicAssembleMatrix();
  op_matrix_->global_operator()->AssembleMatrix();
  //std::cout << "Functional matrix: " << *op_matrix_->global_operator()->A() << "\n";
  //std::cout << "Functional rhs: " << *op_matrix_->global_operator()->rhs()->ViewComponent("cell") << "\n";

  // create diffusion operator if capillary pressure is present
  if (include_capillary_)
  {
    Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
    op1_matrix_->global_operator()->Init();
    //op1_matrix_->SetDensity(rho_[1] - rho_[0]);
    op1_matrix_->Setup(Kptr, diff_coef, Teuchos::null);
    op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
    op1_matrix_->ApplyBCs(true, true);
    //op1_matrix_->global_operator()->SymbolicAssembleMatrix();
    //op1_matrix_->global_operator()->AssembleMatrix();
    //std::cout << "Functional matrix: " << *op1_matrix_->global_operator()->A() << "\n";
    //std::cout << "Functional rhs: " << *op1_matrix_->global_operator()->rhs()->ViewComponent("cell") << "\n";
  }

  // Add source term if any
  Teuchos::RCP<CompositeVector> rhs = op_matrix_->global_operator()->rhs();
  // std::cout << "op_matrix_ rhs: " << *rhs->ViewComponent("cell") << "\n";
  if (src_sink_ != NULL) {
    AddSourceTerm(*rhs, *fractional_flow_);
  }

  // create a unit vector
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);

  CompositeVector u_unit(cvs);
  u_unit.PutScalar(1.0);

  op_matrix_->global_operator()->ComputeNegativeResidual(u_unit, *f->Data());
  //std::cout << "Functional Residual before adding accumulation term: " << *f->Data()->ViewComponent("cell") << "\n";
  if (include_capillary_) {
    TreeVector f_tmp(*f);
    f_tmp.PutScalar(0.0);
    op1_matrix_->global_operator()->ComputeNegativeResidual(*capillary_pressure_->Pc(), *f_tmp.Data(), false);
    // subtract the residual from the diffusion operator since the default include the negative sign
    f->Update(-1.0, f_tmp, 1.0);
  }
  //std::cout << "Functional Residual before adding accumulation term: " << *f->Data()->ViewComponent("cell") << "\n";

  // Add time derivative
  Epetra_MultiVector& f_cell = *f->Data()->ViewComponent("cell", true);

  const Epetra_MultiVector& S1_cell = *u_new->Data()->ViewComponent("cell");
  const Epetra_MultiVector& S2_cell = *u_old->Data()->ViewComponent("cell");
  const Epetra_MultiVector& phi_cell = *S_->GetFieldData("porosity")->ViewComponent("cell");

  // Add accumulation term
  double s1, s2, phi, volume;
  for (int c = 0; c < f_cell.MyLength(); c++) {
    s1 = S1_cell[0][c];
    s2 = S2_cell[0][c];
    phi = phi_cell[0][c];

    //double factor = rho_[0] * phi_ * mesh_->cell_volume(c) / dTp;
    // note that we assume density is constant so we divide the equation 
    // by rho on both sides.
    double factor = phi * mesh_->cell_volume(c) / dTp;
    f_cell[0][c] += (s1 - s2) * factor;
  }
  //std::cout << "Functional matrix A: " << *op_matrix_->A_ << "\n";
  //std::cout << "Functional RHS: " << *op_matrix_->rhs_->ViewComponent("cell") << "\n";
  //std::cout << "Functional u_new: " << *u_new->Data()->ViewComponent("cell") << "\n";
  //std::cout << "Functional u_old: " << *u_old->Data()->ViewComponent("cell") << "\n";
  //std::cout << "Functional Residual: " << *f->Data()->ViewComponent("cell") << "\n";

}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Saturation_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                      Teuchos::RCP<TreeVector> pu)
{
  int ierr = 0;
  Amanzi::timer_manager.start("ApplyPreconditioner");
  if (jacobian_type_ == "numerical")
    solver = sfactory.Create(ti_specs_->solver_name, *linear_operator_list_, op1_preconditioner_->global_operator());
  else
    solver = sfactory.Create(ti_specs_->solver_name, *linear_operator_list_, op_preconditioner_->global_operator());

  ierr = solver->ApplyInverse(*u->Data(), *pu->Data());
  ls_itrs_ += solver->num_itrs();
  Amanzi::timer_manager.stop("ApplyPreconditioner");
  return ierr;
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void Saturation_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  //std::cout << "UpdatePreconditioner, u: " << *u->Data()->ViewComponent("cell") << "\n";
  if (jacobian_type_ == "numerical")
  {
    NumericalJacobian(Tp, u, dTp);
  }
  else {
    AnalyticJacobian(Tp, u, dTp);
  }
}


void Saturation_PK::NumericalJacobian(double Tp, Teuchos::RCP<const TreeVector> u, double dTp) 
{
  double t_old, t_new, eps;
  t_old = 0.0;
  t_new = t_old + dTp;
  eps = 1e-12;
  op1_preconditioner_->global_operator()->Init();
  Teuchos::RCP<TreeVector> u_copy = Teuchos::rcp(new TreeVector(*u));

  // Create data structures to use functional
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false); 

  Teuchos::RCP<CompositeVector> deriv = Teuchos::rcp(new CompositeVector(cvs));
  deriv->PutScalar(0.0);

  // compute f(u_0)
  Teuchos::RCP<TreeVector> f_ref = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<CompositeVector> f_ref_cv = Teuchos::rcp(new CompositeVector(cvs));
  f_ref->SetData(f_ref_cv);
  Functional(t_old, t_new, u_copy, u_copy, f_ref);

  Teuchos::RCP<TreeVector> u_diff = Teuchos::rcp(new TreeVector(*u)); // u_diff = u_0 + du
  Teuchos::RCP<TreeVector> f_diff = Teuchos::rcp(new TreeVector()); // f(u_diff)
  Teuchos::RCP<CompositeVector> f_diff_cv = Teuchos::rcp(new CompositeVector(cvs));
  f_diff->SetData(f_diff_cv);

  Epetra_MultiVector& var_diff_c = *u_diff->Data()->ViewComponent("cell");
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
      if (bc_model[f_id] != Operators::OPERATOR_BC_NEUMANN) nfaces_none++;
    }
    for (int f_it = 0; f_it < faces.size(); ++f_it) {
      int f_id = faces[f_it];
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f_id, AmanziMesh::USED, &cells);
      int ncells = cells.size();
      //std::cout << "Face: " << f_id << "; bc type: " << bc_model[f_id] << "\n";

      //Epetra_MultiVector& deriv_c = *deriv->ViewComponent("cell");
      //std::cout << "numerical deriv: " << deriv_c << "\n";
      WhetStone::DenseMatrix Aface(ncells, ncells);
      Aface = 0.0;

      if (bc_model[f_id] != Operators::OPERATOR_BC_NEUMANN)
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
        //WhetStone::DenseMatrix& tmp_matrix = (*local_op_it)->matrices[f_id];
        WhetStone::DenseMatrix& tmp_matrix = op1_preconditioner_->local_matrices()->matrices[f_id];
        ASSERT(Aface.NumRows() == tmp_matrix.NumRows() && Aface.NumCols() == tmp_matrix.NumCols());
        //Aface += tmp_matrix;
        //tmp_matrix = Aface;
        tmp_matrix += Aface;
      }
    }
  }
  op1_preconditioner_->ApplyBCs(true, true);

  op1_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  op1_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "numerical jacobian: " << *op1_preconditioner_->global_operator()->A() << "\n";
  op1_preconditioner_->global_operator()->InitPreconditioner(ti_specs_->preconditioner_name, *preconditioner_list_);
}


void Saturation_PK::AnalyticJacobian(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  // recompute the fractional flow using the new saturation
  Teuchos::RCP<const CompositeVector> water_saturation = u->Data();
  //std::cout << "water_saturation " << *water_saturation->ViewComponent("cell") << "\n";
  frac_flow_->Compute(*water_saturation);
  fractional_flow_ = frac_flow_->Frac_Flow();
  dfw_dS_ = frac_flow_->dF_dS();
  //std::cout << "dfw_dS_ " << *dfw_dS_->ViewComponent("cell") << "\n";

  // Now upwind fractional flow
  FractionalFlowUpwindFn func1 = &FractionalFlow::Value;
  upwind_->Compute(*darcy_flux_, *darcy_flux_, bc_model, bc_value, *fractional_flow_, *fractional_flow_, func1);
  FractionalFlowUpwindFn func2 = &FractionalFlow::Derivative;
  upwind_->Compute(*darcy_flux_, *darcy_flux_, bc_model, bc_value, *dfw_dS_, *dfw_dS_, func2);
  //std::cout << "dfw_dS_ after upwind" << *dfw_dS_->ViewComponent("face") << "\n";
  Teuchos::RCP<CompositeVector> fractional_flow_copy_ = Teuchos::rcp(new CompositeVector(*fractional_flow_));
  //fractional_flow_->Multiply(1.0, *fractional_flow_, *darcy_flux_, 0.0);
  Teuchos::RCP<CompositeVector> dfw_dS_copy = Teuchos::rcp(new CompositeVector(*dfw_dS_));
  dfw_dS_->Multiply(1.0, *dfw_dS_, *darcy_flux_, 0.0);
  //std::cout << "darcy_flux_: " << *darcy_flux_->ViewComponent("face") << "\n";
  //std::cout << "dfw_dS_ face: " << *dfw_dS_->ViewComponent("face") << "\n";

  /*
  Epetra_MultiVector& darcy_flux_f = *fractional_flow_->ViewComponent("face");
  Epetra_MultiVector& fflow_f = *fractional_flow_->ViewComponent("face");
  Epetra_MultiVector& d_fflow_f = *fractional_flow_->ViewComponent("face");
  for (int nf = 0; nf < fflow_f.MyLength(); nf++) {
    fflow_f[0][nf] *= darcy_flux_f[0][nf];
    d_fflow_f[0][nf] *= darcy_flux_f[0][nf];
  }
  */

  //std::cout << "darcy_flux_face: " << *darcy_flux_->ViewComponent("face") << "\n";
  // update advection operator for preconditioner
  //op_->Init();
  op_preconditioner_->global_operator()->Init();
  op_preconditioner_->Setup(*darcy_flux_);
  op_preconditioner_->UpdateMatrices(*dfw_dS_);
  op_preconditioner_->ApplyBCs(op_bc_, true);

  op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_preconditioner_->global_operator()));

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  //cvs.AddComponent("face", AmanziMesh::FACE, 1);

  //CompositeVector porosity(cvs);
  //porosity.PutScalar(phi_);
  CompositeVector porosity(*S_->GetFieldData("porosity"));

  if (dTp > 0.0) {
    op_acc_->AddAccumulationTerm(*u->Data(), porosity, dTp, "cell");
  }

  // add source term
  Teuchos::RCP<CompositeVector> rhs = op_preconditioner_->global_operator()->rhs();
  if (src_sink_ != NULL) {
    AddSourceTerm(*rhs, *fractional_flow_);
  }

  // Compute capillary pressure and its derivative
  if (include_capillary_) {
    rel_perm_n_->Compute(*water_saturation);
    rel_perm_n_->dKdS()->Scale(-1.0); // must scale derivative by -1 since s2 is the primary variable
    RelativePermeabilityUpwindFn func1 = &RelativePermeability::Value;
    upwind_n_->Compute(*darcy_flux_, *darcy_flux_, bc_model, bc_value, *rel_perm_n_->Krel(), *rel_perm_n_->Krel(), func1);
    RelativePermeabilityUpwindFn func2 = &RelativePermeability::Derivative;
    upwind_n_->Compute(*darcy_flux_, *darcy_flux_, bc_model, bc_value, *rel_perm_n_->dKdS(), *rel_perm_n_->dKdS(), func2);
    rel_perm_n_->Krel()->Scale(1.0/mu_[1]);
    rel_perm_n_->dKdS()->Scale(-1.0/mu_[1]);

    capillary_pressure_->Compute(*water_saturation);
    capillary_pressure_->dPc_dS()->Scale(-1.0); // Must scale derivative by -1 since s2 is the primary variable
    //std::cout << "dPc_dS cell: " << *capillary_pressure_->dPc_dS()->ViewComponent("cell") << "\n";
    CapillaryPressureUpwindFn func3 = &CapillaryPressure::Derivative;
    upwind_pc_->Compute(*darcy_flux_, *darcy_flux_, bc_model, bc_value, 
                        *capillary_pressure_->dPc_dS(), *capillary_pressure_->dPc_dS(), func3);
    capillary_pressure_->dPc_dS()->Scale(-1.0);

    Teuchos::RCP<CompositeVector> diff_coef = Teuchos::rcp(new CompositeVector(*rel_perm_n_->Krel()));
    diff_coef->Multiply(1.0, *diff_coef, *fractional_flow_, 0.0);
    diff_coef->Multiply(1.0, *diff_coef, *capillary_pressure_->dPc_dS(), 0.0);

    diff_coef->Scale(-1.0); // scale by -1 since we add - div(grad Pc)
    //std::cout << "diff_coef: " << *diff_coef->ViewComponent("face") << "\n";
    //fractional_flow_copy_->Multiply(1.0, *fractional_flow_copy_, *rel_perm_n_->Krel(), 0.0);
    //fractional_flow_copy_->Multiply(1.0, *fractional_flow_copy_, *capillary_pressure_->dPc_dS(), 0.0);
    //fractional_flow_copy_->Scale(-1.0);

    Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
    Teuchos::ParameterList& op_list = mp_list_->sublist("operators").sublist("diffusion operator").sublist("preconditioner1");
    op_pres_pc_ = Teuchos::rcp(new Operators::OperatorDiffusionFV(op_list, op_preconditioner_->global_operator()));
    //op_pres_pc_->SetGravity(gravity_);
    op_pres_pc_->SetBCs(op_bc_, op_bc_);
    op_pres_pc_->Setup(Kptr, diff_coef, Teuchos::null);
    op_pres_pc_->UpdateMatrices(Teuchos::null, Teuchos::null);
    op_pres_pc_->ApplyBCs(true, true);

    Teuchos::RCP<CompositeVector> tmp_flux_ = Teuchos::rcp(new CompositeVector(*darcy_flux_));
    tmp_flux_->PutScalar(0.0);
    op1_matrix_->UpdateFlux(*capillary_pressure_->dPc_dS(), *tmp_flux_);

    Teuchos::ParameterList olist_adv = mp_list_->sublist("operators").sublist("advection operator");
    op_sum_ = Teuchos::rcp(new Operators::OperatorAdvection(olist_adv, op_pres_pc_->global_operator()));
    op_sum_->Setup(*tmp_flux_);
    op_sum_->UpdateMatrices(*tmp_flux_);
    op_sum_->ApplyBCs(op_bc_, true);

    Teuchos::RCP<CompositeVector> dlambda_fw_ds = Teuchos::rcp(new CompositeVector(*fractional_flow_));
    dlambda_fw_ds->Multiply(1.0, *dlambda_fw_ds, *rel_perm_n_->dKdS(), 0.0);
    dlambda_fw_ds->Multiply(1.0, *rel_perm_n_->Krel(), *dfw_dS_copy, 1.0);
    dlambda_fw_ds->Scale(-1.0);
    tmp_flux_->PutScalar(0.0);
    op1_matrix_->global_operator()->Init();
    //op1_matrix_->SetDensity(rho_[0] - rho_[1]);
    op1_matrix_->Setup(Kptr, dlambda_fw_ds, Teuchos::null);
    op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
    op1_matrix_->UpdateFlux(*capillary_pressure_->Pc(), *tmp_flux_);
    //dlambda_fw_ds->Multiply(1.0, *dlambda_fw_ds, *tmp_flux_, 0.0);
    //std::cout << "dlambda_fw_ds: " << *dlambda_fw_ds->ViewComponent("face") << "\n";

    op_sum1_ = Teuchos::rcp(new Operators::OperatorAdvection(olist_adv, op_sum_->global_operator()));
    op_sum1_->Setup(*tmp_flux_);
    op_sum1_->UpdateMatrices(*tmp_flux_);
    op_sum1_->ApplyBCs(op_bc_, true);
  }

  // finalize preconditioner
  op_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  op_preconditioner_->global_operator()->AssembleMatrix();
  //std::cout << "UpdatePreconditioner: u: " << *u->Data()->ViewComponent("cell") << "\n";
  //std::cout<<"Saturation preconditioner: Assemble Matrix\n";
  //std::cout<<*op_->A()<<"\n";
  //std::cout<<"UpdatePreconditioner: RHS\n";
  //std::cout<<*op_preconditioner_->rhs()->ViewComponent("cell")<<"\n";
  op_preconditioner_->global_operator()->InitPreconditioner(ti_specs_->preconditioner_name, *preconditioner_list_);
}

double Saturation_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) 
{
  double du_inf;
  const Epetra_MultiVector& du_cell = *du->Data()->ViewComponent("cell");

  du_cell.NormInf(&du_inf);
  return du_inf;
}


AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Saturation_PK::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
           Teuchos::RCP<const TreeVector> u,
           Teuchos::RCP<TreeVector> du)
{
  Teuchos::RCP<CompositeVector> s_next = Teuchos::rcp(new CompositeVector(*u->Data()));
  s_next->Update(-1.0, *du->Data(), 1.0);
  ClipSaturation(s_next, 0.0);
  du->Data()->Update(1.0, *u->Data(), -1.0, *s_next, 0.0);
  /*
  std::cout << "ModifyCorrection, value of h: " << h << "\n";
  int alpha = 0;
  int i = 0;
  bool fail = true;
  double norm_res, norm_res_next;
  res->NormInf(&norm_res);
  Teuchos::RCP<TreeVector> u_tmp = Teuchos::rcp(new TreeVector(*u));
  Teuchos::RCP<TreeVector> res_next = Teuchos::rcp(new TreeVector(*u));
  res_next->PutScalar(0.0);
  Teuchos::RCP<CompositeVector> s_next = u_tmp->SubVector(1)->Data();
  while (fail && i < 10){
    s_next->Update(-std::pow(2.0, alpha), *du->SubVector(1)->Data(), 1.0, *u->SubVector(1)->Data(), 0.0);
    ClipSaturation(s_next, 0.0);
    Functional(0.0, h, soln_, u_tmp, res_next);
    res_next->NormInf(&norm_res_next);

    alpha--;
    i++;
    fail = !((norm_res_next - norm_res) < 1e-12);
  }
  du->SubVector(1)->Data()->Update(1.0, *u->SubVector(1)->Data(), -1.0, *s_next, 0.0);
  */
}


void Saturation_PK::ClipSaturation(Teuchos::RCP<CompositeVector> s, double tol)
{
  double max_s, min_s;
  max_s = 1.0 - tol;
  min_s = 0.0 + tol;
  Epetra_MultiVector& s_c = *s->ViewComponent("cell");
  for (int c = 0; c < s_c.MyLength(); c++)
  {
    s_c[0][c] = std::min(std::max(min_s, s_c[0][c]), max_s);
  }
}


} // End namespace Multiphase
} // End namespace Amanzi
