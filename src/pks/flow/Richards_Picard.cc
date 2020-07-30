/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <vector>

#include "RemapUtils.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Calculates steady-state solution assuming that absolute and
* relative permeabilities do not depend explicitly on time.
* This is the experimental method.                                                 
****************************************************************** */
int Richards_PK::AdvanceToSteadyState_Picard(Teuchos::ParameterList& plist)
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();

  // create verbosity object
  VerboseObject* vo = new VerboseObject("Amanzi::Picard", *fp_list_); 

  CompositeVector  solution_old(*solution);
  CompositeVector& solution_new = *solution;
  CompositeVector  residual(*solution);

  Epetra_MultiVector& pold_cell = *solution_old.ViewComponent("cell");
  Epetra_MultiVector& pnew_cell = *solution_new.ViewComponent("cell");

  // update steady state boundary conditions
  double time = S_->time();
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(time, time);
    bcs_[i]->ComputeSubmodel(mesh_);
  }

  // update steady state source conditions
  for (int i = 0; i < srcs.size(); ++i) {
    srcs[i]->Compute(time, time); 
  }

  Teuchos::RCP<const CompositeVector> mu = S_->GetFieldData("viscosity_liquid");

  std::string linear_solver = plist.get<std::string>("linear solver");
  Teuchos::ParameterList& tmp_list = plist.sublist("picard parameters");
  int max_itrs_nonlinear = tmp_list.get<int>("maximum number of iterations", 400);
  double residual_tol_nonlinear = tmp_list.get<double>("convergence tolerance", 1e-8);

  int itrs = 0;
  double L2norm, L2error = 1.0;

  while (L2error > residual_tol_nonlinear && itrs < max_itrs_nonlinear) {
    // update dynamic boundary conditions
    for (int i = 0; i < bcs_.size(); i++) {
      if (bcs_[i]->bc_name() == "seepage") {
        bcs_[i]->Compute(time, time);
        bcs_[i]->ComputeSubmodel(mesh_);
      }
    }
    ComputeOperatorBCs(*solution);

    // update permeabilities
    darcy_flux_copy->ScatterMasterToGhosted("face");

    relperm_->Compute(solution, bc_model, bc_value, krel_);
    upwind_->Compute(*darcy_flux_copy, *solution, bc_model, *krel_);
    Operators::CellToFace_ScaleInverse(mu, krel_);
    krel_->ScaleMasterAndGhosted(molar_rho_);

    relperm_->ComputeDerivative(solution, bc_model, bc_value, dKdP_);
    upwind_->Compute(*darcy_flux_copy, *solution, bc_model, *dKdP_);
    Operators::CellToFace_ScaleInverse(mu, dKdP_);
    dKdP_->ScaleMasterAndGhosted(molar_rho_);

    // create algebraic problem (matrix = preconditioner)
    op_preconditioner_->Init();
    op_preconditioner_diff_->UpdateMatrices(darcy_flux_copy.ptr(), solution.ptr());
    op_preconditioner_diff_->UpdateMatricesNewtonCorrection(darcy_flux_copy.ptr(), Teuchos::null, molar_rho_);
    op_preconditioner_diff_->ApplyBCs(true, true, true);

    Teuchos::RCP<CompositeVector> rhs = op_preconditioner_->rhs();  // export RHS from the matrix class
    AddSourceTerms(*rhs);

    // create preconditioner
    op_preconditioner_->AssembleMatrix();
    op_preconditioner_->UpdatePreconditioner();

    // check convergence of non-linear residual
    op_preconditioner_->ComputeResidual(solution_new, residual);
    residual.Norm2(&L2error);
    rhs->Norm2(&L2norm);
    L2error /= L2norm;

    // solve linear problem
       solver = factory.Create(linear_solver, *linear_operator_list_, op_preconditioner_);

    solver->ApplyInverse(*rhs, *solution);

    int num_itrs_linear = solver->num_itrs();
    double linear_residual = solver->residual();

    // update relaxation
    double relaxation;
    relaxation = CalculateRelaxationFactor(pold_cell, pnew_cell);

    if (vo->getVerbLevel() >= Teuchos::VERB_HIGH) {
      Teuchos::OSTab tab = vo->getOSTab();
      *(vo->os()) << itrs << ": ||r||=" << L2error << " relax=" << relaxation 
                  << " lin_solver(" << linear_residual << ", " << num_itrs_linear << ")" << std::endl;
    }

    solution_new.Update(1.0 - relaxation, solution_old, relaxation);
    solution_old = solution_new;

    itrs++;
  }

  delete vo;
  return 0;
}


/* ******************************************************************
* Calculate relaxation factor.                                                       
****************************************************************** */
double Richards_PK::CalculateRelaxationFactor(const Epetra_MultiVector& uold,
                                              const Epetra_MultiVector& unew)
{ 
  double relaxation = 1.0;

  if (error_control_ & FLOW_TI_ERROR_CONTROL_SATURATION) {
    Epetra_MultiVector dSdP(uold);
    relperm_->Compute_dSdP(uold, dSdP);

    for (int c = 0; c < ncells_owned; c++) {
      double diff = dSdP[0][c] * fabs(unew[0][c] - uold[0][c]);
      if (diff > 3e-2) relaxation = std::min(relaxation, 3e-2 / diff);
    }
  }

  if (error_control_ & FLOW_TI_ERROR_CONTROL_PRESSURE) {
    for (int c = 0; c < ncells_owned; c++) {
      double diff = fabs(unew[0][c] - uold[0][c]);
      double umax = std::max(fabs(unew[0][c]), fabs(uold[0][c]));
      if (diff > 1e-2 * umax) relaxation = std::min(relaxation, 1e-2 * umax / diff);
    }
  }

#ifdef HAVE_MPI
  double relaxation_tmp = relaxation;
  mesh_->get_comm()->MinAll(&relaxation_tmp, &relaxation, 1);  // find the global minimum
#endif

  return relaxation;
}

}  // namespace Flow
}  // namespace Amanzi

