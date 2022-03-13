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
#include "InverseFactory.hh"
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

  Teuchos::RCP<const CompositeVector> mu = S_->GetFieldData(viscosity_liquid_key_);

  std::string linear_solver = plist.get<std::string>("linear solver");
  Teuchos::ParameterList lin_solve_list = linear_operator_list_->sublist(linear_solver);
  AmanziSolvers::setMakeOneIterationCriteria(lin_solve_list);
  auto solver = AmanziSolvers::createIterativeMethod(lin_solve_list, op_preconditioner_);
  solver->InitializeInverse();
  
  Teuchos::ParameterList& tmp_list = plist.sublist("picard parameters");
  int max_itrs_nonlinear = tmp_list.get<int>("maximum number of iterations", 400);
  double residual_tol_nonlinear = tmp_list.get<double>("convergence tolerance", 1e-8);

  int itrs = 0;
  double L2norm, L2error = 1.0;

  while (L2error > residual_tol_nonlinear && itrs < max_itrs_nonlinear) {
    // update dynamic boundary conditions
    for (int i = 0; i < bcs_.size(); i++) {
      if (bcs_[i]->get_bc_name() == "seepage") {
        bcs_[i]->Compute(time, time);
        bcs_[i]->ComputeSubmodel(mesh_);
      }
    }
    ComputeOperatorBCs(*solution);

    // update diffusion coefficients
    // -- function
    darcy_flux_copy->ScatterMasterToGhosted("face");

    pressure_eval_->SetFieldAsChanged(S_.ptr());
    auto alpha = S_->GetFieldData(alpha_key_, alpha_key_);
    S_->GetFieldEvaluator(alpha_key_)->HasFieldChanged(S_.ptr(), "flow");
  
    *alpha_upwind_->ViewComponent("cell") = *alpha->ViewComponent("cell");
    Operators::BoundaryFacesToFaces(bc_model, *alpha, *alpha_upwind_);
    upwind_->Compute(*darcy_flux_copy, *solution, bc_model, *alpha_upwind_);

    // -- derivative 
    Key der_key = Keys::getDerivKey(alpha_key_, pressure_key_);
    S_->GetFieldEvaluator(alpha_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, pressure_key_);
    auto alpha_dP = S_->GetFieldData(der_key);

    *alpha_upwind_dP_->ViewComponent("cell") = *alpha_dP->ViewComponent("cell");
    Operators::BoundaryFacesToFaces(bc_model, *alpha_dP, *alpha_upwind_dP_);
    upwind_->Compute(*darcy_flux_copy, *solution, bc_model, *alpha_upwind_dP_);

    // create algebraic problem (matrix = preconditioner)
    op_preconditioner_->Init();
    op_preconditioner_diff_->UpdateMatrices(darcy_flux_copy.ptr(), solution.ptr());
    op_preconditioner_diff_->UpdateMatricesNewtonCorrection(darcy_flux_copy.ptr(), Teuchos::null, molar_rho_);
    op_preconditioner_diff_->ApplyBCs(true, true, true);

    Teuchos::RCP<CompositeVector> rhs = op_preconditioner_->rhs();  // export RHS from the matrix class
    AddSourceTerms(*rhs);

    // update inverse
    solver->ComputeInverse();

    // check convergence of non-linear residual
    op_preconditioner_->ComputeResidual(solution_new, residual);
    residual.Norm2(&L2error);
    rhs->Norm2(&L2norm);
    L2error /= L2norm;

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
  double patm = *S_->GetScalarData("atmospheric_pressure");

  if (error_control_ & FLOW_TI_ERROR_CONTROL_SATURATION) {
    for (int c = 0; c < ncells_owned; c++) {
      double dSdP = wrm_->second[(*wrm_->first)[c]]->k_relative(patm - uold[0][c]);
      double diff = dSdP * fabs(unew[0][c] - uold[0][c]);
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

