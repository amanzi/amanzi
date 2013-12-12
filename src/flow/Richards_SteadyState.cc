/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author:  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <vector>

#include "LinearOperatorFactory.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
*  Wrapper for advance to steady-state routines.                                                    
****************************************************************** */
int Richards_PK::AdvanceToSteadyState(double T0, double dT0)
{
  // update the axiliary flow state
  S_->GetFieldData("darcy_flux", passwd_)->ScatterMasterToGhosted("face");

  // override internal parameters
  T_physics = ti_specs_sss_.T0 = T0;
  dT = ti_specs_sss_.dT0 = dT0;

  int ierr = 0;
  int ti_method = ti_specs_sss_.ti_method;
  if (ti_method == FLOW_TIME_INTEGRATION_PICARD) {
    ierr = AdvanceToSteadyState_Picard(ti_specs_sss_);
  } else if (ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    ierr = AdvanceToSteadyState_BDF1(ti_specs_sss_);
  }

  Epetra_MultiVector& p = *solution->ViewComponent("cell");
  Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& ws_prev = *S_->GetFieldData("prev_water_saturation", passwd_)->ViewComponent("cell");

  DeriveSaturationFromPressure(p, ws);
  ws_prev = ws;

  return ierr;
}


/* ******************************************************************* 
* Performs one time step of size dT using first-order time integrator.
******************************************************************* */
int Richards_PK::AdvanceToSteadyState_BDF1(TI_Specs& ti_specs)
{
  bool last_step = false;

  int max_itrs = ti_specs.max_itrs;
  double T0 = ti_specs.T0;
  double T1 = ti_specs.T1;
  double dT0 = ti_specs.dT0;

  int itrs = 0;
  while (itrs < max_itrs && T_physics < T1) {
    if (itrs == 0) {  // initialization of BDF1
      Teuchos::RCP<CompositeVector> udot = Teuchos::rcp(new CompositeVector(*solution));
      // I do not know how to calculate du/dt in a robust way.
      // ComputeUDot(T0, solution, udot);
      udot->PutScalar(0.0);
      bdf1_dae->set_initial_state(T0, solution, udot);

      int ierr;
      update_precon(T0, solution, dT0);
    }

    double dTnext;
    try { 
      bdf1_dae->time_step(dT, dTnext, solution);
    } catch (int i) {
      dT /= 2;
      continue;
    }
    bdf1_dae->commit_solution(dT, solution);

    T_physics = bdf1_dae->time();
    dT = dTnext;
    itrs++;

    double Tdiff = T1 - T_physics;
    if (dTnext > Tdiff) {
      dT = Tdiff * 0.99999991;  // To avoid hitting the wrong BC
      last_step = true;
    }
    if (last_step && dT < 1e-3) break;
  }

  ti_specs.num_itrs = itrs;
  return 0;
}


/* ******************************************************************
* Calculates steady-state solution assuming that absolute and
* relative permeabilities do not depend explicitly on time.
* This is the experimental method.                                                 
****************************************************************** */
int Richards_PK::AdvanceToSteadyState_Picard(TI_Specs& ti_specs)
{
  // create verbosity object
  VerboseObject* vo = new VerboseObject("Amanzi::Picard", rp_list_); 

  CompositeVector  solution_old(*solution);
  CompositeVector& solution_new = *solution;
  CompositeVector  residual(*solution);

  Epetra_MultiVector& pold_cells = *solution_old.ViewComponent("cell");
  Epetra_MultiVector& pnew_cells = *solution_new.ViewComponent("cell");

  // update steady state boundary conditions
  double time = T_physics;
  bc_pressure->Compute(time);
  bc_flux->Compute(time);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(time);
  else
    bc_head->ComputeShift(time, shift_water_table_->Values());

  // update steady state source conditons
  if (src_sink != NULL) {
    if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink->ComputeDistribute(time, Kxy->Values()); 
    } else {
      src_sink->ComputeDistribute(time, NULL);
    }
  }

  int max_itrs_nonlinear = ti_specs.max_itrs;
  double residual_tol_nonlinear = ti_specs.residual_tol;
  LinearSolver_Specs& ls_specs = ti_specs.ls_specs;

  int itrs = 0;
  double L2norm, L2error = 1.0;

  while (L2error > residual_tol_nonlinear && itrs < max_itrs_nonlinear) {
    // update dynamic boundary conditions
    bc_seepage->Compute(time);
    ComputeBCs(*solution);

    // update permeabilities
    rel_perm->Compute(*solution, bc_model, bc_values);

    // create algebraic problem
    matrix_->CreateStiffnessMatricesRichards();
    matrix_->CreateRHSVectors();
    matrix_->AddGravityFluxesRichards(rho_, gravity_, K);
    matrix_->ApplyBoundaryConditions(bc_model, bc_values);
    matrix_->Assemble();

    Teuchos::RCP<CompositeVector> rhs = matrix_->rhs();  // export RHS from the matrix class
    if (src_sink != NULL) AddSourceTerms(*rhs);

    // create preconditioner
    preconditioner_->CreateStiffnessMatricesRichards();
    preconditioner_->CreateRHSVectors();
    preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner_->AssembleSchurComplement(bc_model, bc_values);
    preconditioner_->UpdatePreconditioner();

    // check convergence of non-linear residual
    L2error = matrix_->ComputeResidual(solution_new, residual);
    residual.Norm2(&L2error);
    rhs->Norm2(&L2norm);
    L2error /= L2norm;

    // solve linear problem
    AmanziSolvers::LinearOperatorFactory<FlowMatrix, CompositeVector, CompositeVectorSpace> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<FlowMatrix, CompositeVector, CompositeVectorSpace> >
       solver = factory.Create(ls_specs.solver_name, linear_operator_list_, matrix_, preconditioner_);

    solver->ApplyInverse(*rhs, *solution);

    int num_itrs_linear = solver->num_itrs();
    double linear_residual = solver->residual();

    // update relaxation
    double relaxation;
    relaxation = CalculateRelaxationFactor(pold_cells, pnew_cells);

    if (vo->getVerbLevel() >= Teuchos::VERB_HIGH) {
      Teuchos::OSTab tab = vo->getOSTab();
      *(vo->os()) << itrs << ": ||r||=" << L2error << " relax=" << relaxation 
                  << " lin_solver(" << linear_residual << ", " << num_itrs_linear << ")" << endl;
    }

    solution_new.Update(1.0 - relaxation, solution_old, relaxation);
    solution_old = solution_new;

    T_physics += dT;
    itrs++;
  }

  ti_specs.num_itrs = itrs;

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
    rel_perm->DerivedSdP(uold, dSdP);

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

}  // namespace AmanziFlow
}  // namespace Amanzi

