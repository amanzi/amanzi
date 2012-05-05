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

#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
*  Wrapper for advance to steady-state routines.                                                    
****************************************************************** */
int Richards_PK::AdvanceToSteadyState()
{
  flow_status_++;  // indicates intermediate state

  // initialize pressure and saturation at T=0.
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& water_saturation = FS->ref_water_saturation();

  // start iterations
  int ierr = 0;
  if (ti_method_sss == FLOW_TIME_INTEGRATION_PICARD) {
    ierr = AdvanceSteadyState_Picard();
  } else if (ti_method_sss == FLOW_TIME_INTEGRATION_BACKWARD_EULER) {
    ierr = AdvanceSteadyState_BackwardEuler();
  } else if (ti_method_sss == FLOW_TIME_INTEGRATION_BDF1) {
    ierr = AdvanceSteadyState_BDF1();
  } else if (ti_method_sss == FLOW_TIME_INTEGRATION_BDF2) {
    ierr = AdvanceSteadyState_BDF2();
  }

  if (ierr == 0) flow_status_ = FLOW_STATUS_STEADY_STATE_COMPLETE;
  return ierr;
}


/* ******************************************************************* 
* Performs one time step of size dT using first-order time integrator.
******************************************************************* */
int Richards_PK::AdvanceSteadyState_BDF1()
{
  T_internal = T0_sss;
  dT = dT0_sss;
  bool last_step = false;

  int itrs = 0;
  while (itrs < max_itrs_sss && T_internal < T1_sss) {
    if (itrs == 0) {  // initialization of BDF1
      Epetra_Vector udot(*super_map_);
      ComputeUDot(T0_sss, *solution, udot);
      bdf1_dae->set_initial_state(T0_sss, *solution, udot);

      int ierr;
      update_precon(T0_sss, *solution, dT0_sss, ierr);
    }

    double dTnext;
    bdf1_dae->bdf1_step(dT, *solution, dTnext);
    bdf1_dae->commit_solution(dT, *solution);
    bdf1_dae->write_bdf1_stepping_statistics();

    T_internal = bdf1_dae->most_recent_time();
    dT = dTnext;
    itrs++;

    double Tdiff = T1_sss - T_internal;
    if (dTnext > Tdiff) {
      dT = Tdiff * 0.99999991;  // To avoid hitting the wrong BC
      last_step = true;
    }
    if (last_step && dT < 1e-3) break;
  }

  num_nonlinear_steps = itrs;
  return 0;
}


/* ******************************************************************* 
* Performs one time step of size dT using second-order time integrator.
******************************************************************* */
int Richards_PK::AdvanceSteadyState_BDF2()
{
  T_internal = T0_sss;
  dT = dT0_sss;
  bool last_step = false;

  int itrs = 0;
  while (itrs < max_itrs_sss && T_internal < T1_sss) {
    if (itrs == 0) {  // initialization of BDF2
      Epetra_Vector udot(*super_map_);
      ComputeUDot(T0_sss, *solution, udot);
      bdf2_dae->set_initial_state(T0_sss, *solution, udot);

      int ierr;
      update_precon(T0_sss, *solution, dT0_sss, ierr);
    }

    double dTnext;
    bdf2_dae->bdf2_step(dT, 0.0, *solution, dTnext);
    bdf2_dae->commit_solution(dT, *solution);
    bdf2_dae->write_bdf2_stepping_statistics();

    T_internal = bdf2_dae->most_recent_time();
    dT = dTnext;
    itrs++;

    double Tdiff = T1_sss - T_internal;
    if (dTnext > Tdiff) {
      dT = Tdiff * 0.99999991;  // To avoid hitting the wrong BC
      last_step = true;
    }
    if (last_step && dT < 1e-3) break;
  }

  num_nonlinear_steps = itrs;
  return 0;
}


/* ******************************************************************
* Calculates steady-state solution assuming that absolute and
* relative permeabilities do not depend explicitly on time.
* This is the experimental method.                                                 
****************************************************************** */
int Richards_PK::AdvanceSteadyState_Picard()
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;
  Epetra_Vector  residual(*solution);

  if (!is_matrix_symmetric) solver->SetAztecOption(AZ_solver, AZ_gmres);
  solver->SetAztecOption(AZ_output, AZ_none);
  solver->SetAztecOption(AZ_conv, AZ_rhs);

  int itrs = 0;
  double L2norm, L2error = 1.0;

  while (L2error > convergence_tol_sss && itrs < max_itrs_sss) {
    SetAbsolutePermeabilityTensor(K);

    if (!is_matrix_symmetric) {  // Define K and Krel_faces
      CalculateRelativePermeabilityFace(*solution_cells);
      for (int c = 0; c < K.size(); c++) K[c] *= rho / mu;
    } else {  // Define K and Krel_cells, Krel_faces is always one
      CalculateRelativePermeabilityCell(*solution_cells);
      for (int c = 0; c < K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;
    }

    // create algebraic problem
    matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    rhs = matrix->get_rhs();  // export RHS from the matrix class

    // create preconditioner
    int disc_method = AmanziFlow::FLOW_MFD3D_TWO_POINT_FLUX;
    preconditioner->createMFDstiffnessMatrices(disc_method, K, *Krel_faces);
    preconditioner->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, preconditioner);
    preconditioner->applyBoundaryConditions(bc_markers, bc_values);
    preconditioner->assembleGlobalMatrices();
    preconditioner->computeSchurComplement(bc_markers, bc_values);
    preconditioner->update_ML_preconditioner();

    // check convergence of non-linear residual
    L2error = matrix->computeResidual(solution_new, residual);
    residual.Norm2(&L2error);
    rhs->Norm2(&L2norm);
    L2error /= L2norm;

    // call AztecOO solver
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess

    solver->Iterate(max_itrs, convergence_tol);
    int num_itrs = solver->NumIters();
    double linear_residual = solver->TrueResidual();

    // update relaxation
    double relaxation = 0.2;
    for (int c = 0; c < ncells_owned; c++) {
      double diff = fabs(solution_new[c] - solution_old[c]);
      double umax = std::max(fabs(solution_new[c]), fabs(solution_old[c]));
      if (diff > 5e-3 * umax) relaxation = std::min(relaxation, 5e-3 * umax / diff);
    }

    if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
      std::printf("Picard:%4d   Pressure(res=%9.4e, rhs=%9.4e, relax=%8.3e)  solver(%8.3e, %4d)\n",
          itrs, L2error, L2norm, relaxation, linear_residual, num_itrs);
    }

    int ndof = ncells_owned + nfaces_owned;
    for (int c = 0; c < ndof; c++) {
      solution_new[c] = (1.0 - relaxation) * solution_old[c] + relaxation * solution_new[c];
      solution_old[c] = solution_new[c];
    }

    T_internal += dT;
    itrs++;
  }

  num_nonlinear_steps = itrs;
  return 0;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

