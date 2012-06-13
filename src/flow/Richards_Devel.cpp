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
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.
* WARNING: temporary replacement for missing BDF1 time integrator.                                                    
****************************************************************** */
int Richards_PK::AdvanceToSteadyState_BackwardEuler()
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;

  Teuchos::RCP<Epetra_Vector> solution_old_cells = Teuchos::rcp(FS->CreateCellView(solution_old));
  Teuchos::RCP<Epetra_Vector> solution_new_cells = Teuchos::rcp(FS->CreateCellView(solution_new));

  if (! is_matrix_symmetric) solver->SetAztecOption(AZ_solver, AZ_gmres);
  solver->SetAztecOption(AZ_output, AZ_none);

  int max_itrs = ti_specs_sss_.max_itrs;
  double T1 = ti_specs_sss_.T1;

  int itrs = 0, ifail = 0;
  double L2error = 1.0;
  while (L2error > residual_tol_sss && itrs < max_itrs) {
    if (!is_matrix_symmetric) {  // Define K and Krel_faces
      CalculateRelativePermeabilityFace(*solution_cells);
      Krel_cells->PutScalar(1.0);
    } else {  // Define K and Krel_cells, Krel_faces is always one
      CalculateRelativePermeabilityCell(*solution_cells);
      Krel_faces->PutScalar(1.0);
    }

    // update boundary conditions
    double time = T_physics + dT;
    bc_pressure->Compute(time);
    bc_flux->Compute(time);
    bc_head->Compute(time);
    bc_seepage->Compute(time);
    ProcessBoundaryConditions(
        bc_pressure, bc_head, bc_flux, bc_seepage,
        *solution_faces, atm_pressure,
        bc_markers, bc_values);

    // create algebraic problem
    matrix_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
    matrix_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, matrix_);
    AddTimeDerivative_MFDfake(*solution_cells, dT, matrix_);
    matrix_->ApplyBoundaryConditions(bc_markers, bc_values);
    matrix_->AssembleGlobalMatrices();

    // create preconditioner
    preconditioner_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
    preconditioner_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, preconditioner_);
    AddTimeDerivative_MFDfake(*solution_cells, dT, preconditioner_);
    preconditioner_->ApplyBoundaryConditions(bc_markers, bc_values);
    preconditioner_->AssembleGlobalMatrices();
    preconditioner_->ComputeSchurComplement(bc_markers, bc_values);
    preconditioner_->UpdateML_Preconditioner();

    // call AztecOO solver
    rhs = matrix_->rhs();
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess

    solver->Iterate(max_itrs, convergence_tol);
    int num_itrs = solver->NumIters();
    double residual = solver->TrueResidual();

    // error estimates
    double sol_norm = FS->normLpCell(solution_new, 2.0);
    Epetra_Vector solution_diff(*solution_old_cells);
    solution_diff.Update(-1.0, *solution_new_cells, 1.0);
    L2error = ErrorNormPicardExperimental(*solution_old_cells, solution_diff);

    if (L2error > 1.0 && itrs && ifail < 5) {  // itrs=0 allows to avoid bad initial guess.
      dT /= 10;
      solution_new = solution_old;
      if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
        std::printf("Fail:%4d  Pressure(diff=%9.4e, sol=%9.4e)  solver(%8.3e,%3d), T=%9.3e dT=%7.2e\n",
            itrs, L2error, sol_norm, residual, num_itrs, T_physics, dT);
      }
      ifail++;
    } else {
      T_physics += dT;
      solution_old = solution_new;

      if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
        std::printf("Step:%4d  Pressure(diff=%9.4e, sol=%9.4e)  solver(%8.3e,%3d), T=%9.3e dT=%7.2e\n",
            itrs, L2error, sol_norm, residual, num_itrs, T_physics, dT);
      }

      ifail = 0;
      dT = std::min(dT*1.25, dTmax_sss);
      itrs++;
    }

    if (T_physics > T1) break;
  }

  ti_specs_sss_.num_itrs = itrs;
  return 0;
}


/* ******************************************************************
* Adds time derivative to cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::AddTimeDerivative_MFDfake(
    Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix_operator)
{
  std::vector<double>& Acc_cells = matrix_operator->Acc_cells();
  std::vector<double>& Fc_cells = matrix_operator->Fc_cells();

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume / dT_prec;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}


/* ******************************************************************
* Saturation should be in exact balance with Darcy fluxes in order to
* have extrema dimishing property for concentrations. 
* WARNING: is NOT used
****************************************************************** */
void Richards_PK::CalculateConsistentSaturation(const Epetra_Vector& flux, 
                                                const Epetra_Vector& ws_prev, Epetra_Vector& ws)
{
  // create a disctributed flux vector
  Epetra_Vector flux_d(mesh_->face_map(true));
  for (int f = 0; f < nfaces_owned; f++) flux_d[f] = flux[f];
  FS->CombineGhostFace2MasterFace(flux_d);

  const Epetra_Vector& phi = FS->ref_porosity();
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    ws[c] = ws_prev[c];
    double factor = dT / (phi[c] * mesh_->cell_volume(c));
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      ws[c] -= factor * flux_d[f] * dirs[n]; 
    }
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi

