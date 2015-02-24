/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "LinearOperatorFactory.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusion.hh"
#include "OperatorDiffusionFactory.hh"

#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Solve single phase problem using boundary conditions at time T0.
* We populate both matrix and preconditoner here but use only the
* preconditioner. Matrix may be used in external flux calculation. 
* Moving flux calculation here impose restrictions on multiple 
* possible scenarios of data flow.
****************************************************************** */
void Richards_PK::SolveFullySaturatedProblem(
    double T0, CompositeVector& u, const std::string& solver_name)
{
  UpdateSourceBoundaryData(T0, T0, u);
  rel_perm_->Krel()->PutScalar(1.0);
  rel_perm_->dKdP()->PutScalar(0.0);

  // create diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution);
  op_matrix_diff_->ApplyBCs(op_bc_);

  // create diffusion preconditioner
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, solution);
  op_preconditioner_diff_->ApplyBCs(op_bc_);
  op_preconditioner_->AssembleMatrix();
  op_preconditioner_->InitPreconditioner(ti_specs->preconditioner_name, *preconditioner_list_);

  AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> sfactory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace> >
     solver = sfactory.Create(solver_name, *linear_operator_list_, op_matrix_, op_preconditioner_);

  solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);  // Make at least one iteration

  CompositeVector& rhs = *op_matrix_->rhs();
  int ierr = solver->ApplyInverse(rhs, u);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    int code = solver->returned_code();
    double pnorm;
    u.Norm2(&pnorm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "saturated solver (" << solver->name() 
               << "): ||p,lambda||=" << pnorm << " itr=" << num_itrs 
               << " code=" << code << std::endl;
  }
  if (ierr != 0) {
    Errors::Message msg;
    msg << "\nLinear solver returned an unrecoverable error code.\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Enforce constraints at time Tp by solving diagonalized MFD problem.
* Algorithm is based on de-coupling pressure-lambda system.
****************************************************************** */
void Richards_PK::EnforceConstraints(double Tp, CompositeVector& u)
{
  UpdateSourceBoundaryData(Tp, Tp, u);

  CompositeVector utmp(u);
  Epetra_MultiVector& utmp_face = *utmp.ViewComponent("face");
  Epetra_MultiVector& u_face = *u.ViewComponent("face");
  Epetra_MultiVector& u_cell = *u.ViewComponent("cell");

  // update relative permeability coefficients
  darcy_flux_copy->ScatterMasterToGhosted("face");
  rel_perm_->Compute(u);

  RelativePermeabilityUpwindFn func1 = &RelativePermeability::Value;
  upwind_->Compute(*darcy_flux_upwind, u, bc_model, bc_value,
                   *rel_perm_->Krel(), *rel_perm_->Krel(), func1);

  RelativePermeabilityUpwindFn func2 = &RelativePermeability::Derivative;
  upwind_->Compute(*darcy_flux_upwind, u, bc_model, bc_value,
                   *rel_perm_->dKdP(), *rel_perm_->dKdP(), func2);

  // modify relative permeability coefficient for influx faces
  if (ti_specs->inflow_krel_correction) {
    Epetra_MultiVector& k_face = *rel_perm_->Krel()->ViewComponent("face", true);
    AmanziMesh::Entity_ID_List cells;

    for (int f = 0; f < nfaces_wghost; f++) {
      if ((bc_model[f] == Operators::OPERATOR_BC_NEUMANN || 
           bc_model[f] == Operators::OPERATOR_BC_MIXED) && bc_value[f] < 0.0) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        int c = cells[0];

        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        double area = mesh_->face_area(f);
        double Knn = ((K[c] * normal) * normal) / (area * area);
        // double save = 3.0;
        // k_face[0][f] = std::min(1.0, -save * bc_value[f] * mu_ / (Knn * rho_ * rho_ * g_));
        double kr1 = rel_perm_->Value(c, u_cell[0][c]);
        double kr2 = std::min(1.0, -bc_value[f] * mu_ / (Knn * rho_ * rho_ * g_));
        k_face[0][f] = (kr1 + kr2) / 2;
      } 
    }
  }

  // calculate diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution);
  op_matrix_diff_->ApplyBCs(op_bc_);
  op_matrix_diff_->ModifyMatrices(u);

  // calculate diffusion preconditioner
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, solution);
  op_preconditioner_diff_->ApplyBCs(op_bc_);
  op_preconditioner_diff_->ModifyMatrices(u);
  op_preconditioner_->AssembleMatrix();
  op_preconditioner_->InitPreconditioner(ti_specs->preconditioner_name, *preconditioner_list_);

  // solve non-symmetric problem
  AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(ti_specs->solver_name_constraint, *linear_operator_list_, op_matrix_, op_preconditioner_);

  CompositeVector& rhs = *op_preconditioner_->rhs();
  int ierr = solver->ApplyInverse(rhs, utmp);

  u_face = utmp_face;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();
    int code = solver->returned_code();
    double pnorm;
    u.Norm2(&pnorm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "constraints solver (" << solver->name() 
               << "): ||p,lambda||=" << pnorm << " itr=" << num_itrs 
               << " code=" << code << std::endl;
  }
  if (ierr != 0) {
    Errors::Message msg;
    msg << "\nLinear solver returned an unrecoverable error code.\n";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Flow
}  // namespace Amanzi

