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
// When this is used by Init, shouldn't this also get preconditioner_name_ini instead of preconditioner_name? --etc
void Richards_PK::SolveFullySaturatedProblem(
    double T0, CompositeVector& u, const std::string& solver_name)
{
  UpdateSourceBoundaryData(T0, T0, u);
  krel_->PutScalar(molar_rho_ / mu_);
  dKdP_->PutScalar(0.0);

  // create diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true);

  // create diffusion preconditioner
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_preconditioner_diff_->ApplyBCs(true);
  op_preconditioner_->AssembleMatrix();
  op_preconditioner_->InitPreconditioner(preconditioner_name_, *preconditioner_list_);

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
void Richards_PK::EnforceConstraints(double Tp, Teuchos::RCP<CompositeVector> u)
{
  UpdateSourceBoundaryData(Tp, Tp, *u);

  CompositeVector utmp(*u);
  Epetra_MultiVector& utmp_face = *utmp.ViewComponent("face");
  Epetra_MultiVector& u_face = *u->ViewComponent("face");
  Epetra_MultiVector& u_cell = *u->ViewComponent("cell");

  // update relative permeability coefficients and upwind it
  darcy_flux_copy->ScatterMasterToGhosted("face");

  relperm_->Compute(u, krel_);
  RelPermUpwindFn func1 = &RelPerm::Compute;
  upwind_->Compute(*darcy_flux_upwind, *u, bc_model, bc_value, *krel_, *krel_, func1);
  krel_->ScaleMasterAndGhosted(molar_rho_ / mu_);

  relperm_->ComputeDerivative(u, dKdP_);
  RelPermUpwindFn func2 = &RelPerm::ComputeDerivative;
  upwind_->Compute(*darcy_flux_upwind, *u, bc_model, bc_value, *dKdP_, *dKdP_, func2);
  dKdP_->ScaleMasterAndGhosted(molar_rho_ / mu_);

  // modify relative permeability coefficient for influx faces
  bool inflow_krel_correction(true);
  if (inflow_krel_correction) {
    Epetra_MultiVector& k_face = *krel_->ViewComponent("face", true);
    AmanziMesh::Entity_ID_List cells;

    for (int f = 0; f < nfaces_owned; f++) {
      if ((bc_model[f] == Operators::OPERATOR_BC_NEUMANN || 
           bc_model[f] == Operators::OPERATOR_BC_MIXED) && bc_value[f] < 0.0) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        int c = cells[0];

        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        double area = mesh_->face_area(f);
        double Knn = ((K[c] * normal) * normal) / (area * area);
        // double save = 3.0;
        // k_face[0][f] = std::min(1.0, -save * bc_value[f] * mu_ / (Knn * rho_ * rho_ * g_));
        double value = bc_value[f] / flux_units_;
        double kr1 = relperm_->Compute(c, u_cell[0][c]);
        double kr2 = std::min(1.0, -value * mu_ / (Knn * rho_ * rho_ * g_));
        k_face[0][f] = (molar_rho_ / mu_) * (kr1 + kr2) / 2;
      } 
    }

    krel_->ScatterMasterToGhosted("face");
  }

  // calculate diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true);
  op_matrix_diff_->ModifyMatrices(*u);

  // calculate diffusion preconditioner
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_preconditioner_diff_->ApplyBCs(true);
  op_preconditioner_diff_->ModifyMatrices(*u);
  op_preconditioner_->AssembleMatrix();
  op_preconditioner_->InitPreconditioner(preconditioner_name_, *preconditioner_list_);

  // solve non-symmetric problem
  AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace> >
      solver = factory.Create(solver_name_constraint_, *linear_operator_list_, op_matrix_, op_preconditioner_);

  CompositeVector& rhs = *op_preconditioner_->rhs();
  int ierr = solver->ApplyInverse(rhs, utmp);

  u_face = utmp_face;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();
    int code = solver->returned_code();
    double pnorm;
    u->Norm2(&pnorm);

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

