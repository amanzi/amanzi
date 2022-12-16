/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include "InverseFactory.hh"
#include "Mesh_Algorithms.hh"
#include "OperatorDefs.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "RemapUtils.hh"

#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Solve single phase problem using boundary conditions at time t_old.
* We populate both matrix and preconditoner here but use only the
* preconditioner. Matrix may be used in external flux calculation.
* Moving flux calculation here impose restrictions on multiple
* possible scenarios of data flow.
****************************************************************** */
void
Richards_PK::SolveFullySaturatedProblem(double t_old,
                                        CompositeVector& u,
                                        const std::string& solver_name)
{
  UpdateSourceBoundaryData(t_old, t_old, u);

  const auto& mu = S_->GetPtr<CompositeVector>(viscosity_liquid_key_, Tags::DEFAULT);
  alpha_upwind_->PutScalarMasterAndGhosted(molar_rho_);
  Operators::CellToFace_ScaleInverse(mu, alpha_upwind_);

  alpha_upwind_dP_->PutScalarMasterAndGhosted(0.0);

  // create diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true, true);

  // create diffusion preconditioner
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_preconditioner_diff_->ApplyBCs(true, true, true);

  Teuchos::ParameterList plist = linear_operator_list_->sublist(solver_name);
  AmanziSolvers::setMakeOneIterationCriteria(plist);
  auto solver = AmanziSolvers::createIterativeMethod(plist, op_matrix_, op_preconditioner_);
  solver->InitializeInverse();
  solver->ComputeInverse();

  CompositeVector& rhs = *op_matrix_->rhs();
  int ierr = solver->ApplyInverse(rhs, u);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    int num_itrs = solver->num_itrs();
    int code = solver->returned_code();
    double pnorm;
    u.Norm2(&pnorm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "saturated solver (" << solver->name() << "): ||p,lambda||=" << pnorm
               << " itr=" << num_itrs << " code=" << code << std::endl;

    CompositeVector r(rhs);
    op_matrix_->ComputeResidual(*solution, r);
    double residual;
    r.Norm2(&residual);
    *vo_->os() << "true l2 residual: ||r||=" << residual << std::endl;
  }

  // catastrophic failure
  if (ierr < 0) {
    Errors::Message msg;
    msg << "Richards_LinearOperator error: " << solver->returned_code_string();
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Enforce constraints, using new BCs at time t_new, by solving the
* block-diagonalized problem. Algorithm is based on de-coupling
* pressure-lambda system via zeroing-out off-diagonal blocks.
****************************************************************** */
void
Richards_PK::EnforceConstraints(double t_new, Teuchos::RCP<CompositeVector> u)
{
  std::vector<int>& bc_model = op_bc_->bc_model();

  UpdateSourceBoundaryData(t_new, t_new, *u);

  CompositeVector utmp(*u);
  Epetra_MultiVector& utmp_face = *utmp.ViewComponent("face");
  Epetra_MultiVector& u_face = *u->ViewComponent("face");

  // update diffusion coefficient
  // -- function
  vol_flowrate_copy->ScatterMasterToGhosted("face");

  pressure_eval_->SetChanged();
  auto& alpha = S_->GetW<CompositeVector>(alpha_key_, alpha_key_);
  S_->GetEvaluator(alpha_key_).Update(*S_, "flow");

  *alpha_upwind_->ViewComponent("cell") = *alpha.ViewComponent("cell");
  Operators::BoundaryFacesToFaces(bc_model, alpha, *alpha_upwind_);
  upwind_->Compute(*vol_flowrate_copy, bc_model, *alpha_upwind_);

  // -- derivative
  S_->GetEvaluator(alpha_key_).UpdateDerivative(*S_, passwd_, pressure_key_, Tags::DEFAULT);
  auto& alpha_dP = S_->GetDerivativeW<CompositeVector>(
    alpha_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, alpha_key_);

  *alpha_upwind_dP_->ViewComponent("cell") = *alpha_dP.ViewComponent("cell");
  Operators::BoundaryFacesToFaces(bc_model, alpha_dP, *alpha_upwind_dP_);
  upwind_->Compute(*vol_flowrate_copy, bc_model, *alpha_upwind_dP_);

  // modify relative permeability coefficient for influx faces
  UpwindInflowBoundary(u);
  // UpwindInflowBoundary_New(u);

  // calculate diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true, true);
  op_matrix_diff_->ModifyMatrices(*u);

  // calculate diffusion preconditioner
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_preconditioner_diff_->ApplyBCs(true, true, true);
  op_preconditioner_diff_->ModifyMatrices(*u);

  // solve non-symmetric problem
  Teuchos::ParameterList lin_op_list = linear_operator_list_->sublist(solver_name_constraint_);
  auto solver = AmanziSolvers::createIterativeMethod(lin_op_list, op_matrix_, op_preconditioner_);
  solver->InitializeInverse();
  solver->ComputeInverse();

  CompositeVector& rhs = *op_preconditioner_->rhs();
  int ierr = solver->ApplyInverse(rhs, utmp);

  u_face = utmp_face;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    int code = solver->returned_code();
    double pnorm;
    u->Norm2(&pnorm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "constraints solver (" << solver->name() << "): ||p,lambda||=" << pnorm
               << " itr=" << num_itrs << " code=" << code << std::endl;
  }

  // catastrophic failure
  if (ierr < 0) {
    Errors::Message msg;
    msg << "Richards::EnforceConstraints error: " << solver->returned_code_string();
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Calculates rel perm on the upwind boundary using a FV model.
****************************************************************** */
void
Richards_PK::UpwindInflowBoundary(Teuchos::RCP<const CompositeVector> u)
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();

  const auto& mu_c = *S_->Get<CompositeVector>(viscosity_liquid_key_).ViewComponent("cell");
  const auto& u_cell = *u->ViewComponent("cell");
  auto& k_face = *alpha_upwind_->ViewComponent("face", true);

  for (int f = 0; f < nfaces_owned; f++) {
    if ((bc_model[f] == Operators::OPERATOR_BC_NEUMANN ||
         bc_model[f] == Operators::OPERATOR_BC_MIXED) &&
        bc_value[f] < 0.0) {
      int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);

      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);
      double Knn = ((K[c] * normal) * normal) / (area * area);
      // old version
      // double save = 3.0;
      // k_face[0][f] = std::min(1.0, -save * bc_value[f] * mu_c[0][c] / (Knn * rho_ * rho_ * g_));
      // k_face[0][f] *= rho_ / mu_c[0][c];
      double value = bc_value[f] / flux_units_;
      double kr1 = wrm_->second[(*wrm_->first)[c]]->k_relative(atm_pressure_ - u_cell[0][c]);
      double kr2 = std::min(1.0, -value * mu_c[0][c] / (Knn * rho_ * rho_ * g_));
      k_face[0][f] = (molar_rho_ / mu_c[0][c]) * (kr1 + kr2) / 2;
    }
  }

  alpha_upwind_->ScatterMasterToGhosted("face");
}


/* ******************************************************************
* Calculates rel perm on the upwind boundary using a FV model.
****************************************************************** */
void
Richards_PK::UpwindInflowBoundary_New(Teuchos::RCP<const CompositeVector> u)
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();

  const auto& mu_c = *S_->Get<CompositeVector>(viscosity_liquid_key_).ViewComponent("cell");
  auto& k_face = *alpha_upwind_->ViewComponent("face", true);

  for (int f = 0; f < nfaces_owned; f++) {
    if ((bc_model[f] == Operators::OPERATOR_BC_NEUMANN ||
         bc_model[f] == Operators::OPERATOR_BC_MIXED) &&
        bc_value[f] < 0.0) {
      int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);

      double value = DeriveBoundaryFaceValue(f, *u, wrm_->second[(*wrm_->first)[c]]);
      double kr = wrm_->second[(*wrm_->first)[c]]->k_relative(atm_pressure_ - value);
      k_face[0][f] = kr * (molar_rho_ / mu_c[0][c]);
      // (*u->ViewComponent("face"))[0][f] = value;
    }
  }

  alpha_upwind_->ScatterMasterToGhosted("face");
}

} // namespace Flow
} // namespace Amanzi
