/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MechanicsElasticity PK

*/

#include "WhetStoneDefs.hh"

#include "MechanicsElasticity_PK.hh"

namespace Amanzi {
namespace Mechanics {

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void
MechanicsElasticity_PK::FunctionalResidual(double t_old,
                                           double t_new,
                                           Teuchos::RCP<TreeVector> u_old,
                                           Teuchos::RCP<TreeVector> u_new,
                                           Teuchos::RCP<TreeVector> f)
{
  double dtp = t_new - t_old;

  // refresh data
  UpdateSourceBoundaryData_(t_old, t_new);

  // assemble residual using linear operator
  op_matrix_elas_->global_operator()->Init();
  op_matrix_elas_->UpdateMatrices();
  op_matrix_elas_->ApplyBCs(true, true, true);

  // Teuchos::RCP<CompositeVector> rhs = op_matrix_->rhs();
  // AddSourceTerms(*rhs);

  // compute negative residual, A u - f
  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *f->Data());
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.
****************************************************************** */
int
MechanicsElasticity_PK::ApplyPreconditioner(
  Teuchos::RCP<const TreeVector> X, Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  return op_preconditioner_->ApplyInverse(*X->Data(), *Y->Data());
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void
MechanicsElasticity_PK::UpdatePreconditioner(
  double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  double t_old = tp - dtp;
  // Teuchos::RCP<const CompositeVector> uu = u->SubVector(0)->Data();

  // refresh data
  UpdateSourceBoundaryData_(t_old, tp);

  // populate elastic operator
  auto global_op = op_preconditioner_elas_->global_operator();
  global_op->Init();
  op_preconditioner_elas_->UpdateMatrices();
  op_preconditioner_elas_->ApplyBCs(true, true, true);

  op_preconditioner_->ComputeInverse();
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.
* This is a wrapper for various error control methods.
****************************************************************** */
double
MechanicsElasticity_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  if (!u->Data()->HasComponent("node")) return 0.0;

  const Epetra_MultiVector& uv = *u->Data()->ViewComponent("node");
  const Epetra_MultiVector& duv = *du->Data()->ViewComponent("node");

  double error(0.0);

  for (int v = 0; v < nnodes_owned_; v++) {
    for (int k = 0; k < dim; ++k) {
      double tmp = fabs(duv[k][v]) / (fabs(uv[k][v]) + 1.0);
      error = std::max(error, tmp);
    }
  }

  return error;
}


/* ******************************************************************
* Calculate sources and boundary conditions for operators.
****************************************************************** */
void
MechanicsElasticity_PK::UpdateSourceBoundaryData_(double t_old, double t_new)
{
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(t_old, t_new);
    // bcs_[i]->ComputeSubmodel(mesh_);
  }

  // populate global arrays with data
  ComputeOperatorBCs();
}


/* ******************************************************************
* Add a boundary marker to used faces.
* WARNING: we can skip update of ghost boundary faces, b/c they
* should be always owned.
****************************************************************** */
void
MechanicsElasticity_PK::ComputeOperatorBCs()
{
  int d(mesh_->getSpaceDimension()), mv(0), mf(1);

  dirichlet_bc_ = 0;

  for (int i = 0; i < op_bcs_.size(); ++i) {
    std::vector<int>& bc_model = op_bcs_[i]->bc_model();
    for (int n = 0; n < bc_model.size(); n++) { bc_model[n] = Operators::OPERATOR_BC_NONE; }

    if (op_bcs_[i]->type() == WhetStone::DOF_Type::POINT) {
      mv = i;
      mf = 1 - mv;
      std::vector<AmanziGeometry::Point>& bc_value = op_bcs_[i]->bc_value_point();
      for (int n = 0; n < bc_value.size(); n++) { bc_value[n] = AmanziGeometry::Point(d); }
    }
  }

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "no slip" && bcs_[i]->type() == WhetStone::DOF_Type::POINT) {
      std::vector<int>& bc_model = op_bcs_[mv]->bc_model();
      std::vector<AmanziGeometry::Point>& bc_value = op_bcs_[mv]->bc_value_point();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;
        bc_model[n] = Operators::OPERATOR_BC_DIRICHLET;
        for (int k = 0; k < d; ++k) { bc_value[n][k] = it->second[k]; }
        dirichlet_bc_++;
      }
    }

    if (bcs_[i]->get_bc_name() == "no slip" &&
        bcs_[i]->type() == WhetStone::DOF_Type::NORMAL_COMPONENT) {
      std::vector<int>& bc_model = op_bcs_[mf]->bc_model();
      std::vector<double>& bc_value = op_bcs_[mf]->bc_value();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;
        bc_model[n] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[n] = it->second[0];
      }
    }
  }
}

} // namespace Mechanics
} // namespace Amanzi
