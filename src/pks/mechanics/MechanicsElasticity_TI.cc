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
#include "PorosityEvaluator.hh"

#include "MechanicsElasticity_PK.hh"

namespace Amanzi {
namespace Mechanics {

using CV_t = CompositeVector;

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
  UpdateSourceBoundaryData_(t_old, t_new);

  op_matrix_elas_->global_operator()->Init();

  // Add external forces
  auto rhs = op_matrix_->rhs();
  if (use_gravity_) AddGravityTerm_(*rhs);
  if (biot_model_) AddPressureGradient_(*rhs);

  op_matrix_elas_->UpdateMatrices();
  op_matrix_elas_->ApplyBCs(true, true, true);

  // compute negative residual, A u - f
  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *f->Data());
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.
****************************************************************** */
int
MechanicsElasticity_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X,
                                            Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  // return op_preconditioner_->ApplyInverse(*X->Data(), *Y->Data());
  return op_pc_solver_->ApplyInverse(*X->Data(), *Y->Data());
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void
MechanicsElasticity_PK::UpdatePreconditioner(double tp,
                                             Teuchos::RCP<const TreeVector> u,
                                             double dtp)
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

  // -- force component-wise coarsening
  /*
  auto block_indices = Teuchos::rcp(new std::vector<int>(nnodes_owned_ * dim_));
  for (int n = 0; n < nnodes_owned_ * dim_; ++n) (*block_indices)[n] = n % dim_;
  auto block_ids = std::make_pair(dim_, block_indices);
  op_preconditioner_elas_->global_operator()->set_coloring(dim_, block_indices);
  */

  op_preconditioner_->ComputeInverse();
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.
* This is a wrapper for various error control methods.
****************************************************************** */
double
MechanicsElasticity_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<const TreeVector> du)
{
  if (!u->Data()->HasComponent("node")) return 0.0;

  const auto& uv = *u->Data()->ViewComponent("node");
  const auto& duv = *du->Data()->ViewComponent("node");

  double error(0.0);

  for (int v = 0; v < nnodes_owned_; v++) {
    for (int k = 0; k < dim_; ++k) {
      double tmp = fabs(duv[k][v]) / (fabs(uv[k][v]) + 1.0);
      error = std::max(error, tmp);
    }
  }

#ifdef HAVE_MPI
  double buf = error;
  du->Comm()->MaxAll(&buf, &error, 1);
#endif

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
* We have four BCs on two locations (Node and Face) times two types
* (Point or Scalar).
* NOTE: The order of BCs is fixed and important: N+P, N+S, F+P, F+S
****************************************************************** */
void
MechanicsElasticity_PK::ComputeOperatorBCs()
{
  int d(mesh_->getSpaceDimension());
  dirichlet_bc_ = 0;

  // two nodal BCs are the first on the list
  for (int i = 0; i < op_bcs_.size(); ++i) {
    std::vector<int>& bc_model = op_bcs_[i]->bc_model();
    for (int n = 0; n < bc_model.size(); n++) { bc_model[n] = Operators::OPERATOR_BC_NONE; }
  }

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "no slip" && bcs_[i]->type() == WhetStone::DOF_Type::POINT &&
        bcs_[i]->kind() == AmanziMesh::Entity_kind::NODE) {
      std::vector<int>& bc_model = op_bcs_[0]->bc_model();
      std::vector<AmanziGeometry::Point>& bc_value = op_bcs_[0]->bc_value_point();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;
        bc_model[n] = Operators::OPERATOR_BC_DIRICHLET;
        for (int k = 0; k < d; ++k) { bc_value[n][k] = it->second[k]; }
        dirichlet_bc_++;
      }
    }

    if (bcs_[i]->get_bc_name() == "no slip" && bcs_[i]->type() == WhetStone::DOF_Type::POINT &&
        bcs_[i]->kind() == AmanziMesh::Entity_kind::FACE) {
      std::vector<int>& bc_model = op_bcs_[3]->bc_model();
      std::vector<double>& bc_value = op_bcs_[3]->bc_value();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;

        double dot(0.0);
        const auto& normal = mesh_->getFaceNormal(n);
        for (int k = 0; k < d; ++k) { dot += it->second[k] * normal[k]; }

        bc_model[n] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[n] = dot / mesh_->getFaceArea(n);
      }
    }

    if (bcs_[i]->get_bc_name() == "kinematic" &&
        bcs_[i]->type() == WhetStone::DOF_Type::NORMAL_COMPONENT) {
      std::vector<int>& bc_model = op_bcs_[3]->bc_model();
      std::vector<double>& bc_value = op_bcs_[3]->bc_value();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;
        bc_model[n] = Operators::OPERATOR_BC_KINEMATIC;
        bc_value[n] = it->second[0];
      }
    }
  }

  // facial BCs make the second group on the list
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "traction") {
      std::vector<int>& bc_model = op_bcs_[2]->bc_model();
      std::vector<AmanziGeometry::Point>& bc_value = op_bcs_[2]->bc_value_point();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;
        bc_model[n] = Operators::OPERATOR_BC_NORMAL_STRESS;
        for (int k = 0; k < d; ++k) { bc_value[n][k] = it->second[k]; }
      }
    }
  }
}


/* ******************************************************************
* Add forcing term due to gravity
****************************************************************** */
void
MechanicsElasticity_PK::AddGravityTerm_(CompositeVector& rhs)
{
  int d = mesh_->getSpaceDimension();
  double g = (S_->Get<AmanziGeometry::Point>("gravity"))[d - 1];

  const auto& rho_c = *S_->Get<CV_t>(particle_density_key_, Tags::DEFAULT).ViewComponent("cell");
  auto& rhs_v = *rhs.ViewComponent("node");

  rhs.PutScalarGhosted(0.0);

  for (int c = 0; c < ncells_owned_; ++c) {
    auto nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();

    double add = g * rho_c[0][c] * mesh_->getCellVolume(c) / nnodes;
    for (int n = 0; n < nnodes; ++n) rhs_v[d - 1][nodes[n]] += add;
  }

  rhs.GatherGhostedToMaster("node");
}


/* ******************************************************************
* Add pressure gradient
****************************************************************** */
void
MechanicsElasticity_PK::AddPressureGradient_(CompositeVector& rhs)
{
  int d = mesh_->getSpaceDimension();
  const auto& p = S_->Get<CV_t>("pressure", Tags::DEFAULT);

  auto eval = Teuchos::rcp_dynamic_cast<Flow::PorosityEvaluator>(
    S_->GetEvaluatorPtr("porosity", Tags::DEFAULT));

  if (p.HasComponent("face")) {
    const auto& p_f = *p.ViewComponent("face", true);
    auto& rhs_v = *rhs.ViewComponent("node");

    p.ScatterMasterToGhosted("face");
    rhs.PutScalarGhosted(0.0);

    for (int c = 0; c < ncells_owned_; ++c) {
      // pressure gradient
      double vol = mesh_->getCellVolume(c);
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
      const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
      int nfaces = faces.size();

      AmanziGeometry::Point grad(d);
      for (int i = 0; i < nfaces; ++i) {
        int f = faces[i];
        const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
        const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
        grad += normal * (p_f[0][f] * dirs[i]);
      }
      grad /= vol;

      double b = eval->getBiotCoefficient(c);

      // pressure gradient
      auto nodes = mesh_->getCellNodes(c);
      int nnodes = nodes.size();

      for (int i = 0; i < d; ++i) {
        double add = b * grad[i] * vol / nnodes;
        for (int n = 0; n < nnodes; ++n) rhs_v[i][nodes[n]] -= add;
      }
    }
  }

  rhs.GatherGhostedToMaster("node");
}

} // namespace Mechanics
} // namespace Amanzi
