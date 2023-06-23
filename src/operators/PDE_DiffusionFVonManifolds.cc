/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// Amanzi
#include "Mesh_Algorithms.hh"


// Operators
#include "Op.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"
#include "Op_Face_Cell.hh"
#include "OperatorDefs.hh"
#include "Operator_Cell.hh"
#include "PDE_DiffusionFVonManifolds.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization
****************************************************************** */
void
PDE_DiffusionFVonManifolds::Init_(Teuchos::ParameterList& plist)
{
  local_op_schema_ = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL;

  if (global_op_ == Teuchos::null) {
    global_op_schema_ = OPERATOR_SCHEMA_DOFS_CELL;

    auto cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);
    cvs->AddComponent("cell", AmanziMesh::CELL, 1);

    global_op_ = Teuchos::rcp(new Operator_Cell(cvs, plist, global_op_schema_));
  } else {
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // create an local op and register it with the global operator
  std::string name = "Diffusion: FACE_CELL Manifolds";
  local_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
  global_op_->OpPushBack(local_op_);

  // upwind options
  Errors::Message msg;
  std::string uwname = plist.get<std::string>("nonlinear coefficient", "upwind: face");
  if (uwname == "none") {
    little_k_ = OPERATOR_LITTLE_K_NONE;
  } else if (uwname == "upwind: face") {
    little_k_ = OPERATOR_LITTLE_K_UPWIND;
  } else {
    msg << "DiffusionFVonManifolds: unknonw or unsupported upwind scheme.";
    Exceptions::amanzi_throw(msg);
  }

  // other options
  gravity_ = plist.get<bool>("gravity");

  // solution-independent data
  auto cvs = Operators::CreateManifoldCVS(mesh_);
  beta_ = Teuchos::rcp(new CompositeVector(*cvs, true));
}


/* ******************************************************************
* Setup method: tensor coefficient
****************************************************************** */
void
PDE_DiffusionFVonManifolds::SetTensorCoefficient(
  const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K)
{
  beta_initialized_ = false;
  K_ = K;
}


/* ******************************************************************
* Setup method: upwinded scalar coefficient
****************************************************************** */
void
PDE_DiffusionFVonManifolds::SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                                 const Teuchos::RCP<const CompositeVector>& dkdp)
{
  k_ = k;
  dkdp_ = dkdp;

  // we can use only one component, either cell of face
  if (k_ != Teuchos::null) {
    beta_initialized_ = false;
    AMANZI_ASSERT(k_->HasComponent("face"));
    AMANZI_ASSERT(k_->Map().Ghosted());

    if (dkdp_ != Teuchos::null) {
      AMANZI_ASSERT(dkdp_->HasComponent("face"));
      AMANZI_ASSERT(dkdp_->Map().Ghosted());
    }
  }
}


/* ******************************************************************
* Populate face-based nxn matrices, n=1 for boundary faces.
****************************************************************** */
void
PDE_DiffusionFVonManifolds::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                           const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!beta_initialized_) ComputeBeta_();

  const auto& beta_f = *beta_->ViewComponent("face", true);
  const auto& fmap = *beta_->Map().Map("face", true);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();

  int d = mesh_->space_dimension();
  auto& rhs_c = *global_op_->rhs()->ViewComponent("cell");
  global_op_->rhs()->PutScalarGhosted(0.0);

  // preparing upwind data
  Teuchos::RCP<const Epetra_MultiVector> k_f = Teuchos::null;
  if (k_ != Teuchos::null) {
    k_->ScatterMasterToGhosted("face");
    k_f = k_->ViewComponent("face", true);
  }

  // updating matrix blocks
  AmanziMesh::Entity_ID_List cells;
  WhetStone::DenseVector v(2), av(2);

  for (int f = 0; f != nfaces_owned; ++f) {
    int g = fmap.FirstPointInElement(f);
    int ndofs = fmap.ElementSize(f);

    WhetStone::DenseMatrix Aface(ndofs, ndofs);
    Aface = 0.0;

    double ti, tj, sum(0.0);
    for (int i = 0; i != ndofs; ++i) {
      ti = beta_f[0][g + i] * (k_f.get() ? (*k_f)[0][g + i] : 1.0);
      sum += ti;
    }

    if (ndofs == 1) {
      Aface(0, 0) = sum;
    } else {
      if (sum > 0.0) sum = 1.0 / sum;

      for (int i = 0; i != ndofs; ++i) {
        ti = beta_f[0][g + i] * (k_f.get() ? (*k_f)[0][g + i] : 1.0);
        for (int j = i + 1; j != ndofs; ++j) {
          tj = beta_f[0][g + j] * (k_f.get() ? (*k_f)[0][g + j] : 1.0);
          Aface(i, i) += ti * ti * sum;
          Aface(j, j) += tj * tj * sum;
          Aface(i, j) = -ti * tj * sum;
          Aface(j, i) = -ti * tj * sum;
        }
      }
    }
    local_op_->matrices[f] = Aface;

    if (gravity_) {
      if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        // skip
      } else if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        double factor = rho_ * norm(g_);
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int c = cells[0];

        double zc = (mesh_->cell_centroid(c))[d - 1];
        double zf = (mesh_->face_centroid(f))[d - 1];
        rhs_c[0][cells[0]] += factor * Aface(0, 0) * (zc - zf);
      } else {
        double factor = rho_ * norm(g_);
        v.Reshape(ndofs);
        av.Reshape(ndofs);

        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        for (int n = 0; n < ndofs; ++n) {
          double zc = (mesh_->cell_centroid(cells[n]))[d - 1];
          v(n) = zc * factor;
        }

        Aface.Multiply(v, av, false);

        for (int n = 0; n < ndofs; ++n) { rhs_c[0][cells[n]] += av(n); }
      }
    }
  }

  global_op_->rhs()->GatherGhostedToMaster("cell", Add);
}


/* ******************************************************************
* Special implementation of boundary conditions.
****************************************************************** */
void
PDE_DiffusionFVonManifolds::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  AMANZI_ASSERT(bcs_trial_.size() > 0);
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  Teuchos::RCP<const Epetra_MultiVector> k_f = Teuchos::null;
  if (k_ != Teuchos::null) k_f = k_->ViewComponent("face", true);

  auto& rhs_c = *global_op_->rhs()->ViewComponent("cell", true);

  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] != OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int c = cells[0];

      if (bc_model[f] == OPERATOR_BC_DIRICHLET && primary) {
        rhs_c[0][c] += bc_value[f] * local_op_->matrices[f](0, 0);
      } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        local_op_->matrices_shadow[f] = local_op_->matrices[f];
        local_op_->matrices[f](0, 0) = 0.0;

        if (primary) rhs_c[0][c] -= bc_value[f] * mesh_->face_area(f);
      }
    }
  }
}


/* ******************************************************************
* Calculate mass flux from cell-centered data
****************************************************************** */
void
PDE_DiffusionFVonManifolds::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& solution,
                                       const Teuchos::Ptr<CompositeVector>& mass_flux)
{
  const auto& beta_f = *beta_->ViewComponent("face", true);
  const auto& fmap = *beta_->Map().Map("face", true);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  solution->ScatterMasterToGhosted("cell");

  Teuchos::RCP<const Epetra_MultiVector> k_f = Teuchos::null;
  if (k_ != Teuchos::null) {
    k_->ScatterMasterToGhosted("face");
    k_f = k_->ViewComponent("face", true);
  }

  solution->ScatterMasterToGhosted("cell");
  const auto& p = *solution->ViewComponent("cell", true);
  auto& flux = *mass_flux->ViewComponent("face", false);

  int dir, d(mesh_->space_dimension());
  AmanziMesh::Entity_ID_List cells;
  WhetStone::DenseVector ti(2), pi(2);

  for (int f = 0; f < nfaces_owned; ++f) {
    int g = fmap.FirstPointInElement(f);
    int ndofs = fmap.ElementSize(f);

    if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      double value = bc_value[f];
      int c = getFaceOnBoundaryInternalCell(*mesh_, f);
      getFaceNormalExterior(*mesh_, f, &dir);

      double tmp = p[0][c] - value;
      if (gravity_) {
        double factor = rho_ * norm(g_);
        double zc = (mesh_->cell_centroid(c))[d - 1];
        double zf = (mesh_->face_centroid(f))[d - 1];
        tmp += factor * (zf - zc);
      }

      flux[0][g] = dir * beta_f[0][g] * tmp;
      if (k_f.get()) flux[0][g] *= (*k_f)[0][g];

    } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
      double value = bc_value[f];
      double area = mesh_->face_area(f);
      getFaceNormalExterior(*mesh_, f, &dir);
      flux[0][g] = dir * value * area;

    } else {
      double factor = rho_ * norm(g_);
      ti.Reshape(ndofs);
      pi.Reshape(ndofs);

      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

      double sum(0.0), pf(0.0);
      for (int i = 0; i != ndofs; ++i) {
        ti(i) = beta_f[0][g + i] * (k_f.get() ? (*k_f)[0][g + i] : 1.0);
        sum += ti(i);

        double zc = (mesh_->cell_centroid(cells[i]))[d - 1];
        pi(i) = p[0][cells[i]] - factor * zc;
        pf += ti(i) * pi(i);
      }
      if (sum > 0.0) pf /= sum;

      for (int i = 0; i < ndofs; ++i) { flux[0][g + i] = -ti(i) * (pf - pi(i)); }
    }
  }
}


/* ******************************************************************
* Compute transmissibilities on faces 
****************************************************************** */
void
PDE_DiffusionFVonManifolds::ComputeBeta_()
{
  auto& beta_f = *beta_->ViewComponent("face", true);
  const auto& fmap = *beta_->Map().Map("face", true);

  int d = mesh_->space_dimension();
  AmanziGeometry::Point a(d);
  AmanziMesh::Entity_ID_List cells;

  auto k_f = (k_.get()) ? k_->ViewComponent("face") : Teuchos::null;
  ;
  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;

  for (int f = 0; f < nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    int g = fmap.FirstPointInElement(f);
    for (int i = 0; i < ncells; i++) {
      int c = cells[i];
      a = mesh_->face_centroid(f) - mesh_->cell_centroid(c);
      beta_f[0][g + i] = mesh_->face_area(f) / norm(a);

      if (K_.get()) {
        const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c);
        double perm = (((*K_)[c] * a) * normal);
        double dxn = a * normal;
        beta_f[0][g + i] *= perm / dxn;
      } else if (k_f.get()) {
        beta_f[0][g + i] *= (*k_f)[0][g + i];
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi
