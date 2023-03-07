/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

// Amanzi
#include "Mesh_Algorithms.hh"
#include "MeshFactory.hh"
#include "nlfv.hh"
#include "ParallelCommunication.hh"

#include "Op_Face_Cell.hh"
#include "OperatorDefs.hh"
#include "Operator_Cell.hh"
#include "PDE_DiffusionNLFVonManifolds.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Compute harmonic averaging points (function of geometry and tensor)
* and the positive decomposition of face conormals. The face-based
* data from left and right cells are ordered by the global cells ids.
****************************************************************** */
void
PDE_DiffusionNLFVonManifolds::InitStencils_()
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();

  // allocate persistent memory
  auto cvs_aux = CreateManifoldCVS(mesh_);
  auto mmap = cvs_aux->Map("face", false);
  auto gmap = cvs_aux->Map("face", true);

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("hap", AmanziMesh::FACE, dim_)
    ->AddComponent("gamma", AmanziMesh::FACE, mmap, gmap, 1)
    ->AddComponent("weight", AmanziMesh::FACE, mmap, gmap, dim_ - 1)
    ->AddComponent("flux_data", AmanziMesh::FACE, mmap, gmap, dim_ - 1);
  stencil_data_ = Teuchos::rcp(new CompositeVector(cvs));

  Epetra_MultiVector& hap = *stencil_data_->ViewComponent("hap", true);
  Epetra_MultiVector& gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);

  stencil_data_->PutScalarMasterAndGhosted(0.0);

  stencil_faces_.resize(dim_ - 1);
  stencil_cells_.resize(dim_ - 1);
  for (int i = 0; i < dim_ - 1; ++i) {
    stencil_faces_[i] = Teuchos::rcp(new Epetra_IntVector(*gmap));
    stencil_cells_[i] = Teuchos::rcp(new Epetra_IntVector(*gmap));

    stencil_faces_[i]->PutValue(0);
    stencil_cells_[i]->PutValue(0);
  }

  // allocate temporary memory for distributed tensor
  int dim2 = dim_ * dim_;
  CompositeVectorSpace cvs_tmp;
  cvs_tmp.SetMesh(mesh_)->SetGhosted(true)->AddComponent("tensor", AmanziMesh::CELL, dim2);
  Teuchos::RCP<CompositeVector> cv_tmp = Teuchos::rcp(new CompositeVector(cvs_tmp));
  Epetra_MultiVector& Ktmp = *cv_tmp->ViewComponent("tensor", true);

  // instantiate variables to access supporting tools
  WhetStone::NLFV nlfv(mesh_);

  // distribute diffusion tensor
  WhetStone::DenseVector data(dim2);
  for (int c = 0; c < ncells_owned; ++c) {
    WhetStone::TensorToVector((*K_)[c], data);
    for (int i = 0; i < dim2; ++i) { Ktmp[i][c] = data(i); }
  }
  cv_tmp->ScatterMasterToGhosted();

  // calculate harmonic averaging points (HAPs)
  int dir;
  std::vector<double> hap_weight;
  WhetStone::Tensor T(dim_, 2);
  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point p(dim_);
  std::vector<double> Tf;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    Tf.resize(ncells);

    OrderCellsByGlobalId(*mesh_, cells);
    if (ncells > 1) {
      for (int n = 0; n < ncells; ++n) {
        int c = cells[n];
        const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c, &dir);

        // create data for computing HAP
        for (int i = 0; i < dim2; ++i) data(i) = Ktmp[i][c];
        VectorToTensor(data, T);
        Tf[n] = T(0, 0);
      }
      nlfv.HarmonicAveragingPoint(f, cells, Tf, p, hap_weight);
    } else {
      p = mesh_->face_centroid(f);
      hap_weight[0] = 0.0;
    }

    // factor going to stencil should be (1 - weight)
    int g = gmap->FirstPointInElement(f);
    int ndofs = gmap->ElementSize(f);
    for (int i = 0; i < ndofs; ++i) gamma[0][g + i] = 1.0 - hap_weight[i];

    for (int i = 0; i < dim_; ++i) hap[i][f] = p[i];
  }

  stencil_data_->ScatterMasterToGhosted("hap");
  stencil_data_->ScatterMasterToGhosted("gamma");

  // calculate coefficients in positive decompositions of conormals
  AmanziGeometry::Point conormal(dim_), v(dim_);
  std::vector<AmanziGeometry::Point> tau;

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    // calculate list of candidate vectors
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    tau.clear();
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c, &dir);
        v = (*K_)[c] * normal;
      } else {
        for (int i = 0; i < dim_; ++i) v[i] = hap[i][f] - xc[i];
      }
      tau.push_back(v);
    }

    // decompose co-normals (we need memory allocation like in 3D)
    int ierr, ids[dim_];
    double ws[dim_];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c, &dir);
      conormal = (*K_)[c] * normal;

      ierr = nlfv.PositiveDecompositionManifold(n, tau, conormal, ws, ids);
      AMANZI_ASSERT(ierr == 0);

      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      OrderCellsByGlobalId(*mesh_, cells);

      int g = gmap->FirstPointInElement(f);
      int ndofs = gmap->ElementSize(f);
      int pos(0);
      for (int n = 0; n < ndofs; ++n) {
        if (cells[n] == c) {
          pos = n; 
          break;
        }
      }

      for (int i = 0; i < dim_ - 1; i++) {
        weight[i][g + pos] = ws[i];
        (*stencil_faces_[i])[g + pos] = faces[ids[i]];
        (*stencil_cells_[i])[g + pos] = cell_get_face_adj_cell(*mesh_, c, faces[ids[i]]);
      }
    }
  }

  // distribute stencils
  stencil_data_->GatherGhostedToMaster("weight");
  stencil_data_->ScatterMasterToGhosted("weight");

  ParallelCommunication pp(mesh_);
  for (int i = 0; i < dim_ - 1; ++i) {
    pp.CombineGhostFace2MasterFace(*stencil_faces_[i], (Epetra_CombineMode)Add);
    pp.CombineGhostFace2MasterFace(*stencil_cells_[i], (Epetra_CombineMode)Add);

    pp.CopyMasterFace2GhostFace(*stencil_faces_[i]);
    pp.CopyMasterFace2GhostFace(*stencil_cells_[i]);
  }

  stencil_initialized_ = true;
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces and 1x1 matrices
* on boundary faces. We avoid round-off operations since the stencils
* already incorporate them.
****************************************************************** */
void
PDE_DiffusionNLFVonManifolds::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                             const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!stencil_initialized_) InitStencils_();
  if (k_ != Teuchos::null) k_->ScatterMasterToGhosted("face");

  u->ScatterMasterToGhosted("cell");

  Epetra_MultiVector& hap_gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);
  Epetra_MultiVector& flux_data = *stencil_data_->ViewComponent("flux_data", true);

  // allocate auxiliary matrix structure
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 2);
  CompositeVector matrix_cv(cvs), sideflux_cv(cvs);

  Epetra_MultiVector& matrix = *matrix_cv.ViewComponent("face", true);
  Epetra_MultiVector& sideflux = *sideflux_cv.ViewComponent("face", true);

  // calculate one-sides flux corrections. Since a flux stencil can
  // occupy (dim_ + 1) cells, we need parallel communications.
  OneSidedFluxCorrections_(1, *u, sideflux_cv);

  // un-rolling little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) k_face = k_->ViewComponent("face");

  // split each stencil between different local matrices
  int c1, c2, c3, c4, k1, k2;
  AmanziMesh::Entity_ID_List cells, cells_tmp;

  matrix_cv.PutScalarMasterAndGhosted(0.0);
  flux_data.PutScalar(0.0);

  for (int c = 0; c < ncells_owned; ++c) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      OrderCellsByGlobalId(*mesh_, cells, c1, c2);
      k1 = (c1 == c) ? 0 : 1;
      k2 = k1 * dim_;

      // calculate little_k on the current face
      double kf(1.0);
      if (k_face.get()) kf = (*k_face)[0][f];

      // Calculate solution-dependent weigths using corrections to the
      // two-point flux. Note mu does not depend on one-sided flux.
      double gamma, g1, g2, gg(-1.0), w1, w2(0.0), tpfa, mu(1.0);

      gamma = hap_gamma[0][f];
      if (ncells == 2) {
        w1 = weight[0][f] * gamma;
        w2 = weight[dim_][f] * (1.0 - gamma);

        g1 = sideflux[0][f];
        g2 = sideflux[1][f];
        gg = g1 * g2;


        g1 = fabs(g1);
        g2 = fabs(g2);
        mu = (g1 + g2 == 0.0) ? 0.5 : g2 / (g1 + g2);
      } else {
        w1 = weight[0][f];
      }

      tpfa = mu * w1 + (1.0 - mu) * w2;
      matrix[k1][f] += kf * tpfa;
      flux_data[k2][f] = kf * tpfa;

      // remaining terms of one-sided flux in cell c. Now we need
      // to select mu depending on the one-sided flux.
      if (gg < 0.0) {
        if (c1 != c) mu = 1.0 - mu;
        for (int i = 1; i < dim_; i++) {
          int f1 = (*stencil_faces_[i + k2])[f];
          mesh_->face_get_cells(f1, AmanziMesh::Parallel_type::ALL, &cells_tmp);

          gamma = hap_gamma[0][f1];
          OrderCellsByGlobalId(*mesh_, cells_tmp, c3, c4);

          k1 = 0;
          if (c3 != c) {
            gamma = 1.0 - gamma;
            k1 = 1;
          }

          double tmp = ncells * weight[i + k2][f] * gamma * mu;
          matrix[k1][f1] += kf * tmp;
          flux_data[k2 + i][f] = kf * tmp;
        }
      }
    }
  }

  //stencil_data_->ScatterMasterToGhosted("flux_data");
  stencil_data_->GatherGhostedToMaster("flux_data");
  matrix_cv.GatherGhostedToMaster();

  // populate local matrices
  for (int f = 0; f < nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    WhetStone::DenseMatrix Aface(ncells, ncells);

    if (ncells == 2) {
      k1 = OrderCellsByGlobalId(*mesh_, cells, c3, c4);
      k2 = 1 - k1;
      Aface(0, 0) = matrix[k1][f];
      Aface(0, 1) = -matrix[k1][f];

      Aface(1, 0) = -matrix[k2][f];
      Aface(1, 1) = matrix[k2][f];
    } else {
      Aface(0, 0) = matrix[0][f];
    }

    local_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Modify operator by adding an upwind approximation of the Newton
* correction term.
****************************************************************** */
void
PDE_DiffusionNLFVonManifolds::UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    double scalar_factor)
{
  // ignore correction if no flux provided.
  if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

  if (k_->HasComponent("face")) k_->ScatterMasterToGhosted("face");
  if (dkdp_->HasComponent("face")) dkdp_->ScatterMasterToGhosted("face");

  // Correction is not required
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_NONE) return;

  // only works on upwinded methods
  if (little_k_ == OPERATOR_UPWIND_NONE) return;

  const Epetra_MultiVector& kf = *k_->ViewComponent("face");
  const Epetra_MultiVector& dkdp_f = *dkdp_->ViewComponent("face");
  const Epetra_MultiVector& flux_f = *flux->ViewComponent("face");

  // populate the local matrices
  double v, vmod;
  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface.PutScalar(0.0);

    // We use the upwind discretization of the generalized flux.
    v = std::abs(kf[0][f]) > 0.0 ? flux_f[0][f] * dkdp_f[0][f] / kf[0][f] : 0.0;
    vmod = std::abs(v);

    // prototype for future limiters
    vmod *= scalar_factor;

    // We use the upwind discretization of the generalized flux.
    int i, dir, c1;
    c1 = cells[0];
    mesh_->face_normal(f, false, c1, &dir);
    i = (v * dir >= 0.0) ? 0 : 1;

    if (ncells == 2) {
      Aface(i, i) = vmod;
      Aface(1 - i, i) = -vmod;
    } else if (i == 0) {
      Aface(0, 0) = vmod;
    }

    jac_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Modify operator by adding an upwind approximation of the Newton
* correction term.
****************************************************************** */
void
PDE_DiffusionNLFVonManifolds::UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& factor)
{
  // ignore correction if no flux provided.
  if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

  k_->ScatterMasterToGhosted("face");
  dkdp_->ScatterMasterToGhosted("face");

  // Correction is not required
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_NONE) return;

  // only works on upwinded methods
  if (little_k_ == OPERATOR_UPWIND_NONE) return;

  const Epetra_MultiVector& kf = *k_->ViewComponent("face");
  const Epetra_MultiVector& dkdp_f = *dkdp_->ViewComponent("face");
  const Epetra_MultiVector& flux_f = *flux->ViewComponent("face");
  const Epetra_MultiVector& factor_cell = *factor->ViewComponent("cell");

  // populate the local matrices
  double v, vmod;
  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface.PutScalar(0.0);

    // We use the upwind discretization of the generalized flux.
    v = std::abs(kf[0][f]) > 0.0 ? flux_f[0][f] * dkdp_f[0][f] / kf[0][f] : 0.0;
    vmod = std::abs(v);

    double scalar_factor = 0.;
    for (int j = 0; j < ncells; j++) scalar_factor += factor_cell[0][cells[j]];
    scalar_factor *= 1. / ncells;

    // prototype for future limiters
    vmod *= scalar_factor;

    // We use the upwind discretization of the generalized flux.
    int i, dir, c1;
    c1 = cells[0];
    mesh_->face_normal(f, false, c1, &dir);
    i = (v * dir >= 0.0) ? 0 : 1;

    if (ncells == 2) {
      Aface(i, i) = vmod;
      Aface(1 - i, i) = -vmod;
    } else if (i == 0) {
      Aface(0, 0) = vmod;
    }

    jac_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Calculate one-sided fluxes (i0=0) or flux corrections (i0=1).
****************************************************************** */
void
PDE_DiffusionNLFVonManifolds::OneSidedFluxCorrections_(int i0,
                                                       const CompositeVector& u,
                                                       CompositeVector& flux_cv)
{
  // un-rolling composite vectors
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);
  Epetra_MultiVector& flux = *flux_cv.ViewComponent("face", true);

  Epetra_MultiVector& hap_gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);

  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();

  // un-rolling little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) k_face = k_->ViewComponent("face");

  int c1, c2, c3, k1, k2;
  double gamma, tmp;
  AmanziMesh::Entity_ID_List cells, cells_tmp;

  flux_cv.PutScalarMasterAndGhosted(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

      OrderCellsByGlobalId(*mesh_, cells, c1, c2);
      k1 = (c1 == c) ? 0 : 1;
      k2 = k1 * dim_;

      // scalar (nonlinear) coefficient
      double kf(1.0);
      if (k_face.get()) kf = (*k_face)[0][f];

      double sideflux(0.0), neumann_flux(0.0);
      for (int i = i0; i < dim_; ++i) {
        int f1 = (*stencil_faces_[i + k2])[f];
        c3 = (*stencil_cells_[i + k2])[f];
        if (c3 >= 0) {
          mesh_->face_get_cells(f1, AmanziMesh::Parallel_type::ALL, &cells_tmp);
          OrderCellsByGlobalId(*mesh_, cells_tmp, c1, c2);

          gamma = hap_gamma[0][f1];
          if (c1 != c) gamma = 1.0 - gamma;

          tmp = weight[i + k2][f] * gamma;
          sideflux += tmp * (uc[0][c] - uc[0][c3]);
        } else if (bc_model[f1] == OPERATOR_BC_DIRICHLET) {
          tmp = weight[i + k2][f];
          // mesh_->face_normal(f1, false, c, &dir);
          sideflux += tmp * (uc[0][c] - MapBoundaryValue_(f1, bc_value[f1])); // * dir;
        } else if (bc_model[f1] == OPERATOR_BC_NEUMANN) {
          tmp = weight[i + k2][f];
          neumann_flux += tmp * bc_value[f1] * mesh_->face_area(f1);
        }
      }

      flux[k1][f] = kf * sideflux + neumann_flux;
    }
  }

  flux_cv.GatherGhostedToMaster();
  flux_cv.ScatterMasterToGhosted();
}


/* ******************************************************************
* Calculate one-sided fluxes (i0=0) or flux corrections (i0=1).
****************************************************************** */
void
PDE_DiffusionNLFVonManifolds::OneSidedWeightFluxes_(
    int i0, const CompositeVector& u, CompositeVector& flux_cv)
{
  // un-rolling composite vectors
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);
  Epetra_MultiVector& flux = *flux_cv.ViewComponent("face", true);

  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  Epetra_MultiVector& flux_data = *stencil_data_->ViewComponent("flux_data", true);

  int c1, c2, c3, k1, k2;
  AmanziMesh::Entity_ID_List cells, cells_tmp;

  flux_cv.PutScalarMasterAndGhosted(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

      OrderCellsByGlobalId(*mesh_, cells, c1, c2);
      k1 = (c1 == c) ? 0 : 1;
      k2 = k1 * dim_;

      double sideflux(0.0), neumann_flux(0.0);
      for (int i = i0; i < dim_; ++i) {
        c3 = (*stencil_cells_[i + k2])[f];
        int f1 = (*stencil_faces_[i + k2])[f];

        if (c3 >= 0) {
          sideflux += flux_data[i + k2][f] * (uc[0][c] - uc[0][c3]);
        } else if (bc_model[f1] == OPERATOR_BC_DIRICHLET) {
          sideflux += flux_data[i + k2][f] * (uc[0][c] - MapBoundaryValue_(f1, bc_value[f1]));
        } else if (bc_model[f1] == OPERATOR_BC_NEUMANN) {
          neumann_flux += flux_data[i + k2][f] * bc_value[f1] * mesh_->face_area(f1);
        }
      }
      flux[k1][f] = sideflux + neumann_flux;
    }
  }

  flux_cv.GatherGhostedToMaster();
  flux_cv.ScatterMasterToGhosted();
}


/* ******************************************************************
* Matrix-based implementation of boundary conditions.
****************************************************************** */
void
PDE_DiffusionNLFVonManifolds::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell", true);

  // un-rolling little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) k_face = k_->ViewComponent("face");

  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] != OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int c = cells[0];

      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
        rhs_cell[0][c] += Aface(0, 0) * bc_value[f];
      } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
        local_op_->matrices_shadow[f] = Aface;

        double kf(1.0);
        if (k_face.get()) kf = (*k_face)[0][f];

        rhs_cell[0][c] -= (Aface(0, 0) / kf) * bc_value[f] * mesh_->face_area(f);
        Aface = 0.0;
      }
    }
  }

  return;
}


/* ******************************************************************
* Calculate flux using cell-centered data.
* **************************************************************** */
void
PDE_DiffusionNLFVonManifolds::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                         const Teuchos::Ptr<CompositeVector>& flux)
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  // const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  if (k_ != Teuchos::null) k_->ScatterMasterToGhosted("face");

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 2);
  CompositeVector wgt_sideflux_cv(cvs);

  Epetra_MultiVector& wgt_sideflux = *wgt_sideflux_cv.ViewComponent("face", true);
  Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);

  u->ScatterMasterToGhosted("cell");
  // OneSidedFluxCorrections_(0, u, sideflux_cv);
  OneSidedWeightFluxes_(0, *u, wgt_sideflux_cv);

  int c1, c2, dir;
  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; ++f) {
    if (bc_model[f] == OPERATOR_BC_DIRICHLET || bc_model[f] == OPERATOR_BC_NEUMANN) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      mesh_->face_normal(f, false, cells[0], &dir);
      flux_data[0][f] = wgt_sideflux[0][f] * dir;
      // } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
      //   flux_data[0][f] = bc_value[f] * mesh_->face_area(f);
    } else if (bc_model[f] == OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      OrderCellsByGlobalId(*mesh_, cells, c1, c2);
      mesh_->face_normal(f, false, c1, &dir);

      double wg1 = wgt_sideflux[0][f];
      double wg2 = wgt_sideflux[1][f];

      if (cells.size() == 2)
        flux_data[0][f] = 0.5 * (wg1 - wg2) * dir;
      else
        flux_data[0][f] = dir * wg1;
    }
  }
}

} // namespace Operators
} // namespace Amanzi
