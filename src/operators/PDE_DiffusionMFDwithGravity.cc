// PDE_DiffusionMFDwithGravity prescribes an elliptic operator with gravity using MFD family of discretizations.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "OperatorDefs.hh"
#include "UniqueLocalIndex.hh"

#include "PDE_DiffusionMFDwithGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Add a gravity term to the diffusion operator.
****************************************************************** */
void PDE_DiffusionMFDwithGravity::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  // optimize div[K k(grad u + b g u)] to div[K1 (grad u + b1 g u)]
  if (!gravity_term_initialized_) {
    if (is_scalar_) rho_ *= const_k_;
    gravity_term_initialized_ = true;
  }

  PDE_DiffusionMFD::UpdateMatrices(flux, u);
  AddGravityToRHS_();
}


/* ******************************************************************
* Add a gravity term to the RHS of the operator
****************************************************************** */
void PDE_DiffusionMFDwithGravity::AddGravityToRHS_()
{
  // vector or scalar rho?
  const Epetra_MultiVector* rho_c = NULL;
  if (!is_scalar_) rho_c = &*rho_cv_->ViewComponent("cell", false);

  if (global_op_->rhs()->HasComponent("face")) {
    int dim = mesh_->space_dimension();

    // preparing little-k data
    Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
    Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
    if (k_ != Teuchos::null) {
      if (k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");
      if (k_->HasComponent("face")) k_face = k_->ViewComponent("face", true);
      if (k_->HasComponent("grav")) k_face = k_->ViewComponent("grav", true);
    }

    int dir;

    Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");
    Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);
    for (int f = nfaces_owned; f < nfaces_wghost; f++) rhs_face[0][f] = 0.0;

    WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
    Kc(0, 0) = 1.0;
    if (const_K_.rank() > 0) Kc = const_K_;

    // gravity discretization
    bool fv_flag = (gravity_method_ == OPERATOR_GRAVITY_FV) ||
        !(little_k_ & OPERATOR_LITTLE_K_DIVK_BASE);

    for (int c = 0; c < ncells_owned; c++) {
      const auto& faces = mesh_->cell_get_faces(c);
      const auto& dirs = mesh_->cell_get_face_dirs(c);
      int nfaces = faces.size();
      double zc = (mesh_->cell_centroid(c))[dim - 1];

      // building blocks for the gravity term
      double rho = rho_c ? (*rho_c)[0][c] : rho_;
      WhetStone::DenseMatrix& Wff = Wff_cells_[c];

      // Update terms due to nonlinear coefficient
      double kc(1.0);
      std::vector<double> kf(nfaces, 1.0);
      // -- chefs recommendation: SPD discretization with upwind
      if (little_k_ == OPERATOR_LITTLE_K_DIVK ||
          little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN) {
        kc = (*k_cell)[0][c];
        for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];

      // -- new scheme: SPD discretization with upwind and equal spliting
      } else if (little_k_ == OPERATOR_LITTLE_K_DIVK_BASE) {
        if (!fv_flag) {
          for (int n = 0; n < nfaces; n++) kf[n] = std::sqrt((*k_face)[0][faces[n]]);
        } else {
          for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
        }

      // -- the second most popular choice: classical upwind
      } else if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
        for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];

      } else if (little_k_ == OPERATOR_LITTLE_K_STANDARD && k_cell != Teuchos::null) {
        kc = (*k_cell)[0][c];
        for (int n = 0; n < nfaces; n++) kf[n] = kc;
      }

      // add gravity term to the right-hand side vector.
      // -- use always for the finite volume method
      if (fv_flag) { 
        if (K_.get()) Kc = (*K_)[c];
        AmanziGeometry::Point Kcg(Kc * g_);

        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c, &dir);
          double tmp;

          if (gravity_special_projection_) {
            const AmanziGeometry::Point& xcc = GravitySpecialDirection_(f);
            double sign = (normal * xcc) * dir;

            tmp = (Kcg * xcc) * rho * kf[n] * dirs[n];
            tmp *= copysign(norm(normal) / norm(xcc), sign);
          } else {
            tmp = (Kcg * normal) * rho * kf[n];
          }

          rhs_face[0][f] += tmp; 
          rhs_cell[0][c] -= tmp; 
        }
      }

      // Amanzi's first upwind: the family of DIVK methods uses hydraulic
      // head as the primary variable and linear transformation for pressure.
      if (!fv_flag) {
        WhetStone::DenseVector v(nfaces), av(nfaces);
        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          double zf = (mesh_->face_centroid(f))[dim - 1];
          v(n) = -(zf - zc) * kf[n] * rho * norm(g_) / kc;
        }

        Wff.Multiply(v, av, false);

        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          double tmp = av(n) * kf[n];

          rhs_face[0][f] += tmp; 
          rhs_cell[0][c] -= tmp; 
        }
      }
    }

    global_op_->rhs()->GatherGhostedToMaster("face", Epetra_CombineMode(Add));
  }
}


/* ******************************************************************
* WARNING: Since gravity flux is not continuous, we derive it in 
* exactly the same manner as in other routines.
* **************************************************************** */
void PDE_DiffusionMFDwithGravity::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                             const Teuchos::Ptr<CompositeVector>& flux)
{
  // Calculate diffusive part of the flux.
  PDE_DiffusionMFD::UpdateFlux(u, flux);

  // vector or scalar rho?
  const Epetra_MultiVector* rho_c = NULL;
  if (!is_scalar_) rho_c = &*rho_cv_->ViewComponent("cell", false);

  // preparing little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) {
    if (k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");
    if (k_->HasComponent("face")) k_face = k_->ViewComponent("face", true);
    if (k_->HasComponent("grav")) k_face = k_->ViewComponent("grav", true);
  }

  int dim = mesh_->space_dimension();
  Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);
  Epetra_MultiVector grav_flux(flux_data);
  grav_flux.PutScalar(0.0);

  std::vector<int> hits(nfaces_wghost, 0);

  WhetStone::Tensor Kc(dim, 1);
  Kc(0, 0) = 1.0;
  if (const_K_.rank() > 0) Kc = const_K_;

  // gravity discretization
  bool fv_flag = (gravity_method_ == OPERATOR_GRAVITY_FV) ||
      !(little_k_ & OPERATOR_LITTLE_K_DIVK_BASE);

  for (int c = 0; c < ncells_owned; c++) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();
    double zc = mesh_->cell_centroid(c)[dim - 1];

    // building blocks for the gravity term
    double rho = rho_c ? (*rho_c)[0][c] : rho_;
    WhetStone::DenseMatrix& Wff = Wff_cells_[c];

    // Update terms due to nonlinear coefficient
    double kc(1.0);
    std::vector<double> kf(nfaces, 1.0);
    if (little_k_ == OPERATOR_LITTLE_K_DIVK) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    } else if (little_k_ == OPERATOR_LITTLE_K_DIVK_BASE) {
      for (int n = 0; n < nfaces; n++) kf[n] = std::sqrt((*k_face)[0][faces[n]]);
    } else if (little_k_ == OPERATOR_LITTLE_K_STANDARD && k_cell != Teuchos::null) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    } else if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    }

    if (fv_flag) { 
      if (K_.get()) Kc = (*K_)[c];
      AmanziGeometry::Point Kcg(Kc * g_);

      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        if (f < nfaces_owned) {
          const AmanziGeometry::Point& normal = mesh_->face_normal(f);

          if (gravity_special_projection_) {
            const AmanziGeometry::Point& xcc = GravitySpecialDirection_(f);
            double sign = normal * xcc;
            double tmp = copysign(norm(normal) / norm(xcc), sign);
            grav_flux[0][f] += (Kcg * xcc) * rho * kf[n] * tmp;
          } else {
            grav_flux[0][f] += (Kcg * normal) * rho * kf[n];
          }
            
          hits[f]++;
        }
      }
    }

    if (!fv_flag) { 
      WhetStone::DenseVector v(nfaces), av(nfaces);
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        double zf = (mesh_->face_centroid(f))[dim - 1];
        v(n) = -(zf - zc) * kf[n] * rho * norm(g_) / kc;
      }

      Wff.Multiply(v, av, false);

      for (int n = 0; n < nfaces; n++) {
        int dir, f = faces[n];
        if (f < nfaces_owned) {
          mesh_->face_normal(f, false, c, &dir);
            
          double tmp = av(n) * kf[n] * dir;
          grav_flux[0][f] += tmp;

          hits[f]++;
        }
      }
    }
  }

  for (int f=0; f < nfaces_owned; f++) {
    flux_data[0][f] += grav_flux[0][f] / hits[f];
  }
}


/* ******************************************************************
* Add "gravity flux" to the Darcy flux. 
* **************************************************************** */
void PDE_DiffusionMFDwithGravity::UpdateFluxNonManifold(
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<CompositeVector>& flux)
{
  // Calculate diffusive part of the flux.
  PDE_DiffusionMFD::UpdateFluxNonManifold(u, flux);

  // preparing little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) {
    if (k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");
    if (k_->HasComponent("face")) k_face = k_->ViewComponent("face", true);
    if (k_->HasComponent("grav")) k_face = k_->ViewComponent("grav", true);
  }

  int dim = mesh_->space_dimension();
  Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);

  CompositeVector grav(*flux);
  Epetra_MultiVector& grav_data = *grav.ViewComponent("face", true);
  grav_data.PutScalar(0.0);

  int ndofs_owned = flux->ViewComponent("face")->MyLength();

  const auto& fmap = *flux->Map().Map("face", true);

  WhetStone::Tensor Kc(dim, 1);
  Kc(0, 0) = 1.0;
  if (const_K_.rank() > 0) Kc = const_K_;

  for (int c = 0; c < ncells_owned; c++) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    // Update terms due to nonlinear coefficient
    double kc(1.0);
    std::vector<double> kf(nfaces, 1.0);
    if (little_k_ == OPERATOR_LITTLE_K_DIVK) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    } else if (little_k_ == OPERATOR_LITTLE_K_DIVK_BASE) {
      for (int n = 0; n < nfaces; n++) kf[n] = std::sqrt((*k_face)[0][faces[n]]);
    } else if (little_k_ == OPERATOR_LITTLE_K_STANDARD && k_cell != Teuchos::null) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    } else if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    }

    if (K_.get()) Kc = (*K_)[c];
    AmanziGeometry::Point Kcg(Kc * g_);

    for (int n = 0; n < nfaces; n++) {
      int dir, f = faces[n];
      AmanziGeometry::Point normal = mesh_->face_normal(f, false, c, &dir);
      normal *= dir;

      int g = fmap.FirstPointInElement(f);
      int ndofs = fmap.ElementSize(f);
      if (ndofs > 1) g += Operators::UniqueIndexFaceToCells(*mesh_, f, c);

      if (gravity_special_projection_) {
        const AmanziGeometry::Point& xcc = GravitySpecialDirection_(f);
        double sign = normal * xcc;
        double tmp = copysign(norm(normal) / norm(xcc), sign);
        grav_data[0][g] += (Kcg * xcc) * rho_ * kf[n] * tmp;
      } else {
        grav_data[0][g] += (Kcg * normal) * rho_ * kf[n];
      }
    }
  }

  // if f is on a processor boundary, some g are not initialized
  grav.GatherGhostedToMaster(Add);

  for (int g = 0; g < ndofs_owned; ++g) {
    flux_data[0][g] += grav_data[0][g];
  }
}


/* ******************************************************************
* Put here stuff that has to be done in constructor, i.e. only once.
****************************************************************** */
void PDE_DiffusionMFDwithGravity::Init_(Teuchos::ParameterList& plist)
{
  gravity_special_projection_ = (mfd_primary_ == WhetStone::DIFFUSION_TPFA);

  // gravity discretization
  std::string name = plist.get<std::string>("gravity term discretization", "hydraulic head");
  if (name == "hydraulic head")
    gravity_method_ = OPERATOR_GRAVITY_HH;
  else
    gravity_method_ = OPERATOR_GRAVITY_FV;

  gravity_term_initialized_ = false;
}


/* ******************************************************************
* Compute non-normalized unsigned direction to the next cell needed
* to project gravity vector in the MFD-TPFA discretization method.
****************************************************************** */
AmanziGeometry::Point PDE_DiffusionMFDwithGravity::GravitySpecialDirection_(int f) const
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  int ncells = cells.size();

  if (ncells == 2) {
    return mesh_->cell_centroid(cells[1]) - mesh_->cell_centroid(cells[0]);
  } else {
    return mesh_->face_centroid(f) - mesh_->cell_centroid(cells[0]);
  }
}


/* ******************************************************************
* Return value of the gravity flux on the given face f.
****************************************************************** */
double PDE_DiffusionMFDwithGravity::ComputeGravityFlux(int f) const
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  int c = cells[0];

  double gflux;
  const AmanziGeometry::Point& normal = mesh_->face_normal(f);

  if (K_.get()) {
    gflux = (((*K_)[c] * g_) * normal);
  } else if (const_K_.rank() > 0) {
    gflux = (const_K_ * g_) * normal;
  } else {
    gflux = g_ * normal;
  }

  if (is_scalar_) {
    gflux *= rho_;
  } else {
    const Epetra_MultiVector& rho_c = *rho_cv_->ViewComponent("cell", true);
    gflux *= rho_c[0][c];
  }

  return gflux;
}

}  // namespace Operators
}  // namespace Amanzi
