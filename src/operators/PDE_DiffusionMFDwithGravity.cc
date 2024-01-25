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
  AMANZI_ASSERT(!exclude_primary_terms_);
  PDE_DiffusionMFD::UpdateMatrices(flux, u);
  AddGravityToRHS_();
}


/* ******************************************************************
* Add a gravity term to the RHS of the operator
****************************************************************** */
void PDE_DiffusionMFDwithGravity::AddGravityToRHS_()
{
  if (global_op_->rhs()->hasComponent("face")) {
    global_op_->rhs()->putScalarGhosted(0.);
    int dim = mesh_->getSpaceDimension();

    Preallocate_little_k_(true);

    { // context for views
      auto rho_c = DensityCells(true);

      auto rhs_cell = global_op_->rhs()->viewComponent("cell");
      auto rhs_face = global_op_->rhs()->viewComponent("face", true);

      // gravity discretization
      bool fv_flag = (gravity_method_ == OPERATOR_GRAVITY_FV) ||
        !(little_k_type_ & OPERATOR_LITTLE_K_DIVK_BASE);

      const AmanziMesh::Mesh* mesh = mesh_.get();
      const TensorVector* K = K_.get();

      Kokkos::parallel_for(
          "PDE_DiffusionMFDwithGravity::AddGravityToRHS_",
          ncells_owned,
          KOKKOS_LAMBDA(const int& c) {
            auto [faces, dirs] = mesh->getCellFacesAndDirections(c);
            int nfaces = faces.size();

            // building blocks for the gravity term
            double zc = (mesh->getCellCentroid(c))[dim - 1];
            auto Wff = Wff_cells_[c];
            auto kr = kr_cells_[c];

            // MIGHT HAVE TO TWEAK KR?!?  See #439

            // add gravity term to the right-hand side vector.
            // -- always used for the finite volume method, also for non-DIVK
            //    upwind methods
            if (fv_flag) {
              AmanziGeometry::Point Kcg = K == nullptr ? g_ : K->at(c) * g_;

              for (int n = 0; n < nfaces; n++) {
                int f = faces[n];
                int dir;
                const AmanziGeometry::Point& normal = mesh->getFaceNormal(f, c, &dir);
                double zf = (mesh->getFaceCentroid(f))[dim - 1];
                double tmp;

                if (gravity_special_projection_) {
                  const AmanziGeometry::Point& xcc = Impl::GravitySpecialDirection(mesh,f);
                  double sign = (normal * xcc) * dir;

                  tmp = (Kcg * xcc) * rho_c(c,0) * kr(n) * dirs[n];
                  tmp *= copysign(norm(normal) / norm(xcc), sign);
                } else {
                  tmp = (Kcg * normal) * rho_c(c,0) * kr(n);
                }

                Kokkos::atomic_add(&rhs_face(f,0), tmp);
                rhs_cell(c,0) -= tmp;
              }
            }

            // Amanzi's first upwind: the family of DIVK methods uses hydraulic
            // head as the primary variable and linear transformation for pressure.
            if (!fv_flag) {
              assert(false); // crud -- we used v for kr, and now we
                                    // need all three v, Av, and kr?  Is there
                                    // a way around this?
              // WhetStone::DenseVector v(nfaces), av(nfaces);
              // for (int n = 0; n < nfaces; n++) {
              //   int f = faces[n];
              //   double zf = (mesh_->getFaceCentroid(f))[dim - 1];
              //   v(n) = -(zf - zc) * kf[n] * rho * norm(g_) / kc;
              // }

              // Wff.Multiply(v, av, false);

              // for (int n = 0; n < nfaces; n++) {
              //   int f = faces[n];
              //   double tmp = av(n) * kf[n];

              //   rhs_face[0][f] += tmp;
              //   rhs_cell[0][c] -= tmp;
              // }
            }
          });
    }
    global_op_->rhs()->gatherGhostedToMaster("face", Tpetra::ADD);
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
  int dim = mesh_->getSpaceDimension();

  CompositeVector grav_flux(flux->getMap());
  grav_flux.putScalar(0.0);

  Preallocate_little_k_(true);

  Kokkos::View<int*,DeviceOnlyMemorySpace> hits("hits",nfaces_owned);

  { // context for views
    auto rho_c = DensityCells(true);

    auto grav_flux_v = grav_flux.viewComponent("face", false);

    const AmanziMesh::Mesh* mesh = mesh_.get();
    const TensorVector* K = K_.get();

    // gravity discretization
    bool fv_flag = (gravity_method_ == OPERATOR_GRAVITY_FV) ||
                   !(little_k_type_ & OPERATOR_LITTLE_K_DIVK_BASE);

    Kokkos::parallel_for(
        "PDE_DiffusionMFDwithGravity::UpdateFlux",
        ncells_owned,
        KOKKOS_LAMBDA(const int& c) {
          auto faces = mesh->getCellFaces(c);
          int nfaces = faces.extent(0);

          // building blocks for the gravity term
          double zc = mesh->getCellCentroid(c)[dim - 1];
          auto Wff = Wff_cells_[c];
          auto kr = kr_cells_[c];

          if (fv_flag) {
            AmanziGeometry::Point Kcg = (K == nullptr) ? g_ : K->at(c) * g_;

            for (int n = 0; n < nfaces; n++) {
              int f = faces[n];
              if (f < nfaces_owned) {
                const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);

                if (gravity_special_projection_) {
                  const AmanziGeometry::Point& xcc = Impl::GravitySpecialDirection(mesh, f);
                  double sign = normal * xcc;
                  double tmp = copysign(norm(normal) / norm(xcc), sign);
                  Kokkos::atomic_add(&grav_flux_v(f,0), (Kcg * xcc) * rho_c(c,0) * kr(n) * tmp);
                } else {
                  Kokkos::atomic_add(&grav_flux_v(f,0), (Kcg * normal) * rho_c(c,0) * kr(n));
                }
                Kokkos::atomic_increment(&hits(f));
              }
            }
          }

          if (!fv_flag) {
            assert(false);
            // WhetStone::DenseVector v(nfaces), av(nfaces);
            // for (int n = 0; n < nfaces; n++) {
            //   int f = faces[n];
            //   double zf = (mesh->getFaceCentroid(f))[dim - 1];
            //   v(n) = -(zf - zc) * kf[n] * rho * norm(g_) / kc;
            // }

            // Wff.Multiply(v, av, false);

            // for (int n = 0; n < nfaces; n++) {
            //   int dir, f = faces[n];
            //   if (f < nfaces_owned) {
            //     const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c, &dir);

            //     double tmp = av(n) * kf[n] * dir;
            //     grav_flux[0][f] += tmp;

            //     hits[f]++;
            //   }
            // }
          }
        });
  }

  {
    auto grav_flux_v = grav_flux.viewComponent("face", false);
    auto flux_v = flux->viewComponent("face", false);
    Kokkos::parallel_for(
        "PDE_DiffusionMFDwithGravity::UpdateFlux Average",
        nfaces_owned,
        KOKKOS_LAMBDA(const int& f) {
          flux_v(f,0) += grav_flux_v(f,0) / hits(f);
        });
  }
}


// /* ******************************************************************
// * Add "gravity flux" to the Darcy flux.
// * **************************************************************** */
// void PDE_DiffusionMFDwithGravity::UpdateFluxNonManifold(
//     const Teuchos::Ptr<const CompositeVector>& u,
//     const Teuchos::Ptr<CompositeVector>& flux)
// {
//   // Calculate diffusive part of the flux.
//   PDE_DiffusionMFD::UpdateFluxNonManifold(u, flux);

//   // preparing little-k data
//   Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
//   Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
//   if (k_ != Teuchos::null) {
//     if (k_->hasComponent("cell")) k_cell = k_->viewComponent("cell");
//     if (k_->hasComponent("face")) k_face = k_->viewComponent("face", true);
//     if (k_->hasComponent("grav")) k_face = k_->viewComponent("grav", true);
//   }

//   int dim = mesh_->space_dimension();
//   Epetra_MultiVector& flux_data = *flux->viewComponent("face", true);

//   CompositeVector grav(*flux);
//   Epetra_MultiVector& grav_data = *grav.viewComponent("face", true);
//   grav_data.PutScalar(0.0);

//   int ndofs_owned = flux->viewComponent("face")->MyLength();
//   int ndofs_wghost = flux_data.MyLength();

//   const auto& fmap = *flux->getMap().getMap("face", true);
//   WhetStone::Tensor Kc(dim, 1);
//   Kc(0, 0) = 1.0;

//   for (int c = 0; c < ncells_owned; c++) {
//     mesh_->cell_get_faces(c, &faces);
//     int nfaces = faces.size();
//     double zc = mesh_->getCellCentroid(c)[dim - 1];

//     // Update terms due to nonlinear coefficient
//     double kc(1.0);
//     std::vector<double> kf(nfaces, 1.0);
//     if (little_k_type_ == OPERATOR_LITTLE_K_DIVK) {
//       kc = (*k_cell)[0][c];
//       for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
//     } else if (little_k_type_ == OPERATOR_LITTLE_K_DIVK_BASE) {
//       for (int n = 0; n < nfaces; n++) kf[n] = std::sqrt((*k_face)[0][faces[n]]);
//     } else if (little_k_type_ == OPERATOR_LITTLE_K_STANDARD && k_cell != Teuchos::null) {
//       kc = (*k_cell)[0][c];
//       for (int n = 0; n < nfaces; n++) kf[n] = kc;
//     } else if (little_k_type_ == OPERATOR_LITTLE_K_UPWIND) {
//       for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
//     }

//     if (K_.get()) Kc = (*K_)[c];
//     AmanziGeometry::Point Kcg(Kc * g_);

//     for (int n = 0; n < nfaces; n++) {
//       int dir, f = faces[n];
//       AmanziGeometry::Point normal = mesh_->face_normal(f, false, c, &dir);
//       normal *= dir;

//       int g = fmap.FirstPointInElement(f);
//       int ndofs = fmap.ElementSize(f);
//       if (ndofs > 1) g += Operators::UniqueIndexFaceToCells(*mesh_, f, c);

//       if (gravity_special_projection_) {
//         const AmanziGeometry::Point& xcc = Impl::GravitySpecialDirection(mesh,f);

//         double sign = normal * xcc;
//         double tmp = copysign(norm(normal) / norm(xcc), sign);
//         grav_data[0][g] += (Kcg * xcc) * rho_ * kf[n] * tmp;
//       } else {
//         grav_data[0][g] += (Kcg * normal) * rho_ * kf[n];
//       }
//     }
//   }

//   // if f is on a processor boundary, some g are not initialized
//   grav.gatherGhostedToMaster(Add);

//   for (int g = 0; g < ndofs_owned; ++g) {
//     flux_data[0][g] += grav_data[0][g];
//   }
// }


/* ******************************************************************
* Put here stuff that has to be done in constructor, i.e. only once.
****************************************************************** */
void PDE_DiffusionMFDwithGravity::Init()
{
  PDE_DiffusionMFD::Init();

  gravity_special_projection_ = (mfd_primary_ == WhetStone::DIFFUSION_TPFA);

  // gravity discretization
  std::string name = plist_.get<std::string>("gravity term discretization", "hydraulic head");
  if (name == "hydraulic head")
    gravity_method_ = OPERATOR_GRAVITY_HH;
  else
    gravity_method_ = OPERATOR_GRAVITY_FV;
}


/* ******************************************************************
* Compute non-normalized unsigned direction to the next cell needed
* to project gravity vector in the MFD-TPFA discretization method.
****************************************************************** */


/* ******************************************************************
* Return value of the gravity flux on the given face f.
****************************************************************** */
// double PDE_DiffusionMFDwithGravity::ComputeGravityFlux(int f) const
// {
//   AmanziMesh::Entity_ID_List cells;
//   mesh_->face_get_cells(f, AmanziMesh::Parallel_kind::ALL, &cells);
//   int c = cells[0];

//   double gflux;
//   const AmanziGeometry::Point& normal = mesh_->face_normal(f);

//   if (K_.get()) {
//     gflux = (((*K_)[c] * g_) * normal);
//   } else {
//     gflux = g_ * normal;
//   }

//   if (is_scalar_) {
//     gflux *= rho_;
//   } else {
//     const Epetra_MultiVector& rho_c = *rho_cv_->viewComponent("cell", true);
//     gflux *= rho_c[0][c];
//   }

//   return gflux;
// }

}  // namespace Operators
}  // namespace Amanzi
