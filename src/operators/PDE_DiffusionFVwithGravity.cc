/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Op.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"
#include "Op_Face_Cell.hh"
#include "OperatorDefs.hh"
#include "Operator_Cell.hh"
#include "PDE_DiffusionFVwithGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* This completes initialization of the base class.
****************************************************************** */
void
PDE_DiffusionFVwithGravity::Init()
{
  PDE_DiffusionFV::Init();
  gravity_term_ = Teuchos::rcp(new CompositeVector(*transmissibility_));
}


/* ******************************************************************
* This completes initialization of the base class.
****************************************************************** */
void
PDE_DiffusionFVwithGravity::SetDensity(const Teuchos::RCP<const CompositeVector>& rho)
{
  PDE_DiffusionWithGravity::SetDensity(rho);
  transmissibility_initialized_ = false;
}


void
PDE_DiffusionFVwithGravity::SetDensity(double rho)
{
  PDE_DiffusionWithGravity::SetDensity(rho);
  transmissibility_initialized_ = false;
}


/* ******************************************************************
* Populate face-based matrices.
****************************************************************** */
void
PDE_DiffusionFVwithGravity::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                           const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!transmissibility_initialized_) ComputeTransmissibility_(gravity_term_);

  PDE_DiffusionFV::UpdateMatrices(flux, u);

  // populating right-hand side
  if (local_op_.get()) {
    const auto bc_model = bcs_trial_[0]->bc_model();
    const auto grav_f = gravity_term_->viewComponent("face", true);
    auto rhs_c = global_op_->rhs()->viewComponent("cell");
    const auto k_face = ScalarCoefficientFaces(true);
    const AmanziMesh::MeshCache& m = mesh_->getCache();

    Kokkos::parallel_for(
      "PDE_DiffusionFVwithGravity::UpdateMatrices", ncells_owned, KOKKOS_LAMBDA(const int c) {
        auto [faces, dirs] = m.getCellFacesAndDirections(c);
        int nfaces = faces.size();

        for (int n = 0; n != nfaces; ++n) {
          int f = faces[n];
          if (bc_model(f) != OPERATOR_BC_NEUMANN)
            Kokkos::atomic_add(&rhs_c(c, 0), -dirs(n) * grav_f(f, 0) * k_face(f, 0));
        }
      });
  }
}


/* ******************************************************************
* Special implementation of boundary conditions.
****************************************************************** */
void
PDE_DiffusionFVwithGravity::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  PDE_DiffusionFV::ApplyBCs(primary, eliminate, essential_eqn);


  auto grav_f = gravity_term_->viewComponent("face", true);
  const auto bc_model = bcs_trial_[0]->bc_model();
  Kokkos::parallel_for(
    "PDE_DiffusionFVwithGravity::ApplyBCs", nfaces_owned, KOKKOS_LAMBDA(const int f) {
      if (bc_model(f) == OPERATOR_BC_NEUMANN) grav_f(f, 0) = 0.;
    });
}


/* ******************************************************************
* Calculate flux from cell-centered data
****************************************************************** */
void
PDE_DiffusionFVwithGravity::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& solution,
                                       const Teuchos::Ptr<CompositeVector>& darcy_mass_flux)
{
  PDE_DiffusionFV::UpdateFlux(solution, darcy_mass_flux);

  // add gravity term
  const auto& flux = darcy_mass_flux->getComponent("face", false);
  const auto& grav_f = gravity_term_->getComponent("face", false);
  if (k_.get()) {
    const auto& k_face = k_->getComponent("face", false);
    flux->elementWiseMultiply(1.0, *k_face->getVector(0), *grav_f, 1.0);
  } else {
    flux->update(1.0, *grav_f, 1.0);
  }
}


/* ******************************************************************
* Computation the part of the Jacobian which depends on derivatives
* of the relative permeability wrt to capillary pressure. They must
* be added to the existing matrix structure.
****************************************************************** */
void
PDE_DiffusionFVwithGravity::AnalyticJacobian_(const CompositeVector& u)
{
  u.scatterMasterToGhosted("cell");
  AMANZI_ASSERT(dkdp_.get() && dkdp_->hasComponent("face"));
  dkdp_->scatterMasterToGhosted("face");

  { // scope for views
    const auto bc_model = bcs_trial_[0]->bc_model();
    const auto bc_value = bcs_trial_[0]->bc_value();

    const auto uc = u.viewComponent("cell", true);
    const auto cmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);
    const auto trans_face = transmissibility_->viewComponent("face", true);
    const auto gravity_face = gravity_term_->viewComponent("face", true);

    // dkface_dp(f, 0/1): the 2-DoF, face-based derivative of the upwinded
    // coefficient wrt the pressure of cells[0]/cells[1] of that face,
    // computed by whichever upwinding evaluator was selected.
    const auto dkface_dp = dkdp_->viewComponent("face", true);

    const AmanziMesh::MeshCache& mesh = mesh_->getCache();
    DenseMatrix_Vector& A = jac_op_->A;

    Kokkos::parallel_for(
      "PDE_DiffusionFVwithGravity::AnalyticJacobian_", nfaces_owned, KOKKOS_LAMBDA(const int& f) {
        auto cells = mesh.getFaceCells(f);
        int mcells = cells.size();

        auto Aface = A[f];

        double pres[2] = {0., 0.};
        for (int n = 0; n < mcells; n++) { pres[n] = uc(cells(n), 0); }

        // find the face direction from cell 0 to cell 1
        auto [cfaces, fdirs] = mesh.getCellFacesAndDirections(cells(0));
        int f_index;
        for (f_index = 0; f_index != cfaces.extent(0); ++f_index) {
          if (cfaces(f_index) == f) break;
        }

        double dpres;

        if (mcells == 2) {
          dpres = pres[0] - pres[1];

          Aface(0, 0) = trans_face(f, 0) * dpres * dkface_dp(f, 0) +
                        fdirs(f_index) * gravity_face(f, 0) * dkface_dp(f, 0);
          Aface(0, 1) = trans_face(f, 0) * dpres * dkface_dp(f, 1) +
                        fdirs(f_index) * gravity_face(f, 0) * dkface_dp(f, 1);

          Aface(1, 0) = -Aface(0, 0);
          Aface(1, 1) = -Aface(0, 1);

        } else if (mcells == 1) {
          if (bc_model(f) == OPERATOR_BC_DIRICHLET) {
            pres[1] = bc_value(f);
            dpres = pres[0] - pres[1];
            Aface(0, 0) = trans_face(f, 0) * dpres * dkface_dp(f, 0) +
                          fdirs(f_index) * gravity_face(f, 0) * dkface_dp(f, 0);
          } else {
            Aface(0, 0) = 0.0;
          }
        }
      });
  }
}


/* ******************************************************************
* Compute transmissibilities on faces. Requires K, g, rho.
****************************************************************** */
void
PDE_DiffusionFVwithGravity::ComputeTransmissibility_(Teuchos::RCP<CompositeVector> g_cv)
{
  transmissibility_->putScalar(0.);

  // Compute auxiliary structure. Note that first components of both
  // fields are used symmetrically (no specific order).
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  auto h = cvs.Create();

  {
    if (!K_.get()) {
      CompositeVectorSpace cvs;
      cvs.SetMesh(mesh_)->SetGhosted(false)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      auto K = Teuchos::rcp(new TensorVector(cvs, false));
      for (int c = 0; c != ncells_owned; ++c) { K->set_shape(c, mesh_->getSpaceDimension(), 1); }
      K->Init();
      K_ = K;
    }
    const TensorVector& K = *K_.get();

    auto beta_f = transmissibility_->viewComponent("face", true);
    auto h_f = h->viewComponent("face", true);
    const AmanziMesh::MeshCache& m = mesh_->getCache();
    Kokkos::parallel_for(
      "PDE_DiffusionFVwithGravity::ComputeTransmissibility1",
      ncells_owned,
      KOKKOS_LAMBDA(const int c) {
        auto [faces, bisectors] = m.getCellFacesAndBisectors(c);
        auto Kc = K.at(c);

        for (int i = 0; i < faces.extent(0); i++) {
          auto f = faces(i);
          const AmanziGeometry::Point& a = bisectors(i);
          const AmanziGeometry::Point& normal = m.getFaceNormal(f);
          const double area = m.getFaceArea(f);

          const double h_tmp = AmanziGeometry::norm(a);
          const double s = area / h_tmp;
          const double perm = ((Kc * a) * normal) * s;
          const double dxn = a * normal;

          Kokkos::atomic_add(&h_f(f, 0), h_tmp);
          Kokkos::atomic_add(&beta_f(f, 0), fabs(dxn / perm));
        }
      });
  }

  transmissibility_->gatherGhostedToMaster();
  h->gatherGhostedToMaster();

  {
    auto rho_c = DensityCells(true);
    auto h_f = h->viewComponent("face", false);
    const auto trans_f = transmissibility_->viewComponent("face", false);
    auto grav_f = g_cv->viewComponent("face", false);
    const AmanziMesh::MeshCache& m = mesh_->getCache();
    AmanziGeometry::Point g(g_);

    Kokkos::parallel_for(
      "PDE_DiffusionFVwithGravity::ComputeTransmissibility2",
      nfaces_owned,
      KOKKOS_LAMBDA(const int f) {
        auto cells = m.getFaceCells(f);
        int ncells = cells.size();

        AmanziGeometry::Point a_dist;
        if (ncells == 2) {
          a_dist = m.getCellCentroid(cells[1]) - m.getCellCentroid(cells[0]);
        } else if (ncells == 1) {
          a_dist = m.getFaceCentroid(f) - m.getCellCentroid(cells[0]);
        }
        a_dist *= 1.0 / AmanziGeometry::norm(a_dist);

        trans_f(f, 0) = 1.0 / trans_f(f, 0);

        auto normal = m.getFaceNormal(f);
        double dir = std::copysign(1.0, normal * a_dist);

        double rho =
          ncells == 1 ? rho_c(cells(0), 0) : (rho_c(cells(0), 0) + rho_c(cells(1), 0)) / 2.;
        double grav = (g * a_dist) * rho * dir;
        grav *= h_f(f, 0);

        grav_f(f, 0) = trans_f(f, 0) * grav;
      });
  }

  transmissibility_->scatterMasterToGhosted("face");
  g_cv->scatterMasterToGhosted("face");
  transmissibility_initialized_ = true;
}



} // namespace Operators
} // namespace Amanzi
