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
#include <boost/math/tools/roots.hpp>

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
void PDE_DiffusionFVwithGravity::Init()
{
  PDE_DiffusionFV::Init();
  gravity_term_ = Teuchos::rcp(new CompositeVector(*transmissibility_));
}


/* ******************************************************************
* This completes initialization of the base class.
****************************************************************** */
void PDE_DiffusionFVwithGravity::SetDensity(
    const Teuchos::RCP<const CompositeVector>& rho)
{
  PDE_Diffusion::SetDensity(rho);
  transmissibility_initialized_ = false;
}


void PDE_DiffusionFVwithGravity::SetDensity(double rho) {
  PDE_Diffusion::SetDensity(rho);
  transmissibility_initialized_ = false;
}

  
/* ******************************************************************
* Populate face-based matrices.
****************************************************************** */
void PDE_DiffusionFVwithGravity::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!transmissibility_initialized_) ComputeTransmissibility_(gravity_term_);

  PDE_DiffusionFV::UpdateMatrices(flux, u);

  // populating right-hand side
  if (local_op_.get()) {
    const auto bc_model = bcs_trial_[0]->bc_model();
    const auto grav_f = gravity_term_->ViewComponent("face", true);
    auto rhs_c = global_op_->rhs()->ViewComponent("cell");
    const auto k_face = ScalarCoefficientFaces(true);
    const AmanziMesh::Mesh* m = mesh_.get(); 

    Kokkos::parallel_for(
        "PDE_DiffusionFVwithGravity::UpdateMatrices",
        ncells_owned,
        KOKKOS_LAMBDA(const int c) {
          AmanziMesh::Entity_ID_View faces;
          AmanziMesh::Entity_Dir_View dirs;
          m->cell_get_faces_and_dirs(c, faces, dirs);
          int nfaces = faces.size();

          for (int n = 0; n != nfaces; ++n) {
            int f = faces[n];
            if (bc_model(f) != OPERATOR_BC_NEUMANN)
              Kokkos::atomic_add(&rhs_c(c,0), -dirs(n) * grav_f(f,0) * k_face(f,0));
          }
        });
  }
}


/* ******************************************************************
* Special implementation of boundary conditions.
****************************************************************** */
void PDE_DiffusionFVwithGravity::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  PDE_DiffusionFV::ApplyBCs(primary, eliminate, essential_eqn);


  auto grav_f = gravity_term_->ViewComponent("face", true);
  const auto bc_model = bcs_trial_[0]->bc_model();
  Kokkos::parallel_for(
      "PDE_DiffusionFVwithGravity::ApplyBCs",
      nfaces_owned,
      KOKKOS_LAMBDA(const int f) {
        if (bc_model(f) == OPERATOR_BC_NEUMANN) grav_f(f,0) = 0.;
      });
}


/* ******************************************************************
* Calculate flux from cell-centered data
****************************************************************** */
void PDE_DiffusionFVwithGravity::UpdateFlux(
    const Teuchos::Ptr<const CompositeVector>& solution,
    const Teuchos::Ptr<CompositeVector>& darcy_mass_flux)
{
  PDE_DiffusionFV::UpdateFlux(solution, darcy_mass_flux);

  // add gravity term
  const auto& flux = darcy_mass_flux->GetComponent("face", false);
  const auto& grav_f = gravity_term_->GetComponent("face", false);
  if (k_.get()) {
    const auto& k_face = k_->GetComponent("face", false);
    flux->elementWiseMultiply(1.0, *k_face->getVector(0), *grav_f, 1.0);
  } else {
    flux->update(1.0, *grav_f, 1.0);
  }
}


// /* ******************************************************************
// * Computation of a local submatrix of the analytical Jacobian 
// * (its nonlinear part) on face f.
// ****************************************************************** */
// void PDE_DiffusionFVwithGravity::ComputeJacobianLocal_(
//     int mcells, int f, int face_dir_0to1, int bc_model_f, double bc_value_f,
//     double *pres, double *dkdp_cell, WhetStone::DenseMatrix& Jpp)
// {
//   const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);
//   const Teuchos::Ptr<Epetra_MultiVector> gravity_face =
//       gravity_term_.get() ? gravity_term_->ViewComponent("face", true).ptr() : Teuchos::null;

//   double dKrel_dp[2];
//   double dpres;

//   if (mcells == 2) {
//     dpres = pres[0] - pres[1];  // + grn;
//     if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
//       double flux0to1;
//       flux0to1 = trans_face[0][f] * dpres;
//       if (gravity_face.get()) flux0to1 += face_dir_0to1 * (*gravity_face)[0][f];
//       if (flux0to1  > OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
//         dKrel_dp[0] = dkdp_cell[0];
//         dKrel_dp[1] = 0.0;
//       } else if (flux0to1 < -OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
//         dKrel_dp[0] = 0.0;
//         dKrel_dp[1] = dkdp_cell[1];
//       } else if (fabs(flux0to1) < OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
//         dKrel_dp[0] = 0.5 * dkdp_cell[0];
//         dKrel_dp[1] = 0.5 * dkdp_cell[1];
//       }
//     } else if (little_k_ == OPERATOR_UPWIND_ARITHMETIC_AVERAGE) {
//       dKrel_dp[0] = 0.5 * dkdp_cell[0];
//       dKrel_dp[1] = 0.5 * dkdp_cell[1];
//     } else {
//       AMANZI_ASSERT(0);
//     }

//     Jpp(0, 0) = trans_face[0][f] * dpres * dKrel_dp[0];
//     Jpp(0, 1) = trans_face[0][f] * dpres * dKrel_dp[1];

//     if (gravity_face.get()) {
//       Jpp(0,0) += face_dir_0to1 * (*gravity_face)[0][f] * dKrel_dp[0];
//       Jpp(0,1) += face_dir_0to1 * (*gravity_face)[0][f] * dKrel_dp[1];
//     }
//     Jpp(1, 0) = -Jpp(0, 0);
//     Jpp(1, 1) = -Jpp(0, 1);

//   } else if (mcells == 1) {
//     if (bc_model_f == OPERATOR_BC_DIRICHLET) {                   
//       pres[1] = bc_value_f;
//       dpres = pres[0] - pres[1];  // + grn;
//       Jpp(0, 0) = trans_face[0][f] * dpres * dkdp_cell[0];
//       if (gravity_face.get())
//           Jpp(0,0) += face_dir_0to1 * (*gravity_face)[0][f] * dkdp_cell[0];
//     } else {
//       Jpp(0, 0) = 0.0;
//     }
//   }
// }


/* ******************************************************************
* Compute transmissibilities on faces. Requires K, g, rho.
****************************************************************** */
void PDE_DiffusionFVwithGravity::ComputeTransmissibility_(
   Teuchos::RCP<CompositeVector> g_cv)
{
  transmissibility_->putScalar(0.);

  // Compute auxiliary structure. Note that first components of both 
  // fields are used symmetrically (no specific order).
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("face", AmanziMesh::FACE, 1);
  auto h = cvs.Create();

  {

    if (!K_.get()) {
      CompositeVectorSpace cvs;
      cvs.SetMesh(mesh_)
          ->SetGhosted(false)
          ->AddComponent("cell", AmanziMesh::CELL, 1);
      auto K = Teuchos::rcp(new TensorVector(cvs,false));
      for (int c=0; c!=ncells_owned; ++c) {
        K->set_shape(c, mesh_->space_dimension(), 1);
      }
      K->Init();
      K_ = K;
    }
    const TensorVector* K = K_.get();    
  
    auto beta_f = transmissibility_->ViewComponent("face", true);
    auto h_f = h->ViewComponent("face", true);
    const AmanziMesh::Mesh* m = mesh_.get(); 
    Kokkos::parallel_for(
        "PDE_DiffusionFVwithGravity::ComputeTransmissibility1",
        ncells_owned,
        KOKKOS_LAMBDA(const int c) {

          AmanziMesh::Entity_ID_View faces;
          Kokkos::View<AmanziGeometry::Point*> bisectors;
          m->cell_get_faces_and_bisectors(c, faces, bisectors);

          WhetStone::Tensor<DeviceOnlyMemorySpace> Kc = K->at(c);
          //if (c == 0) std::cout << "Kc = " << Kc(0,0) << std::endl;

          for (int i = 0; i < faces.extent(0); i++) {
            auto f = faces(i);
            const AmanziGeometry::Point& a = bisectors(i);
            const AmanziGeometry::Point& normal = m->face_normal(f);
            const double area = m->face_area(f);

            const double h_tmp = AmanziGeometry::norm(a);
            const double s = area / h_tmp;
            const double perm = ((Kc * a) * normal) * s;
            const double dxn = a * normal;

            Kokkos::atomic_add(&h_f(f,0), h_tmp);
            //if (c == 0) std::cout << "  1/trans_f(" << f << ") += " << dxn << "/" << perm << " = " << fabs(dxn/perm) << std::endl;
            Kokkos::atomic_add(&beta_f(f,0), fabs(dxn / perm));
          }
        });
  }

  transmissibility_->GatherGhostedToMaster();
  transmissibility_->reciprocal(*transmissibility_);
  h->GatherGhostedToMaster();

  {
    auto rho_c = DensityCells(true);
    auto h_f = h->ViewComponent("face", false);
    const auto trans_f = transmissibility_->ViewComponent("face", false);
    auto grav_f = g_cv->ViewComponent("face", false);
    const AmanziMesh::Mesh* m = mesh_.get(); 

    Kokkos::parallel_for(
        "PDE_DiffusionFVwithGravity::ComputeTransmissibility2",
        nfaces_owned,
        KOKKOS_LAMBDA(const int f) {
          AmanziMesh::Entity_ID_View cells;
          m->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          int ncells = cells.size();

          //if (f == 0) std::cout << "rho = " << rho_c(0,0) << std::endl;
          //if (f == 0) std::cout << "g = " << g_ << std::endl;
          AmanziGeometry::Point a_dist;
          if (ncells == 2) {
            a_dist = m->cell_centroid(cells[1]) - m->cell_centroid(cells[0]);
          } else if (ncells == 1) {
            a_dist = m->face_centroid(f) - m->cell_centroid(cells[0]);
          } 
          a_dist *= 1.0 / AmanziGeometry::norm(a_dist);

          const AmanziGeometry::Point& normal = m->face_normal(f);
          double dir = std::copysign(1.0, normal * a_dist);

          double rho = ncells == 1 ? rho_c(cells(0),0) : (rho_c(cells(0),0) + rho_c(cells(1),0))/2.;
          //if (f == 0) std::cout << "grav_f = " << trans_f(f,0) << " * " << h_f(f,0) << " * g * " << a_dist << " * " << rho << std::endl;
          grav_f(f,0) = trans_f(f,0) * h_f(f,0) * (g_ * a_dist) * rho * dir;
        });

  }

  transmissibility_->ScatterMasterToGhosted("face");
  g_cv->ScatterMasterToGhosted("face");
  transmissibility_initialized_ = true;
}


}  // namespace Operators
}  // namespace Amanzi


