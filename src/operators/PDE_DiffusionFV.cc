/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Daniil Svyatskiy (dasvyat@lanl.gov)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <vector>

// Operators
#include "Op.hh"
//#include "Op_SurfaceFace_SurfaceCell.hh"
#include "Op_Face_Cell.hh"
#include "OperatorDefs.hh"
#include "Operator_Cell.hh"
#include "PDE_DiffusionFV.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Initialization
 ****************************************************************** */
void
PDE_DiffusionFV::Init_(Teuchos::ParameterList& plist)
{
  // Define stencil for the FV diffusion method.
  local_op_schema_ = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL;

  // create or check the existing Operator
  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    global_op_schema_ = OPERATOR_SCHEMA_DOFS_CELL;

    // build the CVS from the global schema
    Teuchos::RCP<CompositeVectorSpace> cvs =
      Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);
    cvs->AddComponent("cell", AmanziMesh::CELL, 1);

    global_op_ = Teuchos::rcp(new Operator_Cell(cvs->CreateSpace(), plist, global_op_schema_));

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->getDomainMap()->Mesh();
  }

  // Do we need to exclude the primary terms?
  exclude_primary_terms_ = plist.get<bool>("exclude primary terms", false);

  // create the local Op and register it with the global Operator
  if (!exclude_primary_terms_) {
    if (plist.get<bool>("surface operator", false)) {
      assert(false && "PDE_DiffusionFV: surface operator not impelemented");
      // std::string name = "Diffusion: FACE_CELL Surface";
      // local_op_ = Teuchos::rcp(new Op_SurfaceFace_SurfaceCell(name, mesh_));
      // global_op_->OpPushBack(local_op_);
    } else {
      std::string name = "Diffusion: FACE_CELL";
      local_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
      global_op_->OpPushBack(local_op_);
    }
  }

  // upwind options
  Errors::Message msg;
  std::string uwname =
    plist.get<std::string>("nonlinear coefficient", "upwind: face");
  if (uwname == "none") {
    little_k_ = OPERATOR_LITTLE_K_NONE;

  } else if (uwname == "upwind: face") {
    little_k_ = OPERATOR_LITTLE_K_UPWIND;

  } else if (uwname == "divk: cell-face" ||
             uwname == "divk: cell-grad-face-twin" ||
             uwname == "standard: cell") {
    msg << "DiffusionFV: \"" << uwname << "\" upwinding not supported.";
    Exceptions::amanzi_throw(msg);

  } else {
    msg << "DiffusionFV: unknown upwind scheme specified.";
    Exceptions::amanzi_throw(msg);
  }

  // Do we need to calculate Newton correction terms?
  newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;

  // DEPRECATED INPUT -- remove this error eventually --etc
  if (plist.isParameter("newton correction")) {
    msg << "DiffusionFV: DEPRECATED: \"newton correction\" has been removed in "
           "favor of \"Newton correction\"";
    Exceptions::amanzi_throw(msg);
  }

  std::string jacobian = plist.get<std::string>("Newton correction", "none");
  if (jacobian == "none") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;
  } else if (jacobian == "true Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_TRUE;
  } else if (jacobian == "approximate Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE;
    msg << "DiffusionFV: \"approximate Jacobian\" not supported, use \"true "
           "Jacobian\".";
    Exceptions::amanzi_throw(msg);
  } else {
    msg << "DiffusionFV: invalid parameter \"" << jacobian
        << "\" for option \"Newton correction\" -- valid are: \"none\", \"true "
           "Jacobian\"";
    Exceptions::amanzi_throw(msg);
  }


  if (newton_correction_ != OPERATOR_DIFFUSION_JACOBIAN_NONE) {
    if (plist.get<bool>("surface operator", false)) {
      assert(false && "PDE_DiffusionFV: surface operator not impelemented");
      // std::string name = "Diffusion: FACE_CELL Surface Jacobian terms";
      // jac_op_ = Teuchos::rcp(new Op_SurfaceFace_SurfaceCell(name, mesh_));
    } else {
      std::string name = "Diffusion: FACE_CELL Jacobian terms";
      jac_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
    }

    global_op_->OpPushBack(jac_op_);
  }

  // solution-independent data
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("face", AmanziMesh::FACE, 1);
  transmissibility_ = cvs.Create();
}


/* ******************************************************************
 * Setup methods: scalar coefficients
 ****************************************************************** */


void
PDE_DiffusionFV::SetTensorCoefficient(
    const WhetStone::TensorArray& K)
{
  transmissibility_initialized_ = false;
  K_ = K;
  if (K_.extent(1) != 0) assert((K_.extent(1) == 0) && "Tensor perm not implemented.");
}


/* ******************************************************************
 * Setup methods: krel and deriv -- must be called after calling a
 * setup with K absolute
 ****************************************************************** */
void
PDE_DiffusionFV::SetScalarCoefficient(
  const Teuchos::RCP<const CompositeVector>& k,
  const Teuchos::RCP<const CompositeVector>& dkdp)
{
  transmissibility_initialized_ = false;
  k_ = k;
  dkdp_ = dkdp;

  if (k_ != Teuchos::null) {
    AMANZI_ASSERT(k_->HasComponent("face"));
    // NOTE: it seems that Amanzi passes in a cell based kr which is then
    // ignored, and assumed = 1.  This seems dangerous to me. --etc
    // AMANZI_ASSERT(!k_->HasComponent("cell"));
  }
  if (dkdp_ != Teuchos::null) { AMANZI_ASSERT(dkdp_->HasComponent("cell")); }
}


/* ******************************************************************
 * Populate face-based 2x2 matrices on interior faces and 1x1 matrices
 * on boundary faces.
 ****************************************************************** */
void
PDE_DiffusionFV::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!transmissibility_initialized_) ComputeTransmissibility_();

  if (!exclude_primary_terms_) {
    auto trans_k = transmissibility_->ViewComponent("face", false);
    if (k_ != Teuchos::null && k_->HasComponent("face")) {
      MultiVector_ptr_type trans_f = transmissibility_->GetComponent("face", false);
      MultiVector_type trans_k_vec(*trans_f);

      trans_k_vec.elementWiseMultiply(1., *trans_f->getVector(0),
              *k_->GetComponent("face", false), 0.);
      trans_k = trans_k_vec.getLocalViewDevice();
    }

    auto mesh_p = mesh_.get();
    auto data = local_matrices()->data;
    Kokkos::parallel_for(
        "PDE_DiffusionFV::UpdateMatrices",
        nfaces_owned,
        KOKKOS_LAMBDA(const int& f) {
          data(f,0) = trans_k(f,0);
          data(f,1) = -trans_k(f,0);
          data(f,2) = -trans_k(f,0);
          data(f,3) = trans_k(f,0);
        });
  }

}


/* ******************************************************************
 * Populate face-based 2x2 matrices on interior faces and 1x1 matrices
 * on boundary faces with Newton information
 ****************************************************************** */
void
PDE_DiffusionFV::UpdateMatricesNewtonCorrection(
  const Teuchos::Ptr<const CompositeVector>& flux,
  const Teuchos::Ptr<const CompositeVector>& u, double scalar_factor)
{
  assert(false && "PDE_DiffusionFV::UpdateMatricesNewtonCorrection not implemented.");
  // // Add derivatives to the matrix (Jacobian in this case)
  // if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_TRUE && u.get()) {
  //   AMANZI_ASSERT(u != Teuchos::null);

  //   if (k_ != Teuchos::null) {
  //     if (k_->HasComponent("face")) k_->ScatterMasterToGhosted("face");
  //   }
  //   if (dkdp_ != Teuchos::null) {
  //     if (dkdp_->HasComponent("face")) dkdp_->ScatterMasterToGhosted("face");
  //   }

  //   AnalyticJacobian_(*u);
  // }
}

void
PDE_DiffusionFV::UpdateMatricesNewtonCorrection(
  const Teuchos::Ptr<const CompositeVector>& flux,
  const Teuchos::Ptr<const CompositeVector>& u,
  const Teuchos::Ptr<const CompositeVector>& factor)
{
  assert(false && "PDE_DiffusionFV::UpdateMatricesNewtonCorrection not implemented.");
  // // Add derivatives to the matrix (Jacobian in this case)
  // if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_TRUE && u.get()) {
  //   AMANZI_ASSERT(u != Teuchos::null);

  //   if (k_ != Teuchos::null) {
  //     if (k_->HasComponent("face")) k_->ScatterMasterToGhosted("face");
  //   }
  //   if (dkdp_ != Teuchos::null) {
  //     if (dkdp_->HasComponent("face")) dkdp_->ScatterMasterToGhosted("face");
  //   }

  //   AnalyticJacobian_(*u);
  // }
}


/* ******************************************************************
 * Special implementation of boundary conditions.
 ****************************************************************** */
void
PDE_DiffusionFV::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  assert(false && "PDE_DiffusionFV::ApplyBCs not implemented.");
  // const Epetra_MultiVector& trans_face =
  //   *transmissibility_->ViewComponent("face", true);

  // AMANZI_ASSERT(bcs_trial_.size() > 0);
  // const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  // const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  // Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  // if (k_ != Teuchos::null) { k_face = k_->ViewComponent("face", true); }

  // if (!exclude_primary_terms_) {
  //   Epetra_MultiVector& rhs_cell =
  //     *global_op_->rhs()->ViewComponent("cell", true);

  //   AmanziMesh::Entity_ID_List cells;

  //   for (int f = 0; f < nfaces_owned; f++) {
  //     if (bc_model[f] != OPERATOR_BC_NONE) {
  //       mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  //       int c = cells[0];

  //       if (bc_model[f] == OPERATOR_BC_DIRICHLET && primary) {
  //         rhs_cell[0][c] += bc_value[f] * trans_face[0][f] *
  //                           (k_face.get() ? (*k_face)[0][f] : 1.0);
  //       } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
  //         local_op_->matrices_shadow[f] = local_op_->matrices[f];
  //         local_op_->matrices[f](0, 0) = 0.0;

  //         if (primary) rhs_cell[0][c] -= bc_value[f] * mesh_->face_area(f);
  //       }
  //     }
  //   }
  // }

  // if (jac_op_ != Teuchos::null) {
  //   AMANZI_ASSERT(bc_model.size() == nfaces_wghost);
  //   for (int f = 0; f != nfaces_owned; ++f) {
  //     WhetStone::DenseMatrix& Aface = jac_op_->matrices[f];

  //     if (bc_model[f] == OPERATOR_BC_NEUMANN) {
  //       jac_op_->matrices_shadow[f] = Aface;
  //       Aface *= 0.0;
  //     }
  //   }
  // }
}


/* ******************************************************************
 * Calculate mass flux from cell-centered data
 ****************************************************************** */
void
PDE_DiffusionFV::UpdateFlux(
  const Teuchos::Ptr<const CompositeVector>& solution,
  const Teuchos::Ptr<CompositeVector>& darcy_mass_flux)
{
  assert(false && "PDE_DiffusionFV::UpdateFlux not implemented.");
  // const Epetra_MultiVector& trans_face =
  //   *transmissibility_->ViewComponent("face", true);

  // const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  // const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  // solution->ScatterMasterToGhosted("cell");
  // if (k_ != Teuchos::null) k_->ScatterMasterToGhosted("face");

  // const Teuchos::Ptr<const Epetra_MultiVector> Krel_face =
  //   k_.get() ? k_->ViewComponent("face", false).ptr() : Teuchos::null;

  // const Epetra_MultiVector& p = *solution->ViewComponent("cell", true);
  // Epetra_MultiVector& flux = *darcy_mass_flux->ViewComponent("face", false);

  // AmanziMesh::Entity_ID_List cells, faces;
  // std::vector<int> dirs;

  // std::vector<int> flag(nfaces_wghost, 0);

  // for (int c = 0; c < ncells_owned; c++) {
  //   mesh_->cell_get_faces_and_dirs(c, &faces, dirs);
  //   int nfaces = faces.size();

  //   for (int n = 0; n < nfaces; n++) {
  //     int f = faces[n];

  //     if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
  //       double value = bc_value[f];
  //       flux[0][f] = dirs[n] * trans_face[0][f] * (p[0][c] - value);
  //       if (Krel_face.get()) flux[0][f] *= (*Krel_face)[0][f];

  //     } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
  //       double value = bc_value[f];
  //       double area = mesh_->face_area(f);
  //       flux[0][f] = dirs[n] * value * area;

  //     } else {
  //       if (f < nfaces_owned && !flag[f]) {
  //         mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  //         if (cells.size() <= 1) {
  //           Errors::Message msg(
  //             "Flow PK: These boundary conditions are not supported by FV.");
  //           Exceptions::amanzi_throw(msg);
  //         }
  //         int c1 = cells[0];
  //         int c2 = cells[1];
  //         if (c == c1) {
  //           flux[0][f] = dirs[n] * trans_face[0][f] * (p[0][c1] - p[0][c2]);
  //         } else {
  //           flux[0][f] = dirs[n] * trans_face[0][f] * (p[0][c2] - p[0][c1]);
  //         }
  //         if (Krel_face.get()) flux[0][f] *= (*Krel_face)[0][f];
  //         flag[f] = 1;
  //       }
  //     }
  //   }
  // }
}


/* ******************************************************************
 * Calculate flux from cell-centered data
 ****************************************************************** */
void
PDE_DiffusionFV::UpdateFluxNonManifold(
  const Teuchos::Ptr<const CompositeVector>& u,
  const Teuchos::Ptr<CompositeVector>& flux)
{
  Errors::Message msg;
  msg << "DiffusionFV: missing support for non-manifolds.";
  Exceptions::amanzi_throw(msg);
}


/* ******************************************************************
 * Computation the part of the Jacobian which depends on derivatives
 * of the relative permeability wrt to capillary pressure. They must
 * be added to the existing matrix structure.
 ****************************************************************** */
void
PDE_DiffusionFV::AnalyticJacobian_(const CompositeVector& u)
{
  assert(false && "PDE_DiffusionFV::AnalyticJacobian not implemented.");
  // const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  // const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  // u.ScatterMasterToGhosted("cell");
  // const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);

  // const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  // AmanziMesh::Entity_ID_List cells, faces;

  // double k_rel[2], dkdp[2], pres[2], dist;

  // const Epetra_MultiVector& dKdP_cell = *dkdp_->ViewComponent("cell");
  // Teuchos::RCP<const Epetra_MultiVector> dKdP_face;
  // if (dkdp_->HasComponent("face")) {
  //   dKdP_face = dkdp_->ViewComponent("face", true);
  // }

  // for (int f = 0; f != nfaces_owned; ++f) {
  //   mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  //   int mcells = cells.size();

  //   WhetStone::DenseMatrix Aface(mcells, mcells);

  //   for (int n = 0; n < mcells; n++) {
  //     int c1 = cells[n];
  //     pres[n] = uc[0][c1];
  //     dkdp[n] = dKdP_cell[0][c1];
  //   }

  //   if (mcells == 1) { dkdp[1] = dKdP_face.get() ? (*dKdP_face)[0][f] : 0.; }

  //   // find the face direction from cell 0 to cell 1
  //   AmanziMesh::Entity_ID_List cfaces;
  //   std::vector<int> fdirs;
  //   mesh_->cell_get_faces_and_dirs(cells[0], &cfaces, fdirs);
  //   int f_index = std::find(cfaces.begin(), cfaces.end(), f) - cfaces.begin();
  //   ComputeJacobianLocal_(
  //     mcells, f, fdirs[f_index], bc_model[f], bc_value[f], pres, dkdp, Aface);

  //   jac_op_->matrices[f] = Aface;
  // }
}


// /* ******************************************************************
//  * Computation of a local submatrix of the analytical Jacobian
//  * (its nonlinear part) on face f.
//  ****************************************************************** */
// void
// PDE_DiffusionFV::ComputeJacobianLocal_(int mcells, int f, int face_dir_0to1,
//                                        int bc_model_f, double bc_value_f,
//                                        double* pres, double* dkdp_cell,
//                                        WhetStone::DenseMatrix& Jpp)
// {
//   const Epetra_MultiVector& trans_face =
//     *transmissibility_->ViewComponent("face", true);
//   double dKrel_dp[2];
//   double dpres;

//   if (mcells == 2) {
//     dpres = pres[0] - pres[1]; // + grn;
//     if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
//       double flux0to1;
//       flux0to1 = trans_face[0][f] * dpres;
//       if (flux0to1 > OPERATOR_UPWIND_RELATIVE_TOLERANCE) { // Upwind
//         dKrel_dp[0] = dkdp_cell[0];
//         dKrel_dp[1] = 0.0;
//       } else if (flux0to1 < -OPERATOR_UPWIND_RELATIVE_TOLERANCE) { // Upwind
//         dKrel_dp[0] = 0.0;
//         dKrel_dp[1] = dkdp_cell[1];
//       } else if (fabs(flux0to1) <
//                  OPERATOR_UPWIND_RELATIVE_TOLERANCE) { // Upwind
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

//     Jpp(1, 0) = -Jpp(0, 0);
//     Jpp(1, 1) = -Jpp(0, 1);

//   } else if (mcells == 1) {
//     if (bc_model_f == OPERATOR_BC_DIRICHLET) {
//       pres[1] = bc_value_f;
//       dpres = pres[0] - pres[1];
//       Jpp(0, 0) = trans_face[0][f] * dpres * dkdp_cell[0];
//     } else {
//       Jpp(0, 0) = 0.0;
//     }
//   }
// }


/* ******************************************************************
 * Compute transmissibilities on faces
 ****************************************************************** */
void
PDE_DiffusionFV::ComputeTransmissibility_()
{
  transmissibility_->putScalarMasterAndGhosted(0.);
  auto mesh_p = mesh_.get();

  if (K_.extent(1) != 0) assert((K_.extent(1) == 0) && "Tensor perm not implemented.");

  {
    auto trans_f = transmissibility_->ViewComponent("face", true);
  
    for(int c = 0 ; c < ncells_owned; ++c){
    //Kokkos::parallel_for(
    //    "PDE_DiffusionFV::ComputeTransmissibility",
    //    ncells_owned,
    //    KOKKOS_LAMBDA(const int& c) {
      AmanziMesh::Entity_ID_View faces;
      Kokkos::View<AmanziGeometry::Point*> bisectors;
      mesh_p->cell_get_faces_and_bisectors(c, faces, bisectors);

      WhetStone::Tensor Kc = K_.getTensor(c);
      
      Kokkos::parallel_for(
        faces.extent(0),
        KOKKOS_LAMBDA (const int& i) {
          auto f = faces(i);
          const AmanziGeometry::Point& a = bisectors(i);
          const AmanziGeometry::Point& normal = mesh_p->face_normal(f);

          double s = mesh_p->face_area(f) / AmanziGeometry::norm(a);

          double perm = ((Kc * a) * normal) * s;
          double dxn = a * normal;
          Kokkos::atomic_add(&trans_f(f,0), fabs(dxn / perm));
      });
    }
        //});
  }
  transmissibility_->GatherGhostedToMaster(Tpetra::ADD);
  transmissibility_->reciprocal(*transmissibility_);
  transmissibility_->ScatterMasterToGhosted("face");
  transmissibility_initialized_ = true;
}


// /* ******************************************************************
//  * Return transmissibility value on the given face f.
//  ****************************************************************** */
// double
// PDE_DiffusionFV::ComputeTransmissibility(int f) const
// {
//   const Epetra_MultiVector& trans_face =
//     *transmissibility_->ViewComponent("face", true);
//   return trans_face[0][f];
// }

} // namespace Operators
} // namespace Amanzi
