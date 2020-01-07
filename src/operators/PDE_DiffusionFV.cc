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

// Operators
#include "Op.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"
#include "Op_Face_Cell.hh"
#include "OperatorDefs.hh"
#include "Operator_Cell.hh"
#include "PDE_DiffusionFV.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization
****************************************************************** */
void PDE_DiffusionFV::Init()
{
  // Define stencil for the FV diffusion method.
  local_op_schema_ = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL;

  // create or check the existing Operator
  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    global_op_schema_ = OPERATOR_SCHEMA_DOFS_CELL;

    // build the CVS from the global schema
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);
    cvs->AddComponent("cell", AmanziMesh::CELL, 1);

    global_op_ = Teuchos::rcp(new Operator_Cell(cvs, plist_, global_op_schema_));

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
  }

  // Do we need to exclude the primary terms?
  exclude_primary_terms_ = plist_.get<bool>("exclude primary terms", false);
  
  // create the local Op and register it with the global Operator
  if (!exclude_primary_terms_) {
    if (plist_.get<bool>("surface operator", false)) {
      std::string name = "Diffusion: FACE_CELL Surface";
      local_op_ = Teuchos::rcp(new Op_SurfaceFace_SurfaceCell(name, mesh_));
      global_op_->OpPushBack(local_op_);
    } else {
      std::string name = "Diffusion: FACE_CELL";
      local_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
      global_op_->OpPushBack(local_op_);
    }
  }
  
  // upwind options
  Errors::Message msg;
  std::string uwname = plist_.get<std::string>("nonlinear coefficient", "upwind: face");
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
  if (plist_.isParameter("newton correction")) {
    msg << "DiffusionFV: DEPRECATED: \"newton correction\" has been removed in favor of \"Newton correction\"";
    Exceptions::amanzi_throw(msg);
  }

  std::string jacobian = plist_.get<std::string>("Newton correction", "none");
  if (jacobian == "none") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;
  } else if (jacobian == "true Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_TRUE;
  } else if (jacobian == "approximate Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE;
    msg << "DiffusionFV: \"approximate Jacobian\" not supported, use \"true Jacobian\".";
    Exceptions::amanzi_throw(msg);
  } else {
    msg << "DiffusionFV: invalid parameter \"" << jacobian << "\" for option \"Newton correction\" -- valid are: \"none\", \"true Jacobian\"";
    Exceptions::amanzi_throw(msg);
  }
    

  if (newton_correction_ != OPERATOR_DIFFUSION_JACOBIAN_NONE) {
    if (plist_.get<bool>("surface operator", false)) {
      std::string name = "Diffusion: FACE_CELL Surface Jacobian terms";
      jac_op_ = Teuchos::rcp(new Op_SurfaceFace_SurfaceCell(name, mesh_));
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
void PDE_DiffusionFV::SetTensorCoefficient(const Teuchos::RCP<const std::vector<WhetStone::Tensor> >& K)
{
  transmissibility_initialized_ = false;
  K_ = K;
}


/* ******************************************************************
* Setup methods: krel and deriv -- must be called after calling a
* setup with K absolute
****************************************************************** */
void PDE_DiffusionFV::SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                           const Teuchos::RCP<const CompositeVector>& dkdp)
{
  k_ = k;
  dkdp_ = dkdp;

  if (k_ != Teuchos::null) {
    AMANZI_ASSERT(k_->HasComponent("face"));
    // NOTE: it seems that Amanzi passes in a cell based kr which is then
    // ignored, and assumed = 1.  This seems dangerous to me. --etc
    // AMANZI_ASSERT(!k_->HasComponent("cell"));
  }
  if (dkdp_ != Teuchos::null) {
    AMANZI_ASSERT(dkdp_->HasComponent("cell"));
  }
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces and 1x1 matrices
* on boundary faces.
****************************************************************** */
void PDE_DiffusionFV::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                     const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!transmissibility_initialized_) ComputeTransmissibility_();

  if (!exclude_primary_terms_) {
    auto trans_face = transmissibility_->ViewComponent<AmanziDefaultDevice>("face", true);
    auto k_face = ScalarCoefficientFaces(true);

    // updating matrix blocks
    AmanziMesh::Entity_ID_View cells;
    Kokkos::parallel_for(
        "PDE_DiffusionFV::UpdateMatrices",
        nfaces_owned,
        KOKKOS_LAMBDA(const int f) {
          mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          int ncells = cells.extent(0);

          WhetStone::DenseMatrix Aface(ncells, ncells, Kokkos::subview(local_op_->data, f, Kokkos::ALL));
          Aface = 0.0;

          double tij = trans_face(f,0) * k_face(f,0);

          for (int i = 0; i != ncells; ++i) {
            Aface(i, i) = tij;
            for (int j = i + 1; j != ncells; ++j) {
              Aface(i, j) = -tij;
              Aface(j, i) = -tij;
            }
          }
        });
  }
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces and 1x1 matrices
* on boundary faces with Newton information
****************************************************************** */
void PDE_DiffusionFV::UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    double scalar_factor)
{
  AMANZI_ASSERT(false); // not yet implemented --etc
  
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

void PDE_DiffusionFV::UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& factor)
{
  AMANZI_ASSERT(false); // not yet implemented --etc

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
void PDE_DiffusionFV::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  // prep views
  auto trans_face = transmissibility_->ViewComponent<AmanziDefaultDevice>("face", true);
  AMANZI_ASSERT(bcs_trial_.size() > 0);
  const auto bc_model = bcs_trial_[0]->bc_model();
  const auto bc_value = bcs_trial_[0]->bc_value();
  auto k_face = ScalarCoefficientFaces(false);

  if (!exclude_primary_terms_) {
    const auto rhs_cell = global_op_->rhs()->ViewComponent<AmanziDefaultDevice>("cell", true);

    AmanziMesh::Entity_ID_View cells;

    Kokkos::parallel_for(
        "PDE_DiffusionFV::ApplyBCs",
        nfaces_owned,
        KOKKOS_LAMBDA(const int f) {
          mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          int ncells = cells.extent(0);

          if (bc_model(f) == OPERATOR_BC_DIRICHLET && primary) {
            double tij = trans_face(f,0) * k_face.extent(0) ? k_face(f,0) : 1.0;
            Kokkos::atomic_add(&rhs_cell(cells(0),0), bc_value(f) * tij);
          } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
            local_op_->CopyMasterToShadow(f); // -- NOTE this is probably illegal! --etc
            local_op_->Zero(f);
          }
          
          if (primary) Kokkos::atomic_add(&rhs_cell(cells(0),0), -bc_value(f) * mesh_->face_area(f));
        });
  }

  if (jac_op_ != Teuchos::null) {
    AMANZI_ASSERT(false); // not yet implemented --etc
    // AMANZI_ASSERT(false) // 
    // AMANZI_ASSERT(bc_model.size() == nfaces_wghost);
    // for (int f = 0; f != nfaces_owned; ++f) {
    //   WhetStone::DenseMatrix& Aface = jac_op_->matrices[f];

    //   if (bc_model[f] == OPERATOR_BC_NEUMANN) {
    //     jac_op_->matrices_shadow[f] = Aface;
    //     Aface *= 0.0;
    //   }
    // }
  }
}


/* ******************************************************************
* Calculate mass flux from cell-centered data
****************************************************************** */
void PDE_DiffusionFV::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& solution,
                                 const Teuchos::Ptr<CompositeVector>& darcy_mass_flux)
{
  // prep views
  auto trans_face = transmissibility_->ViewComponent<AmanziDefaultDevice>("face", true);
  AMANZI_ASSERT(bcs_trial_.size() > 0);
  const auto bc_model = bcs_trial_[0]->bc_model();
  const auto bc_value = bcs_trial_[0]->bc_value();

  const auto k_face = ScalarCoefficientFaces(true);
  auto flux = darcy_mass_flux->ViewComponent<AmanziDefaultDevice>("face", false);

  solution->ScatterMasterToGhosted("cell");
  const auto p = solution->ViewComponent<AmanziDefaultDevice>("cell", true);

  AmanziMesh::Entity_ID_View cells, faces;
  AmanziMesh::Entity_Dir_View dirs;
  Kokkos::View<int*> flag("flags", nfaces_wghost); // initialized to 0 by default

  Kokkos::parallel_for(
      "PDE_DiffusionFV::UpdateFlux outer loop",
      ncells_owned,
      KOKKOS_LAMBDA(const int c) {
        mesh_->cell_get_faces_and_dirs(c, faces, dirs);
        int nfaces = faces.size();

        for (int n = 0; n < nfaces; n++) {
          int f = faces(n);

          if (bc_model(f) == OPERATOR_BC_DIRICHLET) {
            // guaranteed single touch?  no internal DIRICHLET BCs required or we need atomics here! --etc
            flux(f,0) = dirs(n) * trans_face(f,0) * (p(c,0) - bc_value(f)) * k_face(f,0);

          } else if (bc_model(f) == OPERATOR_BC_NEUMANN) {
            flux[0][f] = dirs(n) * bc_value(f) * mesh_->face_area(f);
        
          } else {
            // this needs more thought --etc
            if (f < nfaces_owned && Kokkos::atomic_compare_exchange(&flag(f), 0, 1) == 0) {
              mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
              // if (cells.size() <= 1) {
              //   Errors::Message msg("Flow PK: These boundary conditions are not supported by FV.");
              //   Exceptions::amanzi_throw(msg);
              // }
              assert(cells.extent(0) >= 2);
              int c1 = cells[0];
              int c2 = cells[1];
              if (c == c1) {
                flux(f,0) = dirs(n) * trans_face(f,0) * (p(c1,0) - p(c2,0)) * k_face(f,0);
              } else {
                flux(f,0) = dirs(n) * trans_face(f,0) * (p(c2,0) - p(c1,0)) * k_face(f,0);
              }            
            }
          }
        }
      });
}


/* ******************************************************************
* Calculate flux from cell-centered data
****************************************************************** */
void PDE_DiffusionFV::UpdateFluxNonManifold(const Teuchos::Ptr<const CompositeVector>& u,
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
void PDE_DiffusionFV::AnalyticJacobian_(const CompositeVector& u)
{
    AMANZI_ASSERT(false); // not yet implemented --etc
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

  // for (int f=0; f!=nfaces_owned; ++f) {
  //   mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  //   int mcells = cells.size();

  //   WhetStone::DenseMatrix Aface(mcells, mcells);

  //   for (int n = 0; n < mcells; n++) {
  //     int c1 = cells[n];
  //     pres[n] = uc[0][c1];
  //     dkdp[n] = dKdP_cell[0][c1];
  //   }

  //   if (mcells == 1) {
  //     dkdp[1] = dKdP_face.get() ? (*dKdP_face)[0][f] : 0.;
  //   }

  //   // find the face direction from cell 0 to cell 1
  //   AmanziMesh::Entity_ID_List cfaces;
  //   std::vector<int> fdirs;
  //   mesh_->cell_get_faces_and_dirs(cells[0], &cfaces, &fdirs);
  //   int f_index = std::find(cfaces.begin(), cfaces.end(), f) - cfaces.begin();
  //   ComputeJacobianLocal_(mcells, f, fdirs[f_index], bc_model[f], bc_value[f],
  //                         pres, dkdp, Aface);

  //   jac_op_->matrices[f] = Aface;
  // }
}


/* ******************************************************************
* Computation of a local submatrix of the analytical Jacobian 
* (its nonlinear part) on face f.
****************************************************************** */
// void PDE_DiffusionFV::ComputeJacobianLocal_(
//     int mcells, int f, int face_dir_0to1, int bc_model_f, double bc_value_f,
//     double *pres, double *dkdp_cell, WhetStone::DenseMatrix& Jpp)
// {
//   const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);
//   double dKrel_dp[2];
//   double dpres;

//   if (mcells == 2) {
//     dpres = pres[0] - pres[1];  // + grn;
//     if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
//       double flux0to1;
//       flux0to1 = trans_face[0][f] * dpres;
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
void PDE_DiffusionFV::ComputeTransmissibility_()
{
  transmissibility_->putScalar(0.);

  Kokkos::View<double**> K;
  int d, rank;
  if (K_.get()) {
    K = K_->data;
    d = K_->d;
    rank = K_->rank;
  } else {
    Kokkos::resize(K, ncells_owned, 1);
    Kokkos::deep_copy(K, 1.0); // again, not sure this putScalar works/exists... --etc
    d = mesh_->space_dimension();
    rank = 1;
  }

  {
    auto trans_face = transmissibility_->ViewComponent<AmanziDefaultDevice>("face", true);
    Kokkos::parallel_for(
        "PDE_DiffusionFV::ComputeTransmissibility",
        ncells_owned,
        KOKKOS_LAMBDA(const int c) {
          AmanziMesh::Entity_ID_View faces;
          Kokkos::View<AmanziGeometry::Point*> bisectors;
          mesh_->cell_get_faces_and_bisectors(c, faces, bisectors);
          int nfaces = faces.extent(0);

          // ????
          WhetStone::Tensor Kc(d, rank, Kokkos::subview(K, c, Kokkos::ALL));

          for (int i = 0; i < nfaces; i++) {
            int f = faces(i);
            const AmanziGeometry::Point& a = bisectors(i);
            const AmanziGeometry::Point& normal = mesh_->face_normal(f);
            double area = mesh_->face_area(f);

            double h_tmp = norm(a);
            double s = area / h_tmp;
            double perm = ((Kc * a) * normal) * s;
            double dxn = a * normal;
            Kokkos::atomic_add(&trans_face(f,0), fabs(dxn / perm));
          }
        });
  }
  transmissibility_->GatherGhostedToMaster();
  transmissibility_->reciprocal(*transmissibility_);
  transmissibility_->ScatterMasterToGhosted("face");
  transmissibility_initialized_ = true;
}


/* ******************************************************************
* Return transmissibility value on the given face f.
****************************************************************** */
// double PDE_DiffusionFV::ComputeTransmissibility(int f) const
// {
//   const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);
//   return trans_face[0][f];
// }

}  // namespace Operators
}  // namespace Amanzi


