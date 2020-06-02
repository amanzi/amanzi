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
//void PDE_DiffusionFV::Init_(Teuchos::ParameterList& plist)
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

    global_op_ = Teuchos::rcp(new Operator_Cell(cvs->CreateSpace(), plist_, global_op_schema_));

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->getDomainMap()->Mesh();
  }

  // Do we need to exclude the primary terms?
  bool exclude_primary_terms = plist_.get<bool>("exclude primary terms", false);
  
  // create the local Op and register it with the global Operator
  if (!exclude_primary_terms) {
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
    little_k_type_ = OPERATOR_LITTLE_K_NONE;

  } else if (uwname == "upwind: face") {
    little_k_type_ = OPERATOR_LITTLE_K_UPWIND;

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
void PDE_DiffusionFV::SetTensorCoefficient(const Teuchos::RCP<const TensorVector>& K)
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
  transmissibility_initialized_ = false;
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

  // verify that mass matrices were initialized.
  // -- this shouldn't be called here, as Trans has no dependence on
  //    rel perm, and abs perm may not be set yet! --etc
  // if (!transmissibility_initialized_)  ComputeTransmissibility_();
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces and 1x1 matrices
* on boundary faces.
****************************************************************** */
void PDE_DiffusionFV::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                     const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!transmissibility_initialized_) ComputeTransmissibility_();

  if (local_op_.get()) {
    using memory_space = decltype(local_op_->A)::memory_space;
    using DenseMatrix = WhetStone::DenseMatrix<memory_space>;
    
    const auto trans_face = transmissibility_->ViewComponent("face", true);

    // preparing upwind data
    const auto k_face = ScalarCoefficientFaces(true);
    const Amanzi::AmanziMesh::Mesh* m = mesh_.get();
    DenseMatrix_Vector& A = local_op_->A; 

    // updating matrix blocks
    Kokkos::parallel_for(
        "PDE_DiffusionFV::UpdateMatrices",
        nfaces_owned,
        KOKKOS_LAMBDA(const int f) {
          DenseMatrix A_f = A[f]; 
          double tij = trans_face(f,0) * k_face(f,0);

          for (int i = 0; i != A_f.NumRows(); ++i) {
            A_f(i, i) = tij;
            for (int j = i + 1; j != A_f.NumRows(); ++j) {
              A_f(i, j) = -tij;
              A_f(j, i) = -tij;
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
  // Add derivatives to the matrix (Jacobian in this case)
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_TRUE && u.get()) {
    AMANZI_ASSERT(u != Teuchos::null);
    AnalyticJacobian_(*u);
  }
}

void PDE_DiffusionFV::UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& factor)
{
  assert(false); 
  #if 0 
  // Add derivatives to the matrix (Jacobian in this case)
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_TRUE && u.get()) {
    AMANZI_ASSERT(u != Teuchos::null);

    if (k_ != Teuchos::null) {
      if (k_->HasComponent("face")) k_->ScatterMasterToGhosted("face");
    }
    if (dkdp_ != Teuchos::null) {
      if (dkdp_->HasComponent("face")) dkdp_->ScatterMasterToGhosted("face");
    }

    AnalyticJacobian_(*u);
  }
  #endif 
}  

/* ******************************************************************
* Special implementation of boundary conditions.
****************************************************************** */
void PDE_DiffusionFV::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  AMANZI_ASSERT(bcs_trial_.size() > 0);
  const auto bc_model = bcs_trial_[0]->bc_model();
  const auto bc_value = bcs_trial_[0]->bc_value();
  
  if (local_op_.get()) {
    auto local_op = local_op_.get(); 
    // prep views
    auto k_face = ScalarCoefficientFaces(false);
    const auto rhs_cell = global_op_->rhs()->ViewComponent<DefaultDevice>("cell", true);
    auto trans_face = transmissibility_->ViewComponent<DefaultDevice>("face", true);
    const Amanzi::AmanziMesh::Mesh* m = mesh_.get(); 
    Kokkos::parallel_for(
        "PDE_DiffusionFV::ApplyBCs",
        nfaces_owned,
        KOKKOS_LAMBDA(const int f) {
          AmanziMesh::Entity_ID_View cells;
          m->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          int ncells = cells.extent(0);

          if (bc_model(f) == OPERATOR_BC_DIRICHLET && primary) {
            double tij = trans_face(f,0) * k_face(f,0);
            Kokkos::atomic_add(&rhs_cell(cells(0),0), bc_value(f) * tij);
          } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
            local_op->Zero(f);
          
            if (primary) Kokkos::atomic_add(&rhs_cell(cells(0),0), -bc_value(f) * m->face_area(f));
          }
        });
  }

  if (jac_op_ != Teuchos::null) {
    auto& local_op = jac_op_->A;
    
    AMANZI_ASSERT(bc_model.size() == nfaces_wghost);
    Kokkos::parallel_for(
        "PDE_DiffusionFV::ApplyBCs",
        nfaces_owned,
        KOKKOS_LAMBDA(const int f) {
          if (bc_model[f] == OPERATOR_BC_NEUMANN) {
            auto A_f = local_op[f];
            A_f.putScalar(0.);
          }
        });
  }
}


/* ******************************************************************
* Calculate mass flux from cell-centered data
****************************************************************** */
void PDE_DiffusionFV::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& solution,
                                 const Teuchos::Ptr<CompositeVector>& darcy_mass_flux)
{
  // prep views
  auto trans_face = transmissibility_->ViewComponent<DefaultDevice>("face", true);
  AMANZI_ASSERT(bcs_trial_.size() > 0);
  const auto bc_model = bcs_trial_[0]->bc_model();
  const auto bc_value = bcs_trial_[0]->bc_value();

  const auto k_face = ScalarCoefficientFaces(true);
  auto flux = darcy_mass_flux->ViewComponent<DefaultDevice>("face", false);

  solution->ScatterMasterToGhosted("cell");
  const auto p = solution->ViewComponent<DefaultDevice>("cell", true);
  const Amanzi::AmanziMesh::Mesh* m = mesh_.get(); 

  Kokkos::View<int*> flag("flags", nfaces_wghost); // initialized to 0 by default
  Kokkos::parallel_for(
      "PDE_DiffusionFV::UpdateFlux outer loop",
      ncells_owned,
      KOKKOS_LAMBDA(const int c) {
        AmanziMesh::Entity_ID_View cells, faces;
        AmanziMesh::Entity_Dir_View dirs;
        m->cell_get_faces_and_dirs(c, faces, dirs);
        int nfaces = faces.size();

        for (int n = 0; n < nfaces; n++) {
          int f = faces(n);

          if (bc_model(f) == OPERATOR_BC_DIRICHLET) {
            // guaranteed single touch?  no internal DIRICHLET BCs required or we need atomics here! --etc
            flux(f,0) = dirs(n) * trans_face(f,0) * (p(c,0) - bc_value(f)) * k_face(f,0);

          } else if (bc_model(f) == OPERATOR_BC_NEUMANN) {
            // guaranteed single touch?  no internal NEUMANN BCs required or we need atomics here! --etc
            flux(f,0) = dirs(n) * bc_value(f) * m->face_area(f);
        
          } else {
            // this needs more thought --etc
            if (f < nfaces_owned && Kokkos::atomic_compare_exchange(&flag(f), 0, 1) == 0) {
              m->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
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
* Computation the part of the Jacobian which depends on derivatives 
* of the relative permeability wrt to capillary pressure. They must
* be added to the existing matrix structure.
****************************************************************** */
void PDE_DiffusionFV::AnalyticJacobian_(const CompositeVector& u)
{
  u.ScatterMasterToGhosted("cell");
  AMANZI_ASSERT(dkdp_.get() && dkdp_->HasComponent("cell"));
  dkdp_->ScatterMasterToGhosted("cell");

  { // scope for views
    const auto bc_model = bcs_trial_[0]->bc_model();
    const auto bc_value = bcs_trial_[0]->bc_value();

    const auto uc = u.ViewComponent("cell", true);
    const auto cmap_wghost = mesh_->cell_map(true);
    const auto trans_face = transmissibility_->ViewComponent("face", false);

    const auto dKdP_cell = dkdp_->ViewComponent("cell", true);

    const AmanziMesh::Mesh* mesh = mesh_.get();
    DenseMatrix_Vector& A = jac_op_->A; 
    

    Kokkos::parallel_for(
        "PDE_DiffusionFV::AnalyticJacobian_",
        nfaces_owned,
        KOKKOS_LAMBDA(const int& f) {
          AmanziMesh::Entity_ID_View cells;
          mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          int mcells = cells.size();

          auto Aface = A[f];

          double k_rel[2], dkdp[2], pres[2], dist;
          for (int n = 0; n < mcells; n++) {
            pres[n] = uc(cells(n),0);
            dkdp[n] = dKdP_cell(cells(n),0);
          }

          if (mcells == 1) {
            dkdp[1] = 0.;
            //dkdp[1] = dKdP_face.get() ? (*dKdP_face)[0][f] : 0.;
          }

          // find the face direction from cell 0 to cell 1
          AmanziMesh::Entity_ID_View cfaces;

          // Old virtual call... yuck
          // ComputeJacobianLocal_(mcells, f, fdirs[f_index], bc_model[f], bc_value[f],
          //         pres, dkdp, Aface);
          double dKrel_dp[2];
          double dpres;

          if (mcells == 2) {
            dpres = pres[0] - pres[1];
            if (little_k_type_ == OPERATOR_LITTLE_K_UPWIND) {
              double flux0to1;
              flux0to1 = trans_face(f,0) * dpres;
              if (flux0to1  > OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
                dKrel_dp[0] = dkdp[0];
                dKrel_dp[1] = 0.0;
              } else if (flux0to1 < -OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
                dKrel_dp[0] = 0.0;
                dKrel_dp[1] = dkdp[1];
              } else if (fabs(flux0to1) < OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
                dKrel_dp[0] = 0.5 * dkdp[0];
                dKrel_dp[1] = 0.5 * dkdp[1];
              }
            } else if (little_k_type_ == OPERATOR_UPWIND_ARITHMETIC_AVERAGE) {
              dKrel_dp[0] = 0.5 * dkdp[0];
              dKrel_dp[1] = 0.5 * dkdp[1];
            } else {
              assert(false);
            }
            
            Aface(0, 0) = trans_face(f,0) * dpres * dKrel_dp[0];
            Aface(0, 1) = trans_face(f,0) * dpres * dKrel_dp[1];
            Aface(1, 0) = -Aface(0, 0);
            Aface(1, 1) = -Aface(0, 1);

          } else if (mcells == 1) {
            if (bc_model(f) == OPERATOR_BC_DIRICHLET) {
              pres[1] = bc_value(f);
              dpres = pres[0] - pres[1];
              Aface(0, 0) = trans_face(f,0) * dpres * dkdp[0];
            } else {
              Aface(0, 0) = 0.0;
            }
          }
        });
  }
}


#if 0 
/* ******************************************************************
* Computation of a local submatrix of the analytical Jacobian 
* (its nonlinear part) on face f.
****************************************************************** */
void PDE_DiffusionFV::ComputeJacobianLocal_(
    int mcells, int f, int face_dir_0to1, int bc_model_f, double bc_value_f,
    double *pres, double *dkdp_cell, WhetStone::DenseMatrix<>& Jpp)
{
  const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);
}
#endif 

/* ******************************************************************
* Compute transmissibilities on faces 
****************************************************************** */
void PDE_DiffusionFV::ComputeTransmissibility_()
{
  transmissibility_->putScalar(0.);
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

    {
      TensorVector* Kc = K.get();
      Kokkos::parallel_for(
          "PDE_DiffusionFV::ComputeTrans Init Tensor Ceof",
          ncells_owned,
          KOKKOS_LAMBDA(const int& c) {
            Kc->at(c)(0,0) = 1.0;
          });
    }
    K_ = K;
  }

  {
    auto trans_face = transmissibility_->ViewComponent<DefaultDevice>("face", true);
    const Amanzi::AmanziMesh::Mesh* m = mesh_.get();
    const int space_dimension = mesh_->space_dimension(); 

    const TensorVector* K = K_.get();    

    Kokkos::parallel_for(
        "PDE_DiffusionFV::ComputeTransmissibility",
        ncells_owned,
        KOKKOS_LAMBDA(const int c) {

          AmanziMesh::Entity_ID_View faces;
          Kokkos::View<AmanziGeometry::Point*> bisectors;
          m->cell_get_faces_and_bisectors(c, faces, bisectors);

          WhetStone::Tensor<DeviceOnlyMemorySpace> Kc = K->at(c);

          for (int i = 0; i < faces.extent(0); i++) {
            auto f = faces(i);
            const AmanziGeometry::Point& a = bisectors(i);
            const AmanziGeometry::Point& normal = m->face_normal(f);
            const double area = m->face_area(f);

            const double h_tmp = AmanziGeometry::norm(a);
            const double s = area / h_tmp;
            const double perm = ((Kc * a) * normal) * s;
            const double dxn = a * normal;
            Kokkos::atomic_add(&trans_face(f,0), fabs(dxn / perm));
          }
        });
  }
  transmissibility_->GatherGhostedToMaster();
  transmissibility_->reciprocal(*transmissibility_);
  transmissibility_->ScatterMasterToGhosted("face");
  transmissibility_initialized_ = true;
}


}  // namespace Operators
}  // namespace Amanzi

