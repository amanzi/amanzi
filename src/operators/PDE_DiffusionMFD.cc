/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "errors.hh"
#include "LinearOperator.hh"
#include "LinearOperatorFactory.hh"
#include "MatrixFE.hh"
// #include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Diffusion.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "Op.hh"
//#include "Op_Cell_Node.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Face_Cell.hh"
//#include "Op_SurfaceFace_SurfaceCell.hh"

#include "OperatorDefs.hh"
#include "Operator_FaceCell.hh"
//#include "Operator_FaceCellScc.hh"
//#include "Operator_FaceCellSff.hh"
//#include "Operator_Node.hh"
//#include "Operator_ConsistentFace.hh"
#include "UniqueLocalIndex.hh"

#include "PDE_DiffusionMFD.hh"

namespace Amanzi {
namespace Operators {


/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void PDE_DiffusionMFD::SetTensorCoefficient(const Teuchos::RCP<const TensorVector>& K)
{
  std::cout<<"PDE_DiffusionMFD::SetTensorCoefficient"<<std::endl;
  transmissibility_initialized_ = false;
  K_ = K;
  if (local_op_schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
    if (K_ != Teuchos::null && K_.get()) AMANZI_ASSERT(K_->size() == ncells_owned);
    if (!mass_matrices_initialized_) {
      CreateMassMatrices_();
    }
  }
}

/* ******************************************************************
* Initialization of the operator: nonlinear coefficient.
****************************************************************** */
void PDE_DiffusionMFD::SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                            const Teuchos::RCP<const CompositeVector>& dkdp)
{
  std::cout<<"PDE_DiffusionMFD::SetScalarCoefficient"<<std::endl;
  k_ = k;
  dkdp_ = dkdp;

  // compatibility checks
  if (k_ != Teuchos::null) {
    std::cout<<"k_ init"<<std::endl;
    if (little_k_type_ != OPERATOR_LITTLE_K_UPWIND) {
      //AMANZI_ASSERT(k->HasComponent("cell"));
    }

    if (little_k_type_ != OPERATOR_LITTLE_K_STANDARD) {
      AMANZI_ASSERT(k->HasComponent("face"));
    }

    if (little_k_type_ == OPERATOR_LITTLE_K_DIVK_TWIN || 
        little_k_type_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
      AMANZI_ASSERT(k->HasComponent("twin"));
    }

    if (little_k_type_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
      AMANZI_ASSERT(k->HasComponent("grad"));
    }
  }

  // verify that mass matrices were initialized.
  if (!mass_matrices_initialized_) {
    std::cout<<"CreateMassMatrices"<<std::endl;
    CreateMassMatrices_();
  }
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  std::cout<<"PDE_DiffusionMFD::UpdateMatrices"<<std::endl;

  if (k_ != Teuchos::null) k_->ScatterMasterToGhosted();

  if (!exclude_primary_terms_) {
    if (local_op_schema_ & OPERATOR_SCHEMA_DOFS_NODE) {
      assert(false); 
      //UpdateMatricesNodal_();
    } else if ((local_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) &&
               (local_op_schema_ & OPERATOR_SCHEMA_DOFS_FACE)) {
      if (little_k_type_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
        assert(false); 
        //UpdateMatricesMixedWithGrad_(flux);
      } else if (little_k_type_ == OPERATOR_LITTLE_K_NONE) {
        UpdateMatricesMixed_();
      } else {
        assert(false); 
        //UpdateMatricesMixed_little_k_();
      }
    } else if (local_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) {
      assert(false); 
      //UpdateMatricesTPFA_();
    }
  }
}


/* ******************************************************************
* Add stable approximation of Jacobian. It is done typically for 
* the preconditioner.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    double scalar_factor)
{
  assert(false); 
  // add Newton-type corrections
  //if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE) {
  //  if (global_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) {

  //    if (dkdp_ !=  Teuchos::null) dkdp_->ScatterMasterToGhosted();      
  //    AddNewtonCorrectionCell_(flux, u, scalar_factor);
      
  //  } else {
  //    Errors::Message msg("PDE_DiffusionMFD: Newton correction may only be applied to schemas that include CELL dofs.");
  //    Exceptions::amanzi_throw(msg);
  //  }
  //}
}


void PDE_DiffusionMFD::UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& factor)
{
  assert(false); 
  // add Newton-type corrections
  //if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE) {
  //  if (global_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) {

  //    if (dkdp_ !=  Teuchos::null) dkdp_->ScatterMasterToGhosted();      
  //    AddNewtonCorrectionCell_(flux, u, factor);

  //  } else {
  //    Errors::Message msg("PDE_DiffusionMFD: Newton correction may only be applied to schemas that include CELL dofs.");
  //    Exceptions::amanzi_throw(msg);
  //  }
  //}
}  


/* ******************************************************************
* Second-order reconstruction of little k inside mesh cells.
* This member of DIVK-pamily of methods requires to recalcualte all
* mass matrices.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesMixedWithGrad_(
    const Teuchos::Ptr<const CompositeVector>& flux)
{
  AMANZI_ASSERT(!scaled_constraint_);
  assert(false); 
  // preparing little-k data
  //const auto k_cell = k_->ViewComponent("cell");
  //const auto k_face = k_->ViewComponent("face", true);
  //const auto k_grad = k_->ViewComponent("grad");
  //const auto k_twin = k_->ViewComponent("twin", true);
  
  // update matrix blocks
  //int dim = mesh_->space_dimension();
  //WhetStone::MFD3D_Diffusion mfd(mesh_);
  //WhetStone::DenseMatrix<> Wff;

  //AmanziMesh::Entity_ID_View faces, cells;

  //WhetStone::Tensor<> Kc(mesh_->space_dimension(), 1);
  //Kc(0, 0) = 1.0;
  
  //for (int c = 0; c < ncells_owned; c++) {
    // mean value and gradient of nonlinear factor
  //  double kc = k_cell(0,c);
  //  AmanziGeometry::Point kgrad(dim);
  //  for (int i = 0; i < dim; i++) kgrad[i] = k_grad(i,c);
 
    // upwinded values of nonlinear factor
  //  mesh_->cell_get_faces(c, faces);
  //  int nfaces = faces.size();
  //  std::vector<double> kf(nfaces, 1.0); 
  //  if (k_twin.extent(0) != 0) {
  //    for (int n = 0; n < nfaces; n++) kf[n] = k_face(0,faces[n]);
  //  } else {
  //    for (int n = 0; n < nfaces; n++) {
  //      int f = faces[n];
  //      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  //      kf[n] = (c == cells[0]) ? k_face(0,f) : k_twin(0,f);
  //    }
  //  }

  //  if (K_.get()) Kc = K_->at_host(c);
  //  mfd.MassMatrixInverseDivKScaled(c, Kc, kc, kgrad, Wff);

  //  WhetStone::DenseMatrix<> Acell(nfaces + 1, nfaces + 1);

  //  double matsum = 0.0; 
  //  for (int n = 0; n < nfaces; n++) {
  //    double rowsum = 0.0;
  //    for (int m = 0; m < nfaces; m++) {
  //      double tmp = Wff(n, m) * kf[n] * kf[m];
  //      rowsum += tmp;
  //      Acell(n, m) = tmp;
  //    }

  //    Acell(n, nfaces) = -rowsum;
  //    Acell(nfaces, n) = -rowsum;
  //    matsum += rowsum;
  //  }
  //  Acell(nfaces, nfaces) = matsum;
  //  local_op_->matrices[c].assign(Acell);
  //}
}


/* ******************************************************************
* Basic routine for each operator: creation of elemental matrices.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesMixed_()
{
  std::cout<<"PDE_DiffusionMFD::UpdateMatricesMixed_"<<std::endl;
  CSR_Matrix& A = local_op_->A; 

  Kokkos::parallel_for(
    "PDE_DiffusionMFD::UpdateMatricesMixed_",
    ncells_owned, 
    KOKKOS_LAMBDA(const int& c){
      WhetStone::DenseMatrix<DeviceOnlyMemorySpace> Wff(
        Wff_cells_.at(c),Wff_cells_.size(c,0),Wff_cells_.size(c,1));
      int nfaces = Wff.NumRows();

      WhetStone::DenseMatrix<DeviceOnlyMemorySpace> Acell(
        A.at(c),A.size(c,0),A.size(c,1));
      assert(Acell.NumRows() == nfaces && Acell.NumCols() == nfaces); 

      // create stiffness matrix by ellimination of the mass matrix
      double matsum = 0.0;
      for (int n = 0; n < nfaces; n++) {
        double rowsum = 0.0;
        for (int m = 0; m < nfaces; m++) {
          double tmp = Wff(n, m);
          rowsum += tmp;
          Acell(n, m) = tmp;
        }
        Acell(n, nfaces) = -rowsum;
        Acell(nfaces, n) = -rowsum;
        matsum += rowsum;
      }
      Acell(nfaces, nfaces) = matsum;
    }); 
    //local_op_->matrices[c].assign(Acell);
}


/* ******************************************************************
* Basic routine for each operator: creation of elemental matrices.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesMixed_little_k_()
{
  assert(false); 
  //const auto k_cell = k_->ViewComponent("cell");
  //const auto k_face = k_->ViewComponent("face", true);
  //const auto k_twin = k_->ViewComponent("twin", true);

  // update matrix blocks
  //AmanziMesh::Entity_ID_View faces, cells;

  //for (int c = 0; c < ncells_owned; c++) {
  //  mesh_->cell_get_faces(c, faces);
  //  int nfaces = faces.size();

  //  WhetStone::DenseMatrix<> Wff = Wff_cells_[c];
  //  WhetStone::DenseMatrix<> Acell(nfaces + 1, nfaces + 1);

    // Update terms due to nonlinear coefficient
  //  double kc(1.0);
  //  std::vector<double> kf(nfaces, 1.0); 
   
  //  if (k_cell.extent(0) != 0) kc = k_cell(0,c);

    // -- chefs recommendation: SPD discretization with upwind
  //  if (little_k_type_ == OPERATOR_LITTLE_K_DIVK && k_face.extent(0)) {
  //    for (int n = 0; n < nfaces; n++) kf[n] = k_face(0,faces[n]);

    // -- new scheme: SPD discretization with upwind and equal spliting
  //  } else if (little_k_type_ == OPERATOR_LITTLE_K_DIVK_BASE) {
  //    kc = 1.0;
  //    for (int n = 0; n < nfaces; n++) kf[n] = std::sqrt(k_face(0,faces[n]));

    // -- same as above but remains second-order for dicontinuous coefficients
  //  } else if (little_k_type_ == OPERATOR_LITTLE_K_DIVK_TWIN) {
  //    for (int n = 0; n < nfaces; n++) {
  //      int f = faces[n];
  //      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  //      kf[n] = (c == cells[0]) ? k_face(0,f) : k_twin(0,f);
  //    }

    // -- the second most popular choice: classical upwind
  //  } else if (little_k_type_ == OPERATOR_LITTLE_K_UPWIND) {
  //    for (int n = 0; n < nfaces; n++) kf[n] = k_face(0,faces[n]);

  //  } else if (little_k_type_ == OPERATOR_LITTLE_K_STANDARD) {
  //    for (int n = 0; n < nfaces; n++) kf[n] = kc;
  //  }

    // create stiffness matrix by ellimination of the mass matrix
    // -- all methods expect for DIVK-family of methods.
  //  if ((little_k_type_ & OPERATOR_LITTLE_K_DIVK_BASE) == 0) {
      // -- not scaled constraint: kr > 0
  //    if (!scaled_constraint_) {
  //      double matsum = 0.0; 
  //      for (int n = 0; n < nfaces; n++) {
  //        double rowsum = 0.0;
  //        for (int m = 0; m < nfaces; m++) {
  //          double tmp = Wff(n, m) * kf[n];
  //          rowsum += tmp;
  //          Acell(n, m) = tmp;
  //        }

  //        Acell(n, nfaces) = -rowsum;
  //        matsum += rowsum;
  //      }
  //      Acell(nfaces, nfaces) = matsum;

  //      for (int n = 0; n < nfaces; n++) {
  //        double colsum = 0.0;
  //        for (int m = 0; m < nfaces; m++) colsum += Acell(m, n);
  //        Acell(nfaces, n) = -colsum;
  //      }

      // -- scaled constraint: kr >= 0
  //    } else {
  //      double matsum = 0.0;
  //      for (int n = 0; n < nfaces; n++) {
  //        double rowsum = 0.0;
  //        double cur_kf = (kf[n] < scaled_constraint_cutoff_) ? 1.0 : kf[n];
  //        for (int m = 0; m < nfaces; m++) {
  //          double tmp = Wff(n, m) * cur_kf;
  //          rowsum += tmp;
  //          Acell(n, m) = tmp;
  //        }
  //        Acell(n, nfaces) = -rowsum;
  //      }

  //      for (int n = 0; n < nfaces; n++) {
  //        double colsum = 0.0;
  //        for (int m = 0; m < nfaces; m++) colsum += Wff(m, n) * kf[m];
  //        Acell(nfaces, n) = -colsum;
  //        matsum += colsum;
  //      }
  //      Acell(nfaces, nfaces) = matsum;
  //    }
  //  }

    // Amanzi's first upwind: the family of DIVK fmethods
  //  if (little_k_type_ & OPERATOR_LITTLE_K_DIVK_BASE) {
  //    double matsum = 0.0; 
  //    for (int n = 0; n < nfaces; n++) {
  //      double rowsum = 0.0;
  //      for (int m = 0; m < nfaces; m++) {
  //        double tmp = Wff(n, m) * kf[n] * kf[m] / kc;
  //        rowsum += tmp;
  //        Acell(n, m) = tmp;
  //      }
          
  //      Acell(n, nfaces) = -rowsum;
  //      Acell(nfaces, n) = -rowsum;
  //      matsum += rowsum;
  //    }

  //    Acell(nfaces, nfaces) = matsum;
  //  }
    
  //  local_op_->matrices[c].assign(Acell);
  //}
}


/* ******************************************************************
* Calculate elemental stiffness matrices: nodal DOFs.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesNodal_()
{
  AMANZI_ASSERT(!scaled_constraint_);
  assert(false); 
  // update matrix blocks
  //WhetStone::MFD3D_Diffusion mfd(mesh_);
  //mfd.ModifyStabilityScalingFactor(factor_);

  //AmanziMesh::Entity_ID_List nodes;

  //nfailed_primary_ = 0;

  //WhetStone::Tensor<> K(2, 1);
  //K(0, 0) = 1.0;
  
  //for (int c = 0; c < ncells_owned; c++) {
  //  if (K_.get()) K = K_->at_host(c);

  //  mesh_->cell_get_nodes(c, nodes);
  //  int nnodes = nodes.size();

  //  WhetStone::DenseMatrix<> Acell(nnodes, nnodes);

  //  int method = mfd_primary_;
  //  int ok = WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  //  if (method == WhetStone::DIFFUSION_OPTIMIZED_FOR_MONOTONICITY) {
  //    ok = mfd.StiffnessMatrixMMatrix(c, K, Acell);
  //    method = mfd_secondary_;
  //  } else {
  //    ok = mfd.StiffnessMatrix(c, K, Acell);
  //    method = mfd_secondary_;
  //  }

  //  if (ok != WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) {
  //    nfailed_primary_++;
  //    ok = mfd.StiffnessMatrix(c, K, Acell);
  //  }

  //  if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
  //    Errors::Message msg("Stiffness_MFD: unexpected failure of LAPACK in WhetStone.");
  //    Exceptions::amanzi_throw(msg);
  //  }

  //  local_op_->matrices[c].assign(Acell);
  //} 
}

/* ******************************************************************
* Calculate and assemble fluxes using the TPFA scheme.
* This routine does not use little k.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesTPFA_()
{
  assert(false); 
  // populate transmissibilities
  //WhetStone::MFD3D_Diffusion mfd(mesh_);
  //WhetStone::DenseMatrix<> Mff; 

  // \TODO Allocate memory
  //CSR_Matrix Mffv;

  //CSR_Matrix& A = local_op_->A; 
  
  //CompositeVectorSpace cv_space;
  //cv_space.SetMesh(mesh_);
  //cv_space.SetGhosted(true);
  //cv_space.SetComponent("face", AmanziMesh::FACE, 1);

  //Teuchos::RCP<CompositeVector> T = cv_space.Create(); //Teuchos::RCP<CompositeVector>(new CompositeVector(cv_space, true));
  //const auto Ttmp = T->ViewComponent("face", true);

  //WhetStone::Tensor<> Kc(mesh_->space_dimension(), 1);
  //Kc(0, 0) = 1.0;
  //const Amanzi::AmanziMesh::Mesh* m = mesh_.get(); 

  //AmanziMesh::Entity_ID_View cells, faces;
  //Ttmp.PutScalar(0.0);

  //Kokkos::parallel_for(
  //  "PDE_DiffusionMFD::UpdateMatricesTPFA_", 
  //  ncells_owned, 
  //  KOKKOS_LAMBDA(const int& c){
  //    AmanziMesh::Entity_ID_View faces;
  //    mesh_->cell_get_faces(c, faces);
  //    int nfaces = faces.size();

  //    WhetStone::Tensor<DeviceOnlyMemorySpace> Kc = K_->at(c); 
  //    WhetStone::DenseMatrix<DeviceOnlyMemorySpace> Mff(Mffv.at(c),Mffv.size(c,0),Mffv.size(c,1)); 

      //mfd.MassMatrixInverseTPFA(c, Kc, Mff);
   
  //    for (int n = 0; n < nfaces; n++) {
  //      int f = faces[n];
  //      Ttmp(0,f) += 1.0 / Mff(n, n);
  //    }
  //  });

  //Kokkos::parallel_for(
  //  "PDE_DiffusionMFD::UpdateMatricesTPFA_",
  //  nfaces_owned,
  //  KOKKOS_LAMBDA(const int& f){
  //    AmanziMesh::Entity_ID_View cells;
  //    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  //    int ncells = cells.size();
  //    WhetStone::DenseMatrix<DeviceOnlyMemorySpace> Aface(A.at(f),A.size(f,0),A.size(f,1)); 

      //if (Ttmp(0,f) == 0.0) {
      //  Aface = 0.0;
      //  local_op_->matrices[f].assign(Aface);
      //  continue;  // We skip zero transmissibilities
      //}

  //    if (ncells == 2) {
  //      double coef = 1.0 / Ttmp(0,f);
  //      Aface(0, 0) =  coef;
  //      Aface(1, 1) =  coef;
  //      Aface(0, 1) = -coef;
  //      Aface(1, 0) = -coef;
  //    } else {
  //      double coef = 1.0 / Ttmp(0,f);
  //      Aface(0, 0) = coef;
  //    }

  //  });  
}


/* ******************************************************************
* Apply boundary conditions to the local matrices. We always zero-out
* matrix rows for essential test BCs. As to trial BCs, there are
* options: eliminate them or not. Finally we may add the essntial BC
* the the system of equations as the trivial equations.
*
* NOTE 1. Nodal scheme handles only the case trialBC = testBC.
* NOTE 2. Jacobian term handles only trial BCs.
****************************************************************** */
void PDE_DiffusionMFD::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  std::cout<<"PDE_DiffusionMFD::ApplyBCs"<<std::endl;

  if (!exclude_primary_terms_) {
    if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL
                           | OPERATOR_SCHEMA_DOFS_FACE
                           | OPERATOR_SCHEMA_DOFS_CELL)) {
      AMANZI_ASSERT(bcs_trial_.size() == 1);
      AMANZI_ASSERT(bcs_test_.size() == 1);
      ApplyBCs_Mixed_(bcs_trial_[0].ptr(), bcs_test_[0].ptr(), primary, eliminate, essential_eqn);
    
    } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_FACE
                                  | OPERATOR_SCHEMA_DOFS_CELL)) {
      assert(false); 
      //AMANZI_ASSERT(bcs_trial_.size() == 1);
      //AMANZI_ASSERT(bcs_test_.size() == 1);
      //ApplyBCs_Cell_(bcs_trial_[0].ptr(), bcs_test_[0].ptr(), primary, eliminate, essential_eqn);
    
    } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL
                                  | OPERATOR_SCHEMA_DOFS_NODE)) {
      assert(false); 
      //Teuchos::Ptr<const BCs> bc_f, bc_n;
      //for (const auto& bc : bcs_trial_) {
      //  if (bc->kind() == AmanziMesh::FACE) {
      //    bc_f = bc.ptr();
      //  } else if (bc->kind() == AmanziMesh::NODE) {
      //    bc_n = bc.ptr();
      //  }
      //}
      //ApplyBCs_Nodal_(bc_f.ptr(), bc_n.ptr(), primary, eliminate, essential_eqn);
    }
  }

  if (jac_op_ != Teuchos::null) {
    assert(false); 
    //const auto bc_model = bcs_trial_[0]->bc_model();
    //AMANZI_ASSERT(bc_model.size() == nfaces_wghost);
    //for (int f = 0; f != nfaces_owned; ++f) {
    //  WhetStone::DenseMatrix<> Aface = jac_op_->matrices[f];
    //  if (bc_model[f] == OPERATOR_BC_NEUMANN ||
    //      bc_model[f] == OPERATOR_BC_TOTAL_FLUX) {
    //    jac_op_->matrices_shadow[f].assign(Aface);
    //    Aface *= 0.0;
    //  }
    //}
  }
}


/* ******************************************************************
* Apply BCs on face values.
****************************************************************** */
void PDE_DiffusionMFD::ApplyBCs_Mixed_(
    const Teuchos::Ptr<const BCs>& bc_trial,
    const Teuchos::Ptr<const BCs>& bc_test,
    bool primary, bool eliminate, bool essential_eqn)
{
  std::cout<<"PDE_DiffusionMFD::ApplyBCs_Mixed_"<<std::endl;

  // apply diffusion type BCs to FACE-CELL system
  const Amanzi::AmanziMesh::Mesh* m = mesh_.get();
  CSR_Matrix& A = local_op_->A; 

  const auto bc_model_trial = bc_trial->bc_model();
  const auto bc_model_test = bc_test->bc_model();

  const auto bc_value = bc_trial->bc_value();
  const auto bc_mixed = bc_trial->bc_mixed();

  AMANZI_ASSERT(bc_model_trial.size() == nfaces_wghost);
  AMANZI_ASSERT(bc_value.size() == nfaces_wghost);

  //global_op_->rhs()->PutScalarGhosted(0.0);
  const auto rhs_face = global_op_->rhs()->ViewComponent("face", true);
  const auto rhs_cell = global_op_->rhs()->ViewComponent("cell");

  //const auto k_cell = k_->ViewComponent("cell");
  const auto k_face = k_->ViewComponent("face", true);
  std::cout<<"k_face: "<<k_face.extent(0)<<";"<<k_face.extent(1)<<std::endl;

  Kokkos::parallel_for(
    "PDE_DiffusionMFD::ApplyBCs_Mixed_",
    ncells_owned,
    KOKKOS_LAMBDA(const int& c){
      AmanziMesh::Entity_ID_View faces;

      m->cell_get_faces(c, faces);
      int nfaces = faces.size();
      
      // Update terms due to nonlinear coefficient
      double kc(1.0);
      //std::vector<double> kf(nfaces, 1.0);
      double kfa(1.0); 

      if (scaled_constraint_) {
        assert(false); 
        // un-rolling little-k data
        //if (k_cell.extent(0) != 0) kc = k_cell(0,c);
        
        if (little_k_type_ == OPERATOR_LITTLE_K_UPWIND) {
          assert(false); 
          //#if 0 
          //for (int n = 0; n < nfaces; n++) kf[n] = k_face(0,faces[n]);
          //#endif 
          
        } else if (little_k_type_ == OPERATOR_LITTLE_K_STANDARD) {
          assert(false); 
          //#if 0 
          //for (int n = 0; n < nfaces; n++) kf[n] = kc;
          //#endif 
        }
      }
      
      bool flag(true);
      WhetStone::DenseMatrix<DeviceOnlyMemorySpace> Acell(A.at(c),A.size(c,0),A.size(c,1));
      // essential conditions for test functions
      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        if (bc_model_test[f] == OPERATOR_BC_DIRICHLET) { 
          //if (flag) {  // make a copy of elemental matrix
          //  local_op_->matrices_shadow[c].assign(Acell);
          //  flag = false;
          //}
          for (int m = 0; m < nfaces + 1; m++) Acell(n, m) = 0.0; 
        }
      }

      // conditions for trial functions
      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        double value = bc_value[f];

        if (bc_model_trial[f] == OPERATOR_BC_DIRICHLET) { 
          // make a copy of elemental matrix for post-processing
          //if (flag) {
          //  local_op_->matrices_shadow[c].assign(Acell);
          //  flag = false;
          //}

          if (eliminate) { 
            for (int m = 0; m < nfaces; m++) {
              //rhs_face(0,faces[m]) -= Acell(m, n) * value;
              Kokkos::atomic_add(&rhs_face(0,faces[m]),-Acell(m, n) * value); 
              Acell(m, n) = 0.0;
            }

            rhs_cell(0,c) -= Acell(nfaces, n) * value;
            Acell(nfaces, n) = 0.0;
          }

          if (essential_eqn) {
            rhs_face(0,f) = value;
            Acell(n, n) = 1.0;
          } 
        } else if (bc_model_trial[f] == OPERATOR_BC_NEUMANN && primary) {
          assert(false); 
          //#if 0 
          //if (scaled_constraint_) {
          //  if (std::abs(kf[n]) < scaled_constraint_fuzzy_) {
          //    AMANZI_ASSERT(value == 0.0);
          //    rhs_face(0,f) = 0.0;
          //  } else if (kf[n] < scaled_constraint_cutoff_) {
          //    rhs_face(0,f) -= value * mesh_->face_area(f) / kf[n];
          //  } else {
          //    rhs_face(0,f) -= value * mesh_->face_area(f);
          //  }
          //} else {
          //  rhs_face(0,f) -= value * mesh_->face_area(f);
          //}
          //#endif 
        } else if (bc_model_trial[f] == OPERATOR_BC_TOTAL_FLUX && primary) {
          assert(false); 
          //#if 0 
          //if (scaled_constraint_ && kf[n] < scaled_constraint_cutoff_) {
          //  AMANZI_ASSERT(false);
          //} else {
          //  rhs_face(0,f) -= value * mesh_->face_area(f);
          //}
          //#endif 
        } else if (bc_model_trial[f] == OPERATOR_BC_MIXED && primary) {
          //if (flag) {  // make a copy of elemental matrix
          //  local_op_->matrices_shadow[c].assign(Acell);
          //  flag = false;
          //}
          double area = m->face_area(f);
          if (scaled_constraint_) {
            if (std::abs(kfa) < scaled_constraint_fuzzy_) {
              assert((value == 0.0) && (bc_mixed[f] == 0.0));
              //rhs_face(0,f) = 0.0;
              //Kokkos::atomic_store(&rhs_face(0,f),0.0); // NOT FOUND?
              Kokkos::atomic_exchange(&rhs_face(0,f),0.0); 
              
            } else if (kfa < scaled_constraint_cutoff_) {
              //rhs_face(0,f) -= value * area / kfa;
              Kokkos::atomic_add(&rhs_face(0,f),- value * area / kfa); 
              Acell(n, n) += bc_mixed[f] * area / kfa;
            } else {
              //rhs_face(0,f) -= value * area;
              Kokkos::atomic_add(&rhs_face(0,f),- value * area); 
              Acell(n, n) += bc_mixed[f] * area;
            }
          } else {
            //rhs_face(0,f) -= value * area;
            Kokkos::atomic_add(&rhs_face(0,f),- value * area); 
            Acell(n, n) += bc_mixed[f] * area;
          }
        }
      }
    }); 
  //global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* Apply BCs on cell operators
****************************************************************** */
void PDE_DiffusionMFD::ApplyBCs_Cell_(
   const Teuchos::Ptr<const BCs>& bc_trial,
   const Teuchos::Ptr<const BCs>& bc_test,
   bool primary, bool eliminate, bool essential_eqn)
{
  assert(false); 
  // apply diffusion type BCs to CELL system
  //AmanziMesh::Entity_ID_View cells;

  //const auto bc_model = bc_trial->bc_model();
  //const auto bc_value = bc_trial->bc_value();
  //const auto bc_mixed = bc_trial->bc_mixed();

  //AMANZI_ASSERT(bc_model.size() == nfaces_wghost);
  //AMANZI_ASSERT(bc_value.size() == nfaces_wghost);

  //const auto rhs_cell = global_op_->rhs()->ViewComponent("cell");
    
  //for (int f = 0; f != nfaces_owned; ++f) {
  //  WhetStone::DenseMatrix<> Aface = local_op_->matrices[f];
      
  //  if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
  //    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  //    rhs_cell(0,cells[0]) += bc_value[f] * Aface(0, 0);
  //  }
    // Neumann condition contributes to the RHS
  //  else if (bc_model[f] == OPERATOR_BC_NEUMANN && primary) {
  //    local_op_->matrices_shadow[f].assign(Aface);
      
  //    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  //    rhs_cell(0,cells[0]) -= bc_value[f] * mesh_->face_area(f);
  //    Aface *= 0.0;
  //  }
    // solve system of two equations in three unknowns
  //  else if (bc_model[f] == OPERATOR_BC_MIXED && primary) {
  //    local_op_->matrices_shadow[f].assign(Aface);
      
  //    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  //    double area = mesh_->face_area(f);
  //    double factor = area / (1.0 + bc_mixed[f] * area / Aface(0, 0));
  //    rhs_cell(0,cells[0]) -= bc_value[f] * factor;
  //    Aface(0, 0) = bc_mixed[f] * factor;
  //  }
  //}
}


/* ******************************************************************
* Apply BCs on nodal operators
****************************************************************** */
void PDE_DiffusionMFD::ApplyBCs_Nodal_(
    const Teuchos::Ptr<const BCs>& bc_f,
    const Teuchos::Ptr<const BCs>& bc_v,
    bool primary, bool eliminate, bool essential_eqn)
{
  assert(false); 
  //AmanziMesh::Entity_ID_View faces;
  //AmanziMesh::Entity_ID_List nodes, cells; 

  //global_op_->rhs()->PutScalarGhosted(0.0);
  //const auto rhs_node = global_op_->rhs()->ViewComponent("node", true);

  //int nn(0), nm(0);
  //for (int c = 0; c != ncells_owned; ++c) {
  //  bool flag(true);
  //  WhetStone::DenseMatrix<> Acell = local_op_->matrices[c];

    // process boundary integrals
  //  if (bc_f != Teuchos::null) {
  //    const auto bc_model = bc_f->bc_model();
  //    const auto bc_value = bc_f->bc_value();
  //    const auto bc_mixed = bc_f->bc_mixed();

  //    mesh_->cell_get_faces(c, faces);
  //    int nfaces = faces.size();

  //    for (int n = 0; n != nfaces; ++n) {
  //      int f = faces[n];

  //      if (bc_model[f] == OPERATOR_BC_NEUMANN && primary) {
  //        nn++;
  //        double value = bc_value[f];
  //        double area = mesh_->face_area(f);

  //        mesh_->face_get_nodes(f, nodes);
  //        int nnodes = nodes.size();

  //        for (int m = 0; m < nnodes; m++) {
  //          int v = nodes[m];
  //          if (bc_v->bc_model()[v] != OPERATOR_BC_DIRICHLET)
  //            rhs_node(0,v) -= value * area / nnodes;
  //        }
  //      } else if (bc_model[f] == OPERATOR_BC_MIXED && primary) {
  //        nm++;
  //        if (flag) {  // make a copy of cell-based matrix
  //          local_op_->matrices_shadow[c].assign(Acell);
  //          flag = false;
  //        }
  //        double value = bc_value[f];
  //        double area = mesh_->face_area(f);

  //        mesh_->face_get_nodes(f, nodes);
  //        int nnodes = nodes.size();

  //        for (int m = 0; m < nnodes; m++) {
  //          int v = nodes[m];
  //          if (bc_v->bc_model()[v] != OPERATOR_BC_DIRICHLET)
  //            rhs_node(0,v) -= value * area / nnodes;
  //          Acell(n, n) += bc_mixed[f] * area / nnodes;
  //        }
  //      }
  //    }
  //  } 

  //  if (bc_v != Teuchos::null) {
  //    const auto bc_model = bc_v->bc_model();
  //    const auto bc_value = bc_v->bc_value();

  //    mesh_->cell_get_nodes(c, nodes);
  //    int nnodes = nodes.size();

      // essential conditions for test functions
  //    for (int n = 0; n != nnodes; ++n) {
  //      int v = nodes[n];
  //      if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
  //        if (flag) {  // make a copy of elemental matrix
  //          local_op_->matrices_shadow[c].assign(Acell);
  //          flag = false;
  //        }
  //        for (int m = 0; m < nnodes; m++) Acell(n, m) = 0.0;
  //      }
  //    }

  //    for (int n = 0; n != nnodes; ++n) {
  //      int v = nodes[n];
  //      double value = bc_value[v];

  //      if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
  //        if (flag) {  // make a copy of cell-based matrix
  //          local_op_->matrices_shadow[c].assign(Acell);
  //          flag = false;
  //        }
     
  //        if (eliminate) {
  //          for (int m = 0; m < nnodes; m++) {
  //            rhs_node(0,nodes[m]) -= Acell(m, n) * value;
  //            Acell(m, n) = 0.0;
  //          }
  //        }

          // We take into account multiple contributions to matrix diagonal
          // by dividing by the number of cells attached to a vertex.
  //        if (essential_eqn) {
  //          mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, cells);
  //          if (v < nnodes_owned) rhs_node(0,v) = value;
  //          Acell(n, n) = 1.0 / cells.size();
  //        }
  //      }
  //    }
  //  }
  //} 
  //global_op_->rhs()->GatherGhostedToMaster("node", Add);
}


/* ******************************************************************
* Modify operator by adding upwind approximation of Newton corection.
* A special care should be taken later to deal with the case 
* where kf < 0, i.e. for energy when: div (qh) = div (h k grad p), 
* where h is enthalpy and can be negative. I think that the current
* treatment is inadequate.
****************************************************************** */
void PDE_DiffusionMFD::AddNewtonCorrectionCell_(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    double scalar_factor)
{
  assert(false); 
  // hack: ignore correction if no flux provided.
  //if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  //if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

  // only works on upwinded methods
  //if (little_k_type_ == OPERATOR_UPWIND_NONE) return;

  //const auto kf = k_->ViewComponent("face");
  //const auto dkdp_f = dkdp_->ViewComponent("face");
  //const auto flux_f = flux->ViewComponent("face");

  // populate the local matrices
  //double v, vmod;
  //AmanziMesh::Entity_ID_View cells;
  //for (int f = 0; f < nfaces_owned; f++) {
  //  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  //  int ncells = cells.size();
  //  WhetStone::DenseMatrix<> Aface(ncells, ncells);
  //  //Aface.PutScalar(0.0);

    // We use the upwind discretization of the generalized flux.
  //  v = std::abs(kf(0,f)) > 0.0 ? flux_f(0,f) * dkdp_f(0,f) / kf(0,f) : 0.0;
  //  vmod = std::abs(v);

    // prototype for future limiters (external or internal ?)
  //  vmod *= scalar_factor;

    // define the upwind cell, index i in this case
  //  int i, dir, c1;
  //  c1 = cells[0];
  //  const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c1, &dir);
  //  i = (v * dir >= 0.0) ? 0 : 1;

  //  if (ncells == 2) {
  //    Aface(i, i) = vmod;
  //    Aface(1 - i, i) = -vmod;
  //  } else if (i == 0) {
  //    Aface(0, 0) = vmod;
  //  }

  //  jac_op_->matrices[f].assign(Aface);
  //}
}


void PDE_DiffusionMFD::AddNewtonCorrectionCell_(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& factor)
{
  assert(false); 
  // hack: ignore correction if no flux provided.
  //if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  //if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

  // only works on upwinded methods
  //if (little_k_type_ == OPERATOR_UPWIND_NONE) return;

  //const auto kf = k_->ViewComponent("face");
  //const auto dkdp_f = dkdp_->ViewComponent("face");
  //const auto flux_f = flux->ViewComponent("face");
  //const auto factor_cell = factor->ViewComponent("cell");

  // populate the local matrices
  //double v, vmod;
  //AmanziMesh::Entity_ID_View cells;
  //for (int f = 0; f < nfaces_owned; f++) {
  //  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  //  int ncells = cells.size();
  //  WhetStone::DenseMatrix<> Aface(ncells, ncells);
  //  //Aface.PutScalar(0.0);

    // We use the upwind discretization of the generalized flux.
  //  v = std::abs(kf(0,f)) > 0.0 ? flux_f(0,f) * dkdp_f(0,f) / kf(0,f) : 0.0;
  //  vmod = std::abs(v);

  //  double scalar_factor=0.;
  //  for (int j=0; j<ncells; j++) scalar_factor += factor_cell(0,cells[j]);
  //  scalar_factor *= 1./ncells;

    // prototype for future limiters (external or internal ?)
  //  vmod *= scalar_factor;

    // define the upwind cell, index i in this case
  //  int i, dir, c1;
  //  c1 = cells[0];
  //  const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c1, &dir);
  //  i = (v * dir >= 0.0) ? 0 : 1;

  //  if (ncells == 2) {
  //    Aface(i, i) = vmod;
  //    Aface(1 - i, i) = -vmod;
  //  } else if (i == 0) {
  //    Aface(0, 0) = vmod;
  //  }

  //  jac_op_->matrices[f].assign(Aface);
  //}
}  


  
/* ******************************************************************
* Given pressures, reduce the problem to Lagrange multipliers.
****************************************************************** */
void PDE_DiffusionMFD::ModifyMatrices(const CompositeVector& u)
{
  assert(false); 
  //if (local_op_schema_ != (OPERATOR_SCHEMA_BASE_CELL |
  //                         OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
  //  std::cout << "Schema " << global_op_schema_ << " is not supported" << std::endl;
  //  AMANZI_ASSERT(0);
  //}

  // populate the matrix
  //AmanziMesh::Entity_ID_View faces;
  //const auto u_c = u.ViewComponent("cell");

  //global_op_->rhs()->PutScalarGhosted(0.0);

  //const auto rhs_f = global_op_->rhs()->ViewComponent("face", true);
  //for (int c = 0; c != ncells_owned; ++c) {
  //  mesh_->cell_get_faces(c, faces);
  //  int nfaces = faces.size();

  //  WhetStone::DenseMatrix<> Acell = local_op_->matrices[c];

  //  for (int n = 0; n < nfaces; n++) {
  //    int f = faces[n];
  //    rhs_f(0,f) -= Acell(n, nfaces) * u_c(0,c);
  //    Acell(n, nfaces) = 0.0;
  //    Acell(nfaces, n) = 0.0;
  //  }
  //} 
  // Assemble all right-hand sides
  //global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* WARNING: Since diffusive flux may be discontinuous (e.g. for
* Richards equation), we derive it in exactly the same manner as 
* in gravity routines.
* **************************************************************** */
void PDE_DiffusionMFD::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                  const Teuchos::Ptr<CompositeVector>& flux)
{
  std::cout<<"PDS_DiffusionMFD::UpdateFlux"<<std::endl;
  // Initialize intensity in ghost faces.
  //flux->PutScalar(0.0);
  //u->ScatterMasterToGhosted("face");

  //if (k_ != Teuchos::null) {
  //  if (k_->HasComponent("face")) k_->ScatterMasterToGhosted("face");
  //}

  const auto u_cell = u->ViewComponent("cell");
  const auto u_face = u->ViewComponent("face", true);
  const auto flux_data = flux->ViewComponent("face", true);
  const Amanzi::AmanziMesh::Mesh* m = mesh_.get();


  Kokkos::View<int*,DeviceOnlyMemorySpace> hits("hits",nfaces_wghost); 
  //std::vector<int> hits(nfaces_wghost, 0);

  auto local_A = local_op_->A; 
  auto local_Av = local_op_->Av; 
  auto local_v = local_op_->v; 

  Kokkos::parallel_for(
    "PDE_DiffusionMFD::UpdateFlux",
    ncells_owned, 
    KOKKOS_LAMBDA(const int& c){
      AmanziMesh::Entity_ID_View faces;
      Kokkos::View<int*> dirs;
      m->cell_get_faces_and_dirs(c, faces, dirs);
      int nfaces = faces.size();

      WhetStone::DenseVector<Amanzi::DeviceOnlyMemorySpace> lv = getFromCSR<WhetStone::DenseVector>(local_v,c);
      WhetStone::DenseVector<Amanzi::DeviceOnlyMemorySpace> lAv = getFromCSR<WhetStone::DenseVector>(local_Av,c);
      WhetStone::DenseMatrix<Amanzi::DeviceOnlyMemorySpace> lA = getFromCSR<WhetStone::DenseMatrix>(local_A,c);

      //WhetStone::DenseVector<> v(nfaces + 1), av(nfaces + 1);

      for (int n = 0; n < nfaces; n++) {
        lv(n) = u_face(0,faces[n]);
      }
      lv(nfaces) = u_cell(0,c);

      lA.Multiply(lv, lAv, false);

      //if (local_op_->matrices_shadow[c].NumRows() == 0) { 
      //  local_op_->matrices[c].Multiply(v, av, false);
      //} else {
      //  local_op_->matrices_shadow[c].Multiply(v, av, false);
      //}

      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        if (f < nfaces_owned) {
          flux_data(0,f) -= lAv(n) * dirs[n];
          Kokkos::atomic_add(&hits(f), 1);
        }
      }
  });

  Kokkos::parallel_for(
    "PDE_DiffusionMFD::UpdateFlux",
    nfaces_owned, 
    KOKKOS_LAMBDA(const int& f){
      flux_data(0,f) /= hits[f];
    });
    //for (int f = 0; f != nfaces_owned; ++f) {
    //  flux_data(0,f) /= hits[f];
    //}

  //for (int c = 0; c < ncells_owned; c++) {
    //mesh_->cell_get_faces_and_dirs(c, faces, dirs);
    //nt nfaces = faces.size();

    //WhetStone::DenseVector<> v(nfaces + 1), av(nfaces + 1);
    //for (int n = 0; n < nfaces; n++) {
    //  v(n) = u_face(0,faces[n]);
    //}
    //v(nfaces) = u_cell(0,c);

    //if (local_op_->matrices_shadow[c].NumRows() == 0) { 
    //  local_op_->matrices[c].Multiply(v, av, false);
    //} else {
    //  local_op_->matrices_shadow[c].Multiply(v, av, false);
    //}

    //for (int n = 0; n < nfaces; n++) {
    //  int f = faces[n];
    //  if (f < nfaces_owned) {
    //    flux_data(0,f) -= av(n) * dirs[n];
    //    hits[f]++;
    //  }
    //}
  //}

  //for (int f = 0; f != nfaces_owned; ++f) {
  //  flux_data(0,f) /= hits[f];
  //}
}


/* ******************************************************************
* Calculates one-sided (cell-based) fluxes that satisfy proper 
* continuity conditions. This differs from other subroutines
* calculing fluxes due to presence of multiple normals on some
* non-manifold faces.
* **************************************************************** */
void PDE_DiffusionMFD::UpdateFluxNonManifold(
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<CompositeVector>& flux)
{
  assert(false); 
  // Initialize intensity in ghost faces.
  //u->ScatterMasterToGhosted("face");

  //const auto u_cell = u->ViewComponent("cell");
  //const auto u_face = u->ViewComponent("face", true);
  //const auto flux_data = flux->ViewComponent("face", true);

  //flux_data.PutScalar(0.0);

  //int ndofs_owned = flux->ViewComponent("face").extent(0);
  //int ndofs_wghost = flux_data.extent(0);

  //AmanziMesh::Entity_ID_View faces;
  //Kokkos::View<int*> dirs;
  //const auto fmap = flux->getMap()->getMap("face", true);

  //for (int c = 0; c < ncells_owned; c++) {
  //  mesh_->cell_get_faces_and_dirs(c, faces, dirs);
  //  int nfaces = faces.size();

  //  WhetStone::DenseVector<> v(nfaces + 1), av(nfaces + 1);
  //  for (int n = 0; n < nfaces; n++) {
  //    v(n) = u_face(0,faces[n]);
  //  }
  //  v(nfaces) = u_cell(0,c);

  //  if (local_op_->matrices_shadow[c].NumRows() == 0) { 
  //    local_op_->matrices[c].Multiply(v, av, false);
  //  } else {
  //    local_op_->matrices_shadow[c].Multiply(v, av, false);
  //  }

  //  for (int n = 0; n < nfaces; n++) {
  //    int f = faces[n];
  //    int g = fmap.FirstPointInElement(f);

  //    int ndofs = fmap.ElementSize(f);
  //    if (ndofs > 1) g += Operators::UniqueIndexFaceToCells(*mesh_, f, c);

  //    flux_data(0,g) -= av(n) * dirs[n];
  //  }
  //}
  //flux->GatherGhostedToMaster(Add);
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.
****************************************************************** */
void PDE_DiffusionMFD::CreateMassMatrices_()
{
  std::cout<<"PDE_DiffusionMFD::CreateMassMatrices"<<std::endl;
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  WhetStone::DenseMatrix<> Wff;
  bool surface_mesh = (mesh_->manifold_dimension() != mesh_->space_dimension());

  mfd.ModifyStabilityScalingFactor(factor_);

  std::vector<WhetStone::DenseMatrix<>> tmp_Wff_cells; 

  tmp_Wff_cells.resize(ncells_owned);

  WhetStone::Tensor<> Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;

  for (int c = 0; c < ncells_owned; c++) {
    int ok;
    if (K_.get()) Kc = K_->at_host(c);

    // For problems with degenerate coefficients we should skip WhetStone.
    if (Kc.Trace() == 0.0) {
      AMANZI_ASSERT(0);
      // int nfaces = mesh_->cell_get_num_faces(c);
      // Wff.Reshape(nfaces, nfaces);
      // Wff.PutScalar(0.0);
      // ok = WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK;
    } else if (surface_mesh) {
      AMANZI_ASSERT(0);
      // ok = mfd.MassMatrixInverseSurface(c, Kc, Wff);
    } else {
      int method = mfd_primary_;
      ok = WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED;

      // try primary and then secondary discretization methods.
      if (method == WhetStone::DIFFUSION_HEXAHEDRA_MONOTONE) {
        AMANZI_ASSERT(0);
        // ok = mfd.MassMatrixInverseMMatrixHex(c, Kc, Wff);
        // method = mfd_secondary_;
      } else if (method == WhetStone::DIFFUSION_OPTIMIZED_FOR_MONOTONICITY) {
        assert(false); 
        //ok = mfd.MassMatrixInverseMMatrix(c, Kc, Wff);
        //method = mfd_secondary_;
      }

      if (ok != WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) {
        if (method == WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY) {
          assert(false); 
          //ok = mfd.MassMatrixInverseOptimized(c, Kc, Wff);
        } else if (method == WhetStone::DIFFUSION_TPFA) {
          ok = mfd.MassMatrixInverseTPFA(c, Kc, Wff);
        } else if (method == WhetStone::DIFFUSION_SUPPORT_OPERATOR) {
          AMANZI_ASSERT(0);
          // ok = mfd.MassMatrixInverseSO(c, Kc, Wff);
        } else if (method == WhetStone::DIFFUSION_POLYHEDRA_SCALED) {
          if (K_symmetric_) {
            ok = mfd.MassMatrixInverse(c, Kc, Wff);
          } else {
            AMANZI_ASSERT(0);
            // ok = mfd.MassMatrixInverseNonSymmetric(c, Kc, Wff);
          }
        }
      }
    }

    tmp_Wff_cells[c].reshape(Wff.NumRows(),Wff.NumCols()); 
    tmp_Wff_cells[c].assign(Wff);

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("PDE_DiffusionMFD: unexpected failure in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
  }
  // Copy data to GPU 
  Wff_cells_ = CSR_Matrix(ncells_owned);
  for(int c = 0 ; c < ncells_owned; ++c){
    //int n = tmp_Wff_cells[c].NumRows();
    //int m = tmp_Wff_Cells[c].NumCols(); 
    int loc[2] = {tmp_Wff_cells[c].NumRows(),tmp_Wff_cells[c].NumCols()};
    Wff_cells_.set_shape_host(c,loc); 
  }
  Wff_cells_.prefix_sum();  
  // Copy data 
  int k = 0 ; 
  for(int c = 0 ; c < ncells_owned; ++c){
    for(int i = 0 ; i < tmp_Wff_cells[c].Values().extent(0) ; ++i)
      Wff_cells_.set_entries_host(k++,tmp_Wff_cells[c].Values()[i]);
  }
  Wff_cells_.update_entries_device();

  mass_matrices_initialized_ = true;
}


/* ******************************************************************
* Scale elemental inverse mass matrices. Use case is saturated flow.
****************************************************************** */
void PDE_DiffusionMFD::ScaleMassMatrices(double s)
{
  std::cout<<"PDE_DiffusionMFD::ScaleMassMatrices"<<std::endl;
  assert(Wff_cells_.size() == ncells_owned); 

  Kokkos::parallel_for(
    "PDE_DiffusionMFD::ScaleMassMatrices",
    ncells_owned,
    KOKKOS_LAMBDA(const int& c){
      WhetStone::DenseMatrix<DeviceOnlyMemorySpace> lm(
        Wff_cells_.at(c),Wff_cells_.size(c,0),Wff_cells_.size(c,1)); 
      lm *= s; 
    });
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void PDE_DiffusionMFD::ParsePList_()
{
  std::cout<<"PDE_DiffusionMFD::ParsePList_"<<std::endl;

  // Determine discretization
  std::string primary = plist_.get<std::string>("discretization primary");
  std::string secondary = plist_.get<std::string>("discretization secondary", primary);
  K_symmetric_ = (plist_.get<std::string>("diffusion tensor", "symmetric") == "symmetric");

  // Primary discretization methods
  // if (primary == "mfd: monotone for hex") {
  //   mfd_primary_ = WhetStone::DIFFUSION_HEXAHEDRA_MONOTONE;
  // } else if (primary == "mfd: optimized for monotonicity") {
  //   mfd_primary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_MONOTONICITY;
  // } else if (primary == "mfd: two-point flux approximation") {
  //   mfd_primary_ = WhetStone::DIFFUSION_TPFA;
  // } else if (primary == "mfd: optimized for sparsity") {
  //   mfd_primary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY;
  // } else if (primary == "mfd: support operator") {
  //   mfd_primary_ = WhetStone::DIFFUSION_SUPPORT_OPERATOR;
  // } else if (primary == "mfd: default") {
  //   mfd_primary_ = WhetStone::DIFFUSION_POLYHEDRA_SCALED;
  // } else {
  //   Errors::Message msg;
  //   msg << "PDE_DiffusionMFD: primary discretization method \"" << primary << "\" is not supported.";
  //   Exceptions::amanzi_throw(msg);
  // }

  std::cout<<"Prim = "<<primary<<" Sec = "<<secondary<<std::endl;
  if (primary != "mfd: two-point flux approximation") {
    AMANZI_ASSERT(false);
  } else {
    mfd_primary_ = WhetStone::DIFFUSION_TPFA;
    mfd_secondary_ = WhetStone::DIFFUSION_TPFA;
  }
    
  
  // // Secondary discretization methods
  // if (secondary == "mfd: two-point flux approximation") {
  //   mfd_secondary_ = WhetStone::DIFFUSION_TPFA;
  // } else if (secondary == "mfd: optimized for sparsity") {
  //   mfd_secondary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY;
  // } else if (secondary == "mfd: support operator") {
  //   mfd_secondary_ = WhetStone::DIFFUSION_SUPPORT_OPERATOR;
  // } else if (secondary == "mfd: default") {
  //   mfd_secondary_ = WhetStone::DIFFUSION_POLYHEDRA_SCALED;
  // } else {
  //   Errors::Message msg;
  //   msg << "PDE_DiffusionMFD: secondary discretization method \"" << secondary << "\" is not supported.";
  //   Exceptions::amanzi_throw(msg);
  // }

  // Define stencil for the MFD diffusion method.
  std::vector<std::string> names;
  if (plist_.isParameter("schema")) {
    names = plist_.get<Teuchos::Array<std::string> > ("schema").toVector();
  } else {
    names.resize(2);
    names[0] = "face";
    names[1] = "cell";
    plist_.set<Teuchos::Array<std::string> >("schema", names);
  }

  int schema_dofs = 0;
  for (int i = 0; i < names.size(); i++) {
    if (names[i] == "cell") {
      schema_dofs += OPERATOR_SCHEMA_DOFS_CELL;
    } else if (names[i] == "node") {
      schema_dofs += OPERATOR_SCHEMA_DOFS_NODE;
    } else if (names[i] == "face") {
      schema_dofs += OPERATOR_SCHEMA_DOFS_FACE;
    }
  }

  if (schema_dofs == OPERATOR_SCHEMA_DOFS_NODE) {
    local_op_schema_ = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_NODE;
  } else if (schema_dofs == (OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL)) {
    local_op_schema_ = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL;
  } else if (schema_dofs == (OPERATOR_SCHEMA_DOFS_CELL)) {
    local_op_schema_ = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL;
  } else {
    Errors::Message msg;
    msg << "PDE_DiffusionMFD: \"schema\" must be CELL, FACE+CELL, or NODE";
    Exceptions::amanzi_throw(msg);
  }

  // define stencil for the assembled matrix
  schema_prec_dofs_ = 0;
  if (plist_.isParameter("preconditioner schema")) {
    names = plist_.get<Teuchos::Array<std::string> > ("preconditioner schema").toVector();
    for (int i = 0; i < names.size(); i++) {
      std::cout<<"Name "<<i<<" = "<<names[i]<<std::endl;
      if (names[i] == "cell") {
        schema_prec_dofs_ += OPERATOR_SCHEMA_DOFS_CELL;
      } else if (names[i] == "node") {
        schema_prec_dofs_ += OPERATOR_SCHEMA_DOFS_NODE;
      } else if (names[i] == "face") {
        schema_prec_dofs_ += OPERATOR_SCHEMA_DOFS_FACE;
      }
    } 
  } else {
    schema_prec_dofs_ = schema_dofs;
  }
}


/* ******************************************************************
 * Virtual constructor, sets up global op, cvs.
 ****************************************************************** */
void PDE_DiffusionMFD::Init()
{
  std::cout<<"PDE_DiffusionMFD::Init"<<std::endl;

  // create or check the existing Operator
  int global_op_schema = schema_prec_dofs_;  
  if (global_op_ == Teuchos::null) {
    global_op_schema_ = global_op_schema;

    // build the CVS from the global schema
    auto cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);

    if (global_op_schema & OPERATOR_SCHEMA_DOFS_CELL){
      std::cout<<"Adding cell component"<<std::endl;
      cvs->AddComponent("cell", AmanziMesh::CELL, 1);
    }
    if (global_op_schema & OPERATOR_SCHEMA_DOFS_FACE){
      std::cout<<"Adding face component"<<std::endl;
      cvs->AddComponent("face", AmanziMesh::FACE, 1);
    }
    if (global_op_schema & OPERATOR_SCHEMA_DOFS_NODE)
      cvs->AddComponent("node", AmanziMesh::NODE, 1);

    // choose the Operator from the prec schema
    if (schema_prec_dofs_ == OPERATOR_SCHEMA_DOFS_NODE) {
      assert(false); 
      //global_op_ = Teuchos::rcp(new Operator_Node(cvs, plist_));
    } 
    else if (schema_prec_dofs_ == OPERATOR_SCHEMA_DOFS_CELL) {
      assert(false); 
      // cvs->AddComponent("face", AmanziMesh::FACE, 1);
      // global_op_ = Teuchos::rcp(new Operator_FaceCellScc(cvs, plist_));
      global_op_ = Teuchos::rcp(new Operator_Cell(cvs->CreateSpace(), plist_, schema_prec_dofs_));
    } 
    else if (schema_prec_dofs_ == OPERATOR_SCHEMA_DOFS_FACE) {
      assert(false); 
      //cvs->AddComponent("cell", AmanziMesh::CELL, 1);
      //global_op_ = Teuchos::rcp(new Operator_FaceCellSff(cvs, plist_));
    } 
    else if (schema_prec_dofs_ == (OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
      global_op_ = Teuchos::rcp(new Operator_FaceCell(cvs->CreateSpace(), plist_));
    } 
    else {
      Errors::Message msg;
      msg << "PDE_DiffusionMFD: \"preconditioner schema\" must be NODE, CELL, FACE, or FACE+CELL";
      Exceptions::amanzi_throw(msg);
    }

  } else {
    assert(false); 
    // constructor was given an Operator
    //global_op_schema_ = global_op_->schema();
    //mesh_ = global_op_->DomainMap().Mesh();
  }

  // Do we need to exclude the primary terms?
  exclude_primary_terms_ = plist_.get<bool>("exclude primary terms", false);
  
  // create the local Op and register it with the global Operator
  if (!exclude_primary_terms_) {
    if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_NODE)) {
      assert(false); 
      //std::string name = "Diffusion: CELL_NODE";
      //local_op_ = Teuchos::rcp(new Op_Cell_Node(name, mesh_));
    } 
    else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL |
                                  OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL)) {
      std::string name = "Diffusion: CELL_FACE+CELL";
      local_op_ = Teuchos::rcp(new Op_Cell_FaceCell(name, mesh_));
    } 
    else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL)) {
      if (plist_.get<bool>("surface operator", false)) {
        assert(false); 
        //std::string name = "Diffusion: FACE_CELL Surface";
        //local_op_ = Teuchos::rcp(new Op_SurfaceFace_SurfaceCell(name, mesh_));
      } else {
        assert(false); 
        //std::string name = "Diffusion: FACE_CELL";
        //local_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
      }
    }
    else {
      AMANZI_ASSERT(0);
    }
    global_op_->OpPushBack(local_op_);
  }
  
  // scaled constraint -- enables zero value of k on a face
  scaled_constraint_ = plist_.get<bool>("scaled constraint equation", false);
  scaled_constraint_cutoff_ = plist_.get<double>("constraint equation scaling cutoff", 1.0);
  scaled_constraint_fuzzy_ = plist_.get<double>("constraint equation fuzzy number", 1.0e-12);

  // little-k options
  AMANZI_ASSERT(!plist_.isParameter("upwind method"));
  std::string name = plist_.get<std::string>("nonlinear coefficient", "none");
  if (name == "none") {
    little_k_type_ = OPERATOR_LITTLE_K_NONE;
  } else if (name == "upwind: face") {
    little_k_type_ = OPERATOR_LITTLE_K_UPWIND;  // upwind scheme (non-symmetric in general)
  } else if (name == "divk: face") {
    little_k_type_ = OPERATOR_LITTLE_K_DIVK_BASE;  // new SPD upwind scheme
  } else if (name == "divk: cell-face") {
    little_k_type_ = OPERATOR_LITTLE_K_DIVK;  // standard SPD upwind scheme
  } else if (name == "standard: cell") {
    little_k_type_ = OPERATOR_LITTLE_K_STANDARD;  // cell-centered scheme.
  } else if (name == "divk: cell-grad-face-twin") {  
    little_k_type_ = OPERATOR_LITTLE_K_DIVK_TWIN_GRAD;
  } else if (name == "divk: cell-face-twin") {  
    little_k_type_ = OPERATOR_LITTLE_K_DIVK_TWIN;  // for resolved simulation
  } else {
    AMANZI_ASSERT(false);
  }

  // verify input consistency
  if (scaled_constraint_) {
    AMANZI_ASSERT(little_k_type_ != OPERATOR_LITTLE_K_DIVK &&
           little_k_type_ != OPERATOR_LITTLE_K_DIVK_TWIN);
  }

  // Do we need to calculate Newton correction terms?
  std::string jacobian = plist_.get<std::string>("Newton correction", "none");
  if (jacobian == "none") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;
  } else if (jacobian == "true Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_TRUE;
    Errors::Message msg("PDE_DiffusionMFD: \"true Jacobian\" not supported -- maybe you mean \"approximate Jacobian\"?");
    Exceptions::amanzi_throw(msg);
  } else if (jacobian == "approximate Jacobian") {
    // cannot do jacobian terms without cells
    if (!(schema_prec_dofs_ & OPERATOR_SCHEMA_DOFS_CELL)) {
      Errors::Message msg("PDE_DiffusionMFD: incompatible options.  \"approximate Jacobian\" terms require CELL quantities, and the requested preconditioner schema does not include CELL.");
      Exceptions::amanzi_throw(msg);
    }
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE;
    // create a local op
    jac_op_schema_ = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL;
    std::string opname("Jacobian FACE_CELL");
    jac_op_ = Teuchos::rcp(new Op_Face_Cell(opname, mesh_));
    global_op_->OpPushBack(jac_op_);
  } else {
    Errors::Message msg;
    msg << "PDE_DiffusionMFD: invalid parameter \"" << jacobian 
        << "\" for option \"Newton correction\" -- valid are: \"none\", \"approximate Jacobian\"";
    Exceptions::amanzi_throw(msg);
  }

  // miscalleneous variables
  mass_matrices_initialized_ = false;
  K_ = Teuchos::null;
  k_ = Teuchos::null;
  dkdp_ = Teuchos::null;
}


/* ******************************************************************
* Given a set of cell values, update faces using the consistency 
* equations:
*   x_f = Aff^-1 * (y_f - Afc * x_c)
****************************************************************** */
int PDE_DiffusionMFD::UpdateConsistentFaces(CompositeVector& u)
{
  assert(false); 
  //if (consistent_face_op_ == Teuchos::null) {
    // create the op
  //  Teuchos::RCP<CompositeVectorSpace> cface_cvs = Teuchos::rcp(new CompositeVectorSpace());
  //  cface_cvs->SetMesh(mesh_)->SetGhosted()
  //      ->AddComponent("face", AmanziMesh::FACE, 1);

  //  consistent_face_op_ = Teuchos::rcp(new Operator_ConsistentFace(cface_cvs, plist_.sublist("consistent faces")));
  //  consistent_face_op_->OpPushBack(local_op_);
  //  consistent_face_op_->SymbolicAssembleMatrix();
  //  consistent_face_op_->InitializePreconditioner(plist_.sublist("consistent faces").sublist("preconditioner"));
  //}

  // calculate the rhs, given by y_f - Afc * x_c
  //CompositeVector& y = *consistent_face_op_->rhs();
  //const auto y_f = y.ViewComponent("face", true);
  //y_f = global_op_->rhs()->ViewComponent("face", true);
  //consistent_face_op_->rhs()->PutScalarGhosted(0.0);

  // y_f - Afc * x_c
  //const auto x_c = u.ViewComponent("cell", false);
  //AmanziMesh::Entity_ID_View faces;
  //for (int c=0; c!=ncells_owned; ++c) {
  //  mesh_->cell_get_faces(c, faces);
  //  int nfaces = faces.size();

  //  WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

  //  for (int n=0; n!=nfaces; ++n) {
  //    y_f(0,faces[n]) -= Acell(n,nfaces) * x_c(0,c);
  //  }
  //}

  //y.GatherGhostedToMaster("face", Add);

  // x_f = Aff^-1 * ...
  //consistent_face_op_->AssembleMatrix();
  //consistent_face_op_->UpdatePreconditioner();

  //int ierr = 0;
  //if (plist_.sublist("consistent faces").isSublist("linear solver")) {
  //  AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> fac;
  //  Teuchos::RCP<Operator> lin_solver = fac.Create(
  //      plist_.sublist("consistent faces").sublist("linear solver"), consistent_face_op_);

  //  CompositeVector u_f_copy(y);
  //  ierr = lin_solver->ApplyInverse(y, u_f_copy);
  //  u.ViewComponent("face", false) = u_f_copy.ViewComponent("face", false);
  //} else {
  //  CompositeVector u_f_copy(y);
  //  ierr = consistent_face_op_->ApplyInverse(y, u);
  //  u.ViewComponent("face", false) = u_f_copy.ViewComponent("face", false);
  //}
  
  //return (ierr > 0) ? 0 : 1; 
  return 0; 
}
  

/* ******************************************************************
* Calculates transmissibility value on the given BOUNDARY face f.
****************************************************************** */
double PDE_DiffusionMFD::ComputeTransmissibility(int f) const
{
  std::cout<<"PDE_DiffusionMFD::ComputeTransmissibility"<<std::endl;

  WhetStone::MFD3D_Diffusion mfd(mesh_);

  AmanziMesh::Entity_ID_View cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  int c = cells[0];

  if (K_.get()) {
    return mfd.Transmissibility(f, c, K_->at_host(c));
  } else {
    WhetStone::Tensor<> Kc(mesh_->space_dimension(), 1);
    Kc(0, 0) = 1.0;
    return mfd.Transmissibility(f, c, Kc);
  }
}

}  // namespace Operators
}  // namespace Amanzi
