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
#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Diffusion.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "Op.hh"
#include "Op_Cell_Node.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Face_Cell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"

#include "OperatorDefs.hh"
#include "Operator_FaceCell.hh"
#include "Operator_FaceCellScc.hh"
#include "Operator_FaceCellSff.hh"
#include "Operator_Node.hh"
#include "Operator_ConsistentFace.hh"
#include "UniqueLocalIndex.hh"

#include "PDE_DiffusionMFD.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void PDE_DiffusionMFD::SetTensorCoefficient(
    const Teuchos::RCP<const std::vector<WhetStone::Tensor> >& K)
{
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
  k_ = k;
  dkdp_ = dkdp;

  // compatibility checks
  if (k_ != Teuchos::null) {
    if (little_k_ != OPERATOR_LITTLE_K_UPWIND) {
      AMANZI_ASSERT(k->HasComponent("cell"));
    }

    if (little_k_ != OPERATOR_LITTLE_K_STANDARD &&
        little_k_ != OPERATOR_LITTLE_K_NONE) {
      AMANZI_ASSERT(k->HasComponent("face"));
    }

    if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN || 
        little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
      AMANZI_ASSERT(k->HasComponent("twin"));
    }

    if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
      AMANZI_ASSERT(k->HasComponent("grad"));
    }
  }

  // verify that mass matrices were initialized.
  if (!mass_matrices_initialized_) {
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
  if (k_ != Teuchos::null) k_->ScatterMasterToGhosted();

  if (!exclude_primary_terms_) {
    if (local_op_schema_ & OPERATOR_SCHEMA_DOFS_NODE) {
      UpdateMatricesNodal_();
    } else if ((local_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) &&
               (local_op_schema_ & OPERATOR_SCHEMA_DOFS_FACE)) {
      if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
        UpdateMatricesMixedWithGrad_(flux);
      } else if (little_k_ == OPERATOR_LITTLE_K_NONE) {
        UpdateMatricesMixed_();
      } else {
        UpdateMatricesMixed_little_k_();
      }
    } else if (local_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) {
      UpdateMatricesTPFA_();
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
  // add Newton-type corrections
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE) {
    if (global_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) {

      if (dkdp_ !=  Teuchos::null) dkdp_->ScatterMasterToGhosted();      
      AddNewtonCorrectionCell_(flux, u, scalar_factor);
      
    } else {
      Errors::Message msg("PDE_DiffusionMFD: Newton correction may only be applied to schemas that include CELL dofs.");
      Exceptions::amanzi_throw(msg);
    }
  }
}


void PDE_DiffusionMFD::UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& factor)
{
  // add Newton-type corrections
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE) {
    if (global_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) {

      if (dkdp_ !=  Teuchos::null) dkdp_->ScatterMasterToGhosted();      
      AddNewtonCorrectionCell_(flux, u, factor);

    } else {
      Errors::Message msg("PDE_DiffusionMFD: Newton correction may only be applied to schemas that include CELL dofs.");
      Exceptions::amanzi_throw(msg);
    }
  }
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

  // preparing little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_grad = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_twin = Teuchos::null;
  if (k_ != Teuchos::null) {
    k_cell = k_->ViewComponent("cell");
    k_face = k_->ViewComponent("face", true);
    k_grad = k_->ViewComponent("grad");
    if (k_->HasComponent("twin")) k_twin = k_->ViewComponent("twin", true);
  }

  // update matrix blocks
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  WhetStone::DenseMatrix Wff;

  AmanziMesh::Entity_ID_List faces, cells;

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    // mean value and gradient of nonlinear factor
    double kc = (*k_cell)[0][c];
    AmanziGeometry::Point kgrad(dim);
    for (int i = 0; i < dim; i++) kgrad[i] = (*k_grad)[i][c];
 
    // upwinded values of nonlinear factor
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    std::vector<double> kf(nfaces, 1.0); 
    if (k_twin == Teuchos::null) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    } else {
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        kf[n] = (c == cells[0]) ? (*k_face)[0][f] : (*k_twin)[0][f];
      }
    }

    if (K_.get()) Kc = (*K_)[c];
    mfd.MassMatrixInverseDivKScaled(c, Kc, kc, kgrad, Wff);

    WhetStone::DenseMatrix Acell(nfaces + 1, nfaces + 1);

    double matsum = 0.0; 
    for (int n = 0; n < nfaces; n++) {
      double rowsum = 0.0;
      for (int m = 0; m < nfaces; m++) {
        double tmp = Wff(n, m) * kf[n] * kf[m];
        rowsum += tmp;
        Acell(n, m) = tmp;
      }

      Acell(n, nfaces) = -rowsum;
      Acell(nfaces, n) = -rowsum;
      matsum += rowsum;
    }
    Acell(nfaces, nfaces) = matsum;
    local_op_->matrices[c] = Acell;
  }
}


/* ******************************************************************
* Basic routine for each operator: creation of elemental matrices.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesMixed_()
{
  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::DenseMatrix& Wff = Wff_cells_[c];
    int nfaces = Wff.NumRows();
    WhetStone::DenseMatrix Acell(nfaces + 1, nfaces + 1);

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

    local_op_->matrices[c] = Acell;
  }
}


/* ******************************************************************
* Basic routine for each operator: creation of elemental matrices.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesMixed_little_k_()
{
  // un-rolling little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_twin = Teuchos::null;
  if (k_ != Teuchos::null) {
    if (k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");
    if (k_->HasComponent("face")) k_face = k_->ViewComponent("face", true);
    if (k_->HasComponent("twin")) k_twin = k_->ViewComponent("twin", true);
  }

  // update matrix blocks
  AmanziMesh::Entity_ID_List faces, cells;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Wff = Wff_cells_[c];
    WhetStone::DenseMatrix Acell(nfaces + 1, nfaces + 1);

    // Update terms due to nonlinear coefficient
    double kc(1.0);
    std::vector<double> kf(nfaces, 1.0); 
   
    if (k_cell != Teuchos::null && k_cell.get()) kc = (*k_cell)[0][c];

    // -- chefs recommendation: SPD discretization with upwind
    if (little_k_ == OPERATOR_LITTLE_K_DIVK && k_face != Teuchos::null) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];

    // -- new scheme: SPD discretization with upwind and equal spliting
    } else if (little_k_ == OPERATOR_LITTLE_K_DIVK_BASE) {
      kc = 1.0;
      for (int n = 0; n < nfaces; n++) kf[n] = std::sqrt((*k_face)[0][faces[n]]);

    // -- same as above but remains second-order for dicontinuous coefficients
    } else if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN) {
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        kf[n] = (c == cells[0]) ? (*k_face)[0][f] : (*k_twin)[0][f];
      }

    // -- the second most popular choice: classical upwind
    } else if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];

    } else if (little_k_ == OPERATOR_LITTLE_K_STANDARD) {
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    }

    // create stiffness matrix by ellimination of the mass matrix
    // -- all methods expect for DIVK-family of methods.
    if ((little_k_ & OPERATOR_LITTLE_K_DIVK_BASE) == 0) {
      // -- not scaled constraint: kr > 0
      if (!scaled_constraint_) {
        double matsum = 0.0; 
        for (int n = 0; n < nfaces; n++) {
          double rowsum = 0.0;
          for (int m = 0; m < nfaces; m++) {
            double tmp = Wff(n, m) * kf[n];
            rowsum += tmp;
            Acell(n, m) = tmp;
          }

          Acell(n, nfaces) = -rowsum;
          matsum += rowsum;
        }
        Acell(nfaces, nfaces) = matsum;

        for (int n = 0; n < nfaces; n++) {
          double colsum = 0.0;
          for (int m = 0; m < nfaces; m++) colsum += Acell(m, n);
          Acell(nfaces, n) = -colsum;
        }

      // -- scaled constraint: kr >= 0
      } else {
        double matsum = 0.0;
        for (int n = 0; n < nfaces; n++) {
          double rowsum = 0.0;
          double cur_kf = (kf[n] < scaled_constraint_cutoff_) ? 1.0 : kf[n];
          for (int m = 0; m < nfaces; m++) {
            double tmp = Wff(n, m) * cur_kf;
            rowsum += tmp;
            Acell(n, m) = tmp;
          }
          Acell(n, nfaces) = -rowsum;
        }

        for (int n = 0; n < nfaces; n++) {
          double colsum = 0.0;
          for (int m = 0; m < nfaces; m++) colsum += Wff(m, n) * kf[m];
          Acell(nfaces, n) = -colsum;
          matsum += colsum;
        }
        Acell(nfaces, nfaces) = matsum;
      }
    }

    // Amanzi's first upwind: the family of DIVK fmethods
    if (little_k_ & OPERATOR_LITTLE_K_DIVK_BASE) {
      double matsum = 0.0; 
      for (int n = 0; n < nfaces; n++) {
        double rowsum = 0.0;
        for (int m = 0; m < nfaces; m++) {
          double tmp = Wff(n, m) * kf[n] * kf[m] / kc;
          rowsum += tmp;
          Acell(n, m) = tmp;
        }
          
        Acell(n, nfaces) = -rowsum;
        Acell(nfaces, n) = -rowsum;
        matsum += rowsum;
      }

      Acell(nfaces, nfaces) = matsum;
    }
    
    local_op_->matrices[c] = Acell;
  }
}


/* ******************************************************************
* Calculate elemental stiffness matrices: nodal DOFs.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesNodal_()
{
  AMANZI_ASSERT(!scaled_constraint_);

  // update matrix blocks
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  mfd.ModifyStabilityScalingFactor(factor_);

  AmanziMesh::Entity_ID_List nodes;

  nfailed_primary_ = 0;

  WhetStone::Tensor K(2, 1);
  K(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    if (K_.get()) K = (*K_)[c];

    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    WhetStone::DenseMatrix Acell(nnodes, nnodes);

    int method = mfd_primary_;
    int ok = WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED;

    if (method == WhetStone::DIFFUSION_OPTIMIZED_FOR_MONOTONICITY) {
      ok = mfd.StiffnessMatrixMMatrix(c, K, Acell);
      method = mfd_secondary_;
    } else {
      ok = mfd.StiffnessMatrix(c, K, Acell);
      method = mfd_secondary_;
    }

    if (ok != WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) {
      nfailed_primary_++;
      ok = mfd.StiffnessMatrix(c, K, Acell);
    }

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("Stiffness_MFD: unexpected failure of LAPACK in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }

    local_op_->matrices[c] = Acell;
  }
}


/* ******************************************************************
* Calculate and assemble fluxes using the TPFA scheme.
* This routine does not use little k.
****************************************************************** */
void PDE_DiffusionMFD::UpdateMatricesTPFA_()
{
  // populate transmissibilities
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  WhetStone::DenseMatrix Mff;

  CompositeVectorSpace cv_space;
  cv_space.SetMesh(mesh_);
  cv_space.SetGhosted(true);
  cv_space.SetComponent("face", AmanziMesh::FACE, 1);

  Teuchos::RCP<CompositeVector> T = Teuchos::RCP<CompositeVector>(new CompositeVector(cv_space, true));
  Epetra_MultiVector& Ttmp = *T->ViewComponent("face", true);

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;

  AmanziMesh::Entity_ID_List cells, faces;
  Ttmp.PutScalar(0.0);

  for (int c = 0; c < ncells_owned; c++) {
    if (K_.get()) Kc = (*K_)[c];
    if (Kc.isZero()) continue;  // We skip zero matrices

    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    mfd.MassMatrixInverseTPFA(c, Kc, Mff);
   
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      Ttmp[0][f] += 1.0 / Mff(n, n);
    }
  }
  T->GatherGhostedToMaster();
 
  // populate the global matrix
  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);

    if (Ttmp[0][f] == 0.0) {
      Aface = 0.0;
      local_op_->matrices[f] = Aface;
      continue;  // We skip zero transmissibilities
    }

    if (ncells == 2) {
      double coef = 1.0 / Ttmp[0][f];
      Aface(0, 0) =  coef;
      Aface(1, 1) =  coef;
      Aface(0, 1) = -coef;
      Aface(1, 0) = -coef;
    } else {
      double coef = 1.0 / Ttmp[0][f];
      Aface(0, 0) = coef;
    }

    local_op_->matrices[f] = Aface;
  }
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
  if (!exclude_primary_terms_) {
    if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL
                           | OPERATOR_SCHEMA_DOFS_FACE
                           | OPERATOR_SCHEMA_DOFS_CELL)) {
      AMANZI_ASSERT(bcs_trial_.size() == 1);
      AMANZI_ASSERT(bcs_test_.size() == 1);
      ApplyBCs_Mixed_(bcs_trial_[0].ptr(), bcs_test_[0].ptr(), primary, eliminate, essential_eqn);
    
    } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_FACE
                                  | OPERATOR_SCHEMA_DOFS_CELL)) {
      AMANZI_ASSERT(bcs_trial_.size() == 1);
      AMANZI_ASSERT(bcs_test_.size() == 1);
      ApplyBCs_Cell_(bcs_trial_[0].ptr(), bcs_test_[0].ptr(), primary, eliminate, essential_eqn);
    
    } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL
                                  | OPERATOR_SCHEMA_DOFS_NODE)) {
      Teuchos::Ptr<const BCs> bc_f, bc_n;
      for (const auto& bc : bcs_trial_) {
        if (bc->kind() == AmanziMesh::FACE) {
          bc_f = bc.ptr();
        } else if (bc->kind() == AmanziMesh::NODE) {
          bc_n = bc.ptr();
        }
      }
      ApplyBCs_Nodal_(bc_f.ptr(), bc_n.ptr(), primary, eliminate, essential_eqn);
    }
  }

  if (jac_op_ != Teuchos::null) {
    const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
    AMANZI_ASSERT(bc_model.size() == nfaces_wghost);

    for (int f = 0; f != nfaces_owned; ++f) {
      WhetStone::DenseMatrix& Aface = jac_op_->matrices[f];

      if (bc_model[f] == OPERATOR_BC_NEUMANN ||
          bc_model[f] == OPERATOR_BC_TOTAL_FLUX) {
        jac_op_->matrices_shadow[f] = Aface;
        Aface *= 0.0;
      }
    }
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
  // apply diffusion type BCs to FACE-CELL system
  AmanziMesh::Entity_ID_List faces;

  const std::vector<int>& bc_model_trial = bc_trial->bc_model();
  const std::vector<int>& bc_model_test = bc_test->bc_model();

  const std::vector<double>& bc_value = bc_trial->bc_value();
  const std::vector<double>& bc_mixed = bc_trial->bc_mixed();

  AMANZI_ASSERT(bc_model_trial.size() == nfaces_wghost);
  AMANZI_ASSERT(bc_value.size() == nfaces_wghost);

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");

  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) {
    if (k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");
    if (k_->HasComponent("face")) k_face = k_->ViewComponent("face", true);
  }

  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    
    // Update terms due to nonlinear coefficient
    double kc(1.0);
    std::vector<double> kf(nfaces, 1.0);
    if (scaled_constraint_) {
      // un-rolling little-k data
      if (k_cell != Teuchos::null && k_cell.get()) kc = (*k_cell)[0][c];
      
      if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
        for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
        
      } else if (little_k_ == OPERATOR_LITTLE_K_STANDARD) {
        for (int n = 0; n < nfaces; n++) kf[n] = kc;
      }
    }
    
    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
        
    // essential conditions for test functions
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      if (bc_model_test[f] == OPERATOR_BC_DIRICHLET) {
        if (flag) {  // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        for (int m = 0; m < nfaces + 1; m++) Acell(n, m) = 0.0;
      }
    }

    // conditions for trial functions
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      double value = bc_value[f];

      if (bc_model_trial[f] == OPERATOR_BC_DIRICHLET) {
        // make a copy of elemental matrix for post-processing
        if (flag) {
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }

        if (eliminate) { 
          for (int m = 0; m < nfaces; m++) {
            rhs_face[0][faces[m]] -= Acell(m, n) * value;
            Acell(m, n) = 0.0;
          }

          rhs_cell[0][c] -= Acell(nfaces, n) * value;
          Acell(nfaces, n) = 0.0;
        }

        if (essential_eqn) {
          rhs_face[0][f] = value;
          Acell(n, n) = 1.0;
        }

      } else if (bc_model_trial[f] == OPERATOR_BC_NEUMANN && primary) {
        if (scaled_constraint_) {
          if (std::abs(kf[n]) < scaled_constraint_fuzzy_) {
            AMANZI_ASSERT(value == 0.0);
            rhs_face[0][f] = 0.0;
          } else if (kf[n] < scaled_constraint_cutoff_) {
            rhs_face[0][f] -= value * mesh_->face_area(f) / kf[n];
          } else {
            rhs_face[0][f] -= value * mesh_->face_area(f);
          }
        } else {
          rhs_face[0][f] -= value * mesh_->face_area(f);
        }

      } else if (bc_model_trial[f] == OPERATOR_BC_TOTAL_FLUX && primary) {
        if (scaled_constraint_ && kf[n] < scaled_constraint_cutoff_) {
          AMANZI_ASSERT(false);
        } else {
          rhs_face[0][f] -= value * mesh_->face_area(f);
        }

      } else if (bc_model_trial[f] == OPERATOR_BC_MIXED && primary) {
        if (flag) {  // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        double area = mesh_->face_area(f);
        if (scaled_constraint_) {
          if (std::abs(kf[n]) < scaled_constraint_fuzzy_) {
            AMANZI_ASSERT((value == 0.0) && (bc_mixed[f] == 0.0));
            rhs_face[0][f] = 0.0;
          } else if (kf[n] < scaled_constraint_cutoff_) {
            rhs_face[0][f] -= value * area / kf[n];
            Acell(n, n) += bc_mixed[f] * area / kf[n];
          } else {
            rhs_face[0][f] -= value * area;
            Acell(n, n) += bc_mixed[f] * area;
          }
        } else {
          rhs_face[0][f] -= value * area;
          Acell(n, n) += bc_mixed[f] * area;
        }
      }
    }
  }

  global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* Apply BCs on cell operators
****************************************************************** */
void PDE_DiffusionMFD::ApplyBCs_Cell_(
   const Teuchos::Ptr<const BCs>& bc_trial,
   const Teuchos::Ptr<const BCs>& bc_test,
   bool primary, bool eliminate, bool essential_eqn)
{
  // apply diffusion type BCs to CELL system
  AmanziMesh::Entity_ID_List cells;

  const std::vector<int>& bc_model = bc_trial->bc_model();
  const std::vector<double>& bc_value = bc_trial->bc_value();
  const std::vector<double>& bc_mixed = bc_trial->bc_mixed();

  AMANZI_ASSERT(bc_model.size() == nfaces_wghost);
  AMANZI_ASSERT(bc_value.size() == nfaces_wghost);

  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");
    
  for (int f = 0; f != nfaces_owned; ++f) {
    WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
      
    if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      rhs_cell[0][cells[0]] += bc_value[f] * Aface(0, 0);
    }
    // Neumann condition contributes to the RHS
    else if (bc_model[f] == OPERATOR_BC_NEUMANN && primary) {
      local_op_->matrices_shadow[f] = Aface;
      
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      rhs_cell[0][cells[0]] -= bc_value[f] * mesh_->face_area(f);
      Aface *= 0.0;
    }
    // solve system of two equations in three unknowns
    else if (bc_model[f] == OPERATOR_BC_MIXED && primary) {
      local_op_->matrices_shadow[f] = Aface;
      
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      double area = mesh_->face_area(f);
      double factor = area / (1.0 + bc_mixed[f] * area / Aface(0, 0));
      rhs_cell[0][cells[0]] -= bc_value[f] * factor;
      Aface(0, 0) = bc_mixed[f] * factor;
    }
  }
}


/* ******************************************************************
* Apply BCs on nodal operators
****************************************************************** */
void PDE_DiffusionMFD::ApplyBCs_Nodal_(
    const Teuchos::Ptr<const BCs>& bc_f,
    const Teuchos::Ptr<const BCs>& bc_v,
    bool primary, bool eliminate, bool essential_eqn)
{
  AmanziMesh::Entity_ID_List faces, nodes, cells;

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_node = *global_op_->rhs()->ViewComponent("node", true);

  int nn(0), nm(0);
  for (int c = 0; c != ncells_owned; ++c) {
    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

    // process boundary integrals
    if (bc_f != Teuchos::null) {
      const std::vector<int>& bc_model = bc_f->bc_model();
      const std::vector<double>& bc_value = bc_f->bc_value();
      const std::vector<double>& bc_mixed = bc_f->bc_mixed();

      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];

        if (bc_model[f] == OPERATOR_BC_NEUMANN && primary) {
          nn++;
          double value = bc_value[f];
          double area = mesh_->face_area(f);

          mesh_->face_get_nodes(f, &nodes);
          int nnodes = nodes.size();

          for (int m = 0; m < nnodes; m++) {
            int v = nodes[m];
            if (bc_v->bc_model()[v] != OPERATOR_BC_DIRICHLET)
              rhs_node[0][v] -= value * area / nnodes;
          }
        } else if (bc_model[f] == OPERATOR_BC_MIXED && primary) {
          nm++;
          if (flag) {  // make a copy of cell-based matrix
            local_op_->matrices_shadow[c] = Acell;
            flag = false;
          }
          double value = bc_value[f];
          double area = mesh_->face_area(f);

          mesh_->face_get_nodes(f, &nodes);
          int nnodes = nodes.size();

          for (int m = 0; m < nnodes; m++) {
            int v = nodes[m];
            if (bc_v->bc_model()[v] != OPERATOR_BC_DIRICHLET)
              rhs_node[0][v] -= value * area / nnodes;
            Acell(n, n) += bc_mixed[f] * area / nnodes;
          }
        }
      }
    } 

    if (bc_v != Teuchos::null) {
      const std::vector<int>& bc_model = bc_v->bc_model();
      const std::vector<double>& bc_value = bc_v->bc_value();

      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      // essential conditions for test functions
      for (int n = 0; n != nnodes; ++n) {
        int v = nodes[n];
        if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
          if (flag) {  // make a copy of elemental matrix
            local_op_->matrices_shadow[c] = Acell;
            flag = false;
          }
          for (int m = 0; m < nnodes; m++) Acell(n, m) = 0.0;
        }
      }

      for (int n = 0; n != nnodes; ++n) {
        int v = nodes[n];
        double value = bc_value[v];

        if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
          if (flag) {  // make a copy of cell-based matrix
            local_op_->matrices_shadow[c] = Acell;
            flag = false;
          }
     
          if (eliminate) {
            for (int m = 0; m < nnodes; m++) {
              rhs_node[0][nodes[m]] -= Acell(m, n) * value;
              Acell(m, n) = 0.0;
            }
          }

          // We take into account multiple contributions to matrix diagonal
          // by dividing by the number of cells attached to a vertex.
          if (essential_eqn) {
            mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &cells);
            if (v < nnodes_owned) rhs_node[0][v] = value;
            Acell(n, n) = 1.0 / cells.size();
          }
        }
      }
    }
  } 

  global_op_->rhs()->GatherGhostedToMaster("node", Add);
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
  // hack: ignore correction if no flux provided.
  if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

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

    // prototype for future limiters (external or internal ?)
    vmod *= scalar_factor;

    // define the upwind cell, index i in this case
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


void PDE_DiffusionMFD::AddNewtonCorrectionCell_(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& factor)
{
  // hack: ignore correction if no flux provided.
  if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

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

    double scalar_factor=0.;
    for (int j=0; j<ncells; j++) scalar_factor += factor_cell[0][cells[j]];
    scalar_factor *= 1./ncells;

    // prototype for future limiters (external or internal ?)
    vmod *= scalar_factor;

    // define the upwind cell, index i in this case
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
* Given pressures, reduce the problem to Lagrange multipliers.
****************************************************************** */
void PDE_DiffusionMFD::ModifyMatrices(const CompositeVector& u)
{
  if (local_op_schema_ != (OPERATOR_SCHEMA_BASE_CELL |
                           OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
    std::cout << "Schema " << global_op_schema_ << " is not supported" << std::endl;
    AMANZI_ASSERT(0);
  }

  // populate the matrix
  AmanziMesh::Entity_ID_List faces;
  const Epetra_MultiVector& u_c = *u.ViewComponent("cell");

  global_op_->rhs()->PutScalarGhosted(0.0);

  Epetra_MultiVector& rhs_f = *global_op_->rhs()->ViewComponent("face", true);
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      rhs_f[0][f] -= Acell(n, nfaces) * u_c[0][c];
      Acell(n, nfaces) = 0.0;
      Acell(nfaces, n) = 0.0;
    }
  }

  // Assemble all right-hand sides
  global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* WARNING: Since diffusive flux may be discontinuous (e.g. for
* Richards equation), we derive it in exactly the same manner as 
* in gravity routines.
* **************************************************************** */
void PDE_DiffusionMFD::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                  const Teuchos::Ptr<CompositeVector>& flux)
{
  // Initialize intensity in ghost faces.
  flux->PutScalar(0.0);
  u->ScatterMasterToGhosted("face");

  if (k_ != Teuchos::null) {
    if (k_->HasComponent("face")) k_->ScatterMasterToGhosted("face");
  }

  const Epetra_MultiVector& u_cell = *u->ViewComponent("cell");
  const Epetra_MultiVector& u_face = *u->ViewComponent("face", true);
  Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<int> hits(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    WhetStone::DenseVector v(nfaces + 1), av(nfaces + 1);
    for (int n = 0; n < nfaces; n++) {
      v(n) = u_face[0][faces[n]];
    }
    v(nfaces) = u_cell[0][c];

    if (local_op_->matrices_shadow[c].NumRows() == 0) { 
      local_op_->matrices[c].Multiply(v, av, false);
    } else {
      local_op_->matrices_shadow[c].Multiply(v, av, false);
    }

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (f < nfaces_owned) {
        flux_data[0][f] -= av(n) * dirs[n];
        hits[f]++;
      }
    }
  }

  for (int f = 0; f != nfaces_owned; ++f) {
    flux_data[0][f] /= hits[f];
  }
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
  // Initialize intensity in ghost faces.
  u->ScatterMasterToGhosted("face");

  const Epetra_MultiVector& u_cell = *u->ViewComponent("cell");
  const Epetra_MultiVector& u_face = *u->ViewComponent("face", true);
  Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);

  flux_data.PutScalar(0.0);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  const auto& fmap = *flux->Map().Map("face", true);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    WhetStone::DenseVector v(nfaces + 1), av(nfaces + 1);
    for (int n = 0; n < nfaces; n++) {
      v(n) = u_face[0][faces[n]];
    }
    v(nfaces) = u_cell[0][c];

    if (local_op_->matrices_shadow[c].NumRows() == 0) { 
      local_op_->matrices[c].Multiply(v, av, false);
    } else {
      local_op_->matrices_shadow[c].Multiply(v, av, false);
    }

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      int g = fmap.FirstPointInElement(f);

      int ndofs = fmap.ElementSize(f);
      if (ndofs > 1) g += Operators::UniqueIndexFaceToCells(*mesh_, f, c);

      flux_data[0][g] -= av(n) * dirs[n];
    }
  }

  flux->GatherGhostedToMaster(Add);
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.
****************************************************************** */
void PDE_DiffusionMFD::CreateMassMatrices_()
{
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  WhetStone::DenseMatrix Wff;
  bool surface_mesh = (mesh_->manifold_dimension() != mesh_->space_dimension());

  mfd.ModifyStabilityScalingFactor(factor_);

  Wff_cells_.resize(ncells_owned);

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;

  for (int c = 0; c < ncells_owned; c++) {
    int ok;
    if (K_.get()) Kc = (*K_)[c];

    // For problems with degenerate coefficients we should skip WhetStone.
    if (Kc.Trace() == 0.0) {
      int nfaces = mesh_->cell_get_num_faces(c);
      Wff.Reshape(nfaces, nfaces);
      Wff.PutScalar(0.0);
      ok = WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK;
    } else if (surface_mesh) {
      ok = mfd.MassMatrixInverseSurface(c, Kc, Wff);
    } else {
      int method = mfd_primary_;
      ok = WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED;

      // try primary and then secondary discretization methods.
      if (method == WhetStone::DIFFUSION_HEXAHEDRA_MONOTONE) {
        ok = mfd.MassMatrixInverseMMatrixHex(c, Kc, Wff);
        method = mfd_secondary_;
      } else if (method == WhetStone::DIFFUSION_OPTIMIZED_FOR_MONOTONICITY) {
        ok = mfd.MassMatrixInverseMMatrix(c, Kc, Wff);
        method = mfd_secondary_;
      }

      if (ok != WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) {
        if (method == WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY) {
          ok = mfd.MassMatrixInverseOptimized(c, Kc, Wff);
        } else if (method == WhetStone::DIFFUSION_TPFA) {
          ok = mfd.MassMatrixInverseTPFA(c, Kc, Wff);
        } else if (method == WhetStone::DIFFUSION_SUPPORT_OPERATOR) {
          ok = mfd.MassMatrixInverseSO(c, Kc, Wff);
        } else if (method == WhetStone::DIFFUSION_POLYHEDRA_SCALED) {
          if (K_symmetric_) {
            ok = mfd.MassMatrixInverse(c, Kc, Wff);
          } else {
            ok = mfd.MassMatrixInverseNonSymmetric(c, Kc, Wff);
          }
        }
      }
    }

    Wff_cells_[c] = Wff;

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("PDE_DiffusionMFD: unexpected failure in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
  }

  mass_matrices_initialized_ = true;
}


/* ******************************************************************
* Scale elemental inverse mass matrices. Use case is saturated flow.
****************************************************************** */
void PDE_DiffusionMFD::ScaleMassMatrices(double s)
{
  for (int c = 0; c < ncells_owned; c++) {
    Wff_cells_[c] *= s;
  }
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void PDE_DiffusionMFD::ParsePList_(Teuchos::ParameterList& plist)
{
  // Determine discretization
  std::string primary = plist.get<std::string>("discretization primary");
  std::string secondary = plist.get<std::string>("discretization secondary", primary);
  K_symmetric_ = (plist.get<std::string>("diffusion tensor", "symmetric") == "symmetric");

  // Primary discretization methods
  if (primary == "mfd: monotone for hex") {
    mfd_primary_ = WhetStone::DIFFUSION_HEXAHEDRA_MONOTONE;
  } else if (primary == "mfd: optimized for monotonicity") {
    mfd_primary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_MONOTONICITY;
  } else if (primary == "mfd: two-point flux approximation") {
    mfd_primary_ = WhetStone::DIFFUSION_TPFA;
  } else if (primary == "mfd: optimized for sparsity") {
    mfd_primary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY;
  } else if (primary == "mfd: support operator") {
    mfd_primary_ = WhetStone::DIFFUSION_SUPPORT_OPERATOR;
  } else if (primary == "mfd: default") {
    mfd_primary_ = WhetStone::DIFFUSION_POLYHEDRA_SCALED;
  } else {
    Errors::Message msg;
    msg << "PDE_DiffusionMFD: primary discretization method \"" << primary << "\" is not supported.";
    Exceptions::amanzi_throw(msg);
  }

  // Secondary discretization methods
  if (secondary == "mfd: two-point flux approximation") {
    mfd_secondary_ = WhetStone::DIFFUSION_TPFA;
  } else if (secondary == "mfd: optimized for sparsity") {
    mfd_secondary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY;
  } else if (secondary == "mfd: support operator") {
    mfd_secondary_ = WhetStone::DIFFUSION_SUPPORT_OPERATOR;
  } else if (secondary == "mfd: default") {
    mfd_secondary_ = WhetStone::DIFFUSION_POLYHEDRA_SCALED;
  } else {
    Errors::Message msg;
    msg << "PDE_DiffusionMFD: secondary discretization method \"" << secondary << "\" is not supported.";
    Exceptions::amanzi_throw(msg);
  }

  // Define stencil for the MFD diffusion method.
  std::vector<std::string> names;
  if (plist.isParameter("schema")) {
    names = plist.get<Teuchos::Array<std::string> > ("schema").toVector();
  } else {
    names.resize(2);
    names[0] = "face";
    names[1] = "cell";
    plist.set<Teuchos::Array<std::string> >("schema", names);
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
  if (plist.isParameter("preconditioner schema")) {
    names = plist.get<Teuchos::Array<std::string> > ("preconditioner schema").toVector();
    for (int i = 0; i < names.size(); i++) {
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
* Put here stuff that has to be done in constructor.
****************************************************************** */
void PDE_DiffusionMFD::Init(Teuchos::ParameterList& plist)
{
  // create or check the existing Operator
  int global_op_schema = schema_prec_dofs_;  
  if (global_op_ == Teuchos::null) {
    global_op_schema_ = global_op_schema;

    // build the CVS from the global schema
    auto cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);

    if (global_op_schema & OPERATOR_SCHEMA_DOFS_CELL)
      cvs->AddComponent("cell", AmanziMesh::CELL, 1);
    if (global_op_schema & OPERATOR_SCHEMA_DOFS_FACE)
      cvs->AddComponent("face", AmanziMesh::FACE, 1);
    if (global_op_schema & OPERATOR_SCHEMA_DOFS_NODE)
      cvs->AddComponent("node", AmanziMesh::NODE, 1);

    // choose the Operator from the prec schema
    if (schema_prec_dofs_ == OPERATOR_SCHEMA_DOFS_NODE) {
      global_op_ = Teuchos::rcp(new Operator_Node(cvs, plist));
    } 
    else if (schema_prec_dofs_ == OPERATOR_SCHEMA_DOFS_CELL) {
      // cvs->AddComponent("face", AmanziMesh::FACE, 1);
      // global_op_ = Teuchos::rcp(new Operator_FaceCellScc(cvs, plist));
      global_op_ = Teuchos::rcp(new Operator_Cell(cvs, plist, schema_prec_dofs_));
    } 
    else if (schema_prec_dofs_ == OPERATOR_SCHEMA_DOFS_FACE) {
      cvs->AddComponent("cell", AmanziMesh::CELL, 1);
      global_op_ = Teuchos::rcp(new Operator_FaceCellSff(cvs, plist));
    } 
    else if (schema_prec_dofs_ == (OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
      global_op_ = Teuchos::rcp(new Operator_FaceCell(cvs, plist));
    } 
    else {
      Errors::Message msg;
      msg << "PDE_DiffusionMFD: \"preconditioner schema\" must be NODE, CELL, FACE, or FACE+CELL";
      Exceptions::amanzi_throw(msg);
    }

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // Do we need to exclude the primary terms?
  exclude_primary_terms_ = plist.get<bool>("exclude primary terms", false);
  
  // create the local Op and register it with the global Operator
  if (!exclude_primary_terms_) {
    if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_NODE)) {
      std::string name = "Diffusion: CELL_NODE";
      local_op_ = Teuchos::rcp(new Op_Cell_Node(name, mesh_));
    } 
    else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL |
                                  OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL)) {
      std::string name = "Diffusion: CELL_FACE+CELL";
      local_op_ = Teuchos::rcp(new Op_Cell_FaceCell(name, mesh_));
    } 
    else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL)) {
      if (plist.get<bool>("surface operator", false)) {
        std::string name = "Diffusion: FACE_CELL Surface";
        local_op_ = Teuchos::rcp(new Op_SurfaceFace_SurfaceCell(name, mesh_));
      } else {
        std::string name = "Diffusion: FACE_CELL";
        local_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
      }
    }
    else {
      AMANZI_ASSERT(0);
    }
    global_op_->OpPushBack(local_op_);
  }
  
  // scaled constraint -- enables zero value of k on a face
  scaled_constraint_ = plist.get<bool>("scaled constraint equation", false);
  scaled_constraint_cutoff_ = plist.get<double>("constraint equation scaling cutoff", 1.0);
  scaled_constraint_fuzzy_ = plist.get<double>("constraint equation fuzzy number", 1.0e-12);

  // little-k options
  AMANZI_ASSERT(!plist.isParameter("upwind method"));
  std::string name = plist.get<std::string>("nonlinear coefficient", "none");
  if (name == "none") {
    little_k_ = OPERATOR_LITTLE_K_NONE;
  } else if (name == "upwind: face") {
    little_k_ = OPERATOR_LITTLE_K_UPWIND;  // upwind scheme (non-symmetric in general)
  } else if (name == "divk: face") {
    little_k_ = OPERATOR_LITTLE_K_DIVK_BASE;  // new SPD upwind scheme
  } else if (name == "divk: cell-face") {
    little_k_ = OPERATOR_LITTLE_K_DIVK;  // standard SPD upwind scheme
  } else if (name == "standard: cell") {
    little_k_ = OPERATOR_LITTLE_K_STANDARD;  // cell-centered scheme.
  } else if (name == "divk: cell-grad-face-twin") {  
    little_k_ = OPERATOR_LITTLE_K_DIVK_TWIN_GRAD;
  } else if (name == "divk: cell-face-twin") {  
    little_k_ = OPERATOR_LITTLE_K_DIVK_TWIN;  // for resolved simulation
  } else {
    AMANZI_ASSERT(false);
  }

  // verify input consistency
  if (scaled_constraint_) {
    AMANZI_ASSERT(little_k_ != OPERATOR_LITTLE_K_DIVK &&
           little_k_ != OPERATOR_LITTLE_K_DIVK_TWIN);
  }

  // Do we need to calculate Newton correction terms?
  std::string jacobian = plist.get<std::string>("Newton correction", "none");
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
  if (consistent_face_op_ == Teuchos::null) {
    // create the op
    Teuchos::RCP<CompositeVectorSpace> cface_cvs = Teuchos::rcp(new CompositeVectorSpace());
    cface_cvs->SetMesh(mesh_)->SetGhosted()
        ->AddComponent("face", AmanziMesh::FACE, 1);

    consistent_face_op_ = Teuchos::rcp(new Operator_ConsistentFace(cface_cvs, plist_.sublist("consistent faces")));
    consistent_face_op_->OpPushBack(local_op_);
    consistent_face_op_->SymbolicAssembleMatrix();
    consistent_face_op_->InitializePreconditioner(plist_.sublist("consistent faces").sublist("preconditioner"));
  }

  // calculate the rhs, given by y_f - Afc * x_c
  CompositeVector& y = *consistent_face_op_->rhs();
  Epetra_MultiVector& y_f = *y.ViewComponent("face", true);
  y_f = *global_op_->rhs()->ViewComponent("face", true);
  consistent_face_op_->rhs()->PutScalarGhosted(0.0);

  // y_f - Afc * x_c
  const Epetra_MultiVector& x_c = *u.ViewComponent("cell", false);
  AmanziMesh::Entity_ID_List faces;
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

    for (int n=0; n!=nfaces; ++n) {
      y_f[0][faces[n]] -= Acell(n,nfaces) * x_c[0][c];
    }
  }

  y.GatherGhostedToMaster("face", Add);

  // x_f = Aff^-1 * ...
  consistent_face_op_->AssembleMatrix();
  consistent_face_op_->UpdatePreconditioner();

  int ierr = 0;
  if (plist_.sublist("consistent faces").isSublist("linear solver")) {
    AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> fac;
    Teuchos::RCP<Operator> lin_solver = fac.Create(
        plist_.sublist("consistent faces").sublist("linear solver"), consistent_face_op_);

    CompositeVector u_f_copy(y);
    ierr = lin_solver->ApplyInverse(y, u_f_copy);
    *u.ViewComponent("face", false) = *u_f_copy.ViewComponent("face", false);
  } else {
    CompositeVector u_f_copy(y);
    ierr = consistent_face_op_->ApplyInverse(y, u);
    *u.ViewComponent("face", false) = *u_f_copy.ViewComponent("face", false);
  }
  
  return (ierr > 0) ? 0 : 1;
}
  

/* ******************************************************************
* Calculates transmissibility value on the given BOUNDARY face f.
****************************************************************** */
double PDE_DiffusionMFD::ComputeTransmissibility(int f) const
{
  WhetStone::MFD3D_Diffusion mfd(mesh_);

  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  int c = cells[0];

  if (K_.get()) {
    return mfd.Transmissibility(f, c, (*K_)[c]);
  } else {
    WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
    Kc(0, 0) = 1.0;
    return mfd.Transmissibility(f, c, Kc);
  }
}

}  // namespace Operators
}  // namespace Amanzi
