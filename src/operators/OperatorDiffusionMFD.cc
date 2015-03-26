/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Epetra_Vector.h"

#include "errors.hh"
#include "WhetStoneDefs.hh"
#include "mfd3d_diffusion.hh"

#include "PreconditionerFactory.hh"
#include "MatrixFE.hh"
#include "SuperMap.hh"

#include "LinearOperator.hh"
#include "LinearOperatorFactory.hh"

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

#include "OperatorDiffusionMFD.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficients.
****************************************************************** */
void OperatorDiffusionMFD::Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                                 double rho, double mu)
{
  scalar_rho_mu_ = true;
  rho_ = rho;
  mu_ = mu;
  K_ = K;

  if (local_op_schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
    if (K_.get()) ASSERT(K_->size() == ncells_owned);
    CreateMassMatrices_();
  }
}


/* ******************************************************************
* Initialization of the operator, vector coefficients.
****************************************************************** */
void OperatorDiffusionMFD::Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                                 const Teuchos::RCP<const CompositeVector>& rho,
                                 const Teuchos::RCP<const CompositeVector>& mu)
{
  scalar_rho_mu_ = false;
  rho_cv_ = rho;
  mu_cv_ = mu;
  K_ = K;

  if (local_op_schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
    if (K_.get()) ASSERT(K_->size() == ncells_owned);
    CreateMassMatrices_();
  }
}


/* ******************************************************************
* Initialization of the operator.                                           
****************************************************************** */
void OperatorDiffusionMFD::Setup(const Teuchos::RCP<const CompositeVector>& k,
                                 const Teuchos::RCP<const CompositeVector>& dkdp)
{
  k_ = k;
  dkdp_ = dkdp;

  // compatibility
  if (upwind_ == OPERATOR_UPWIND_FLUX || 
      upwind_ == OPERATOR_UPWIND_AMANZI_ARTIFICIAL_DIFFUSION ||
      upwind_ == OPERATOR_UPWIND_AMANZI_DIVK) {
    ASSERT(k->HasComponent("face"));
  }

  if (upwind_ == OPERATOR_UPWIND_AMANZI_SECOND_ORDER) {
    ASSERT(k->HasComponent("face"));
    ASSERT(k->HasComponent("grad"));
  }
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void OperatorDiffusionMFD::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!exclude_primary_terms_) {
    if (local_op_schema_ & OPERATOR_SCHEMA_DOFS_NODE) {
      UpdateMatricesNodal_();
    } else if ((local_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) &&
               (local_op_schema_ & OPERATOR_SCHEMA_DOFS_FACE)) {
      if (upwind_ == OPERATOR_UPWIND_AMANZI_SECOND_ORDER) {
        UpdateMatricesMixedWithGrad_(flux);
      } else {
        UpdateMatricesMixed_(flux);
      }
    } else if (local_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) {
      UpdateMatricesTPFA_();
    }
  }

  // add Newton-type corrections
  // if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE) {
  //   if (global_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) {
  //     AddNewtonCorrectionCell_(flux, u);
  //   } else {
  //     Errors::Message msg("OperatorDiffusion: Newton Correction may only be applied to schemas that include CELL dofs.");
  //     Exceptions::amanzi_throw(msg);
  //   }
  // }
}


/* ******************************************************************
* Add stable approximation of Jacobian. It is done typically for 
* the preconditioner.
****************************************************************** */
void OperatorDiffusionMFD::UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  // add Newton-type corrections
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE) {
    if (global_op_schema_ & OPERATOR_SCHEMA_DOFS_CELL) {
      AddNewtonCorrectionCell_(flux, u);
    } else {
      Errors::Message msg("OperatorDiffusion: Newton correction may only be applied to schemas that include CELL dofs.");
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Second-order upwind. Mass matrices are recalculated.
****************************************************************** */
void OperatorDiffusionMFD::UpdateMatricesMixedWithGrad_(
    const Teuchos::Ptr<const CompositeVector>& flux)
{
  ASSERT(!scaled_constraint_);

  // preparing upwind data
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

  AmanziMesh::Entity_ID_List faces, cells;
  std::vector<int> dirs;

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1); Kc(0,0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    // mean value and gradient of nonlinear factor
    double kc = (*k_cell)[0][c];
    AmanziGeometry::Point kgrad(dim);
    for (int i = 0; i < dim; i++) kgrad[i] = (*k_grad)[i][c];
 
    // upwinded values of nonlinear factor
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    std::vector<double> kf(nfaces, 1.0); 
    if (k_twin == Teuchos::null) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    } else {
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        kf[n] = (c == cells[0]) ? (*k_face)[0][f] : (*k_twin)[0][f];
      }
    }

    WhetStone::DenseMatrix Wff(nfaces, nfaces);
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
* Basic routine of each operator: creation of matrices.
****************************************************************** */
void OperatorDiffusionMFD::UpdateMatricesMixed_(
    const Teuchos::Ptr<const CompositeVector>& flux)
{
  // un-rolling upwind data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_twin = Teuchos::null;
  if (k_ != Teuchos::null) {
    if (k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");
    if (k_->HasComponent("twin")) k_twin = k_->ViewComponent("twin", true);
  }
  if (upwind_ == OPERATOR_UPWIND_FLUX || 
      upwind_ == OPERATOR_UPWIND_AMANZI_ARTIFICIAL_DIFFUSION ||
      upwind_ == OPERATOR_UPWIND_AMANZI_DIVK) {
    k_face = k_->ViewComponent("face", true);
  }

  // update matrix blocks
  AmanziMesh::Entity_ID_List faces, cells;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Wff = Wff_cells_[c];
    WhetStone::DenseMatrix Acell(nfaces + 1, nfaces + 1);

    // Update terms due to nonlinear coefficient
    double kc(1.0);
    std::vector<double> kf(nfaces, 1.0); 
    if (upwind_ == OPERATOR_UPWIND_AMANZI_ARTIFICIAL_DIFFUSION) {
      kc = k_cell.get() ? (*k_cell)[0][c] : 1.0;
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    } else if (upwind_ == OPERATOR_UPWIND_AMANZI_DIVK && k_twin == Teuchos::null) {
      kc = k_cell.get() ? (*k_cell)[0][c] : 1.0;
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    } else if (upwind_ == OPERATOR_UPWIND_AMANZI_DIVK && k_twin != Teuchos::null) {
      kc = k_cell.get() ? (*k_cell)[0][c] : 1.0;
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        kf[n] = (c == cells[0]) ? (*k_face)[0][f] : (*k_twin)[0][f];
      }
    } else if (upwind_ == OPERATOR_UPWIND_NONE && k_cell != Teuchos::null) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    } else if (upwind_ == OPERATOR_UPWIND_FLUX) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    }
      
    if (upwind_ != OPERATOR_UPWIND_AMANZI_DIVK) {
      if (!scaled_constraint_) {
        // not scaled constraint: kr > 0
        double matsum = 0.0;  // elimination of mass matrix
        for (int n = 0; n < nfaces; n++) {
          double rowsum = 0.0;
          for (int m = 0; m < nfaces; m++) {
            double tmp = Wff(n, m) * kf[n];
            rowsum += tmp;
            Acell(n, m) = tmp; ASSERT(std::abs(tmp) < 1.e20);
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

      } else {
        // scaled constraint: kr >= 0
        double matsum = 0.0;  // elimination of mass matrix
        for (int n = 0; n < nfaces; n++) {
          double rowsum = 0.0;
          for (int m = 0; m < nfaces; m++) {
            double tmp = Wff(n, m);
            rowsum += tmp;
            Acell(n, m) = tmp;
          }

          Acell(n, nfaces) = -rowsum;
          matsum += rowsum * kf[n];
        }
        Acell(nfaces, nfaces) = matsum;

        for (int n = 0; n < nfaces; n++) {
          double colsum = 0.0;
          for (int m = 0; m < nfaces; m++) colsum += Acell(m, n) * kf[m];
          Acell(nfaces, n) = -colsum;
        }
      }
    }

    // Amanzi's first upwind: add additional flux 
    if (upwind_ == OPERATOR_UPWIND_AMANZI_ARTIFICIAL_DIFFUSION) {
      ASSERT(!scaled_constraint_);
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        double alpha = (*k_face)[0][f] - kc;
        if (alpha > 0) {
          alpha *= Wff(n, n);
          Acell(n, n) += alpha;
          Acell(n, nfaces) -= alpha;
          Acell(nfaces, n) -= alpha;
          Acell(nfaces, nfaces) += alpha;
        }
      }
    }

    // Amanzi's second upwind: replace the matrix
    if (upwind_ == OPERATOR_UPWIND_AMANZI_DIVK) {
      ASSERT(!scaled_constraint_);
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
* Calculate elemental inverse mass matrices.                                           
****************************************************************** */
void OperatorDiffusionMFD::UpdateMatricesNodal_()
{
  ASSERT(!scaled_constraint_);

  // update matrix blocks
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  mfd.ModifyStabilityScalingFactor(factor_);

  AmanziMesh::Entity_ID_List nodes;

  nfailed_primary_ = 0;

  WhetStone::Tensor K(2,1); K(0,0) = 1.0;
  
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
****************************************************************** */
void OperatorDiffusionMFD::UpdateMatricesTPFA_()
{
  // This does not seem to consider Krel? --etc

  // populate transmissibilities
  WhetStone::MFD3D_Diffusion mfd(mesh_);

  CompositeVectorSpace cv_space;
  cv_space.SetMesh(mesh_);
  cv_space.SetGhosted(true);
  cv_space.SetComponent("face", AmanziMesh::FACE, 1);

  Teuchos::RCP<CompositeVector> T = Teuchos::RCP<CompositeVector>(new CompositeVector(cv_space, true));
  Epetra_MultiVector& Ttmp = *T->ViewComponent("face", true);

  WhetStone::Tensor Kc(mesh_->space_dimension(),1); Kc(0,0) = 1.0;
  AmanziMesh::Entity_ID_List cells, faces;
  Ttmp.PutScalar(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    if (K_.get()) Kc = (*K_)[c];
    if (Kc.isZero()) continue;  // We skip zero matrices

    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix Mff(nfaces, nfaces);
    mfd.MassMatrixInverseTPFA(c, Kc, Mff);
   
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      Ttmp[0][f] += 1.0 / Mff(n, n);
    }
  }
  T->GatherGhostedToMaster();
 
  // populate the global matrix
  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
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
* Apply boundary conditions to the local matrices
****************************************************************** */
void OperatorDiffusionMFD::ApplyBCs(bool primary)
{
  if (!exclude_primary_terms_) {
    if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL
                             | OPERATOR_SCHEMA_DOFS_FACE
                             | OPERATOR_SCHEMA_DOFS_CELL)) {
      ASSERT(bcs_.size() == 1);
      ApplyBCs_Mixed_(*bcs_[0], primary);
    
    } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_FACE
            | OPERATOR_SCHEMA_DOFS_CELL)) {
      ASSERT(bcs_.size() == 1);
      ApplyBCs_Cell_(*bcs_[0], primary);
    
    } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL
            | OPERATOR_SCHEMA_DOFS_NODE)) {
      Teuchos::RCP<BCs> bc_f, bc_n;
      for (std::vector<Teuchos::RCP<BCs> >::iterator bc = bcs_.begin();
           bc != bcs_.end(); ++bc) {
        if ((*bc)->type() == OPERATOR_BC_TYPE_FACE) {
          bc_f = *bc;
        } else if ((*bc)->type() == OPERATOR_BC_TYPE_NODE) {
          bc_n = *bc;
        }
      }
      ApplyBCs_Nodal_(bc_f.ptr(), bc_n.ptr(), primary);
    }
  }

  if (jac_op_ != Teuchos::null) {
    AmanziMesh::Entity_ID_List cells, nodes;

    const std::vector<int>& bc_model = bcs_[0]->bc_model();
    const std::vector<double>& bc_value = bcs_[0]->bc_value();
    const std::vector<double>& bc_mixed = bcs_[0]->bc_mixed();
    ASSERT(bc_model.size() == nfaces_wghost);
    ASSERT(bc_value.size() == nfaces_wghost);

    for (int f = 0; f != nfaces_owned; ++f) {
      WhetStone::DenseMatrix& Aface = jac_op_->matrices[f];

      if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        jac_op_->matrices_shadow[f] = Aface;
        Aface *= 0.0;
      }
    }
  }
}


/* ******************************************************************
* Apply BCs on face values
****************************************************************** */
void OperatorDiffusionMFD::ApplyBCs_Mixed_(BCs& bc, bool primary)
{
  // apply diffusion type BCs to FACE-CELL system
  AmanziMesh::Entity_ID_List faces;

  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<double>& bc_value = bc.bc_value();
  const std::vector<double>& bc_mixed = bc.bc_mixed();
  ASSERT(bc_model.size() == nfaces_wghost);
  ASSERT(bc_value.size() == nfaces_wghost);

  global_op_->rhs()->PutScalarGhosted(0.);
  Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");

  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
        
    bool flag(true);
    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      double value = bc_value[f];

      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        if (flag) {  // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        for (int m = 0; m < nfaces; m++) {
          if (bc_model[faces[m]] != OPERATOR_BC_DIRICHLET)
            rhs_face[0][faces[m]] -= Acell(m, n) * value;
          Acell(n, m) = Acell(m, n) = 0.0;
        }

        if (primary) {
          rhs_face[0][f] = value;
          Acell(n,n) = 1.0;
        }

        rhs_cell[0][c] -= Acell(nfaces, n) * value;
        Acell(nfaces, n) = 0.0;
        Acell(n, nfaces) = 0.0;
      } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        rhs_face[0][f] -= value * mesh_->face_area(f);
      } else if (bc_model[f] == OPERATOR_BC_MIXED) {
        if (flag) {  // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        double area = mesh_->face_area(f);
        rhs_face[0][f] -= value * area;
        Acell(n, n) += bc_mixed[f] * area;
      }
    }
  }

  global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* Apply BCs on cell operators
****************************************************************** */
void OperatorDiffusionMFD::ApplyBCs_Cell_(BCs& bc, bool primary)
{
  // apply diffusion type BCs to CELL system
  AmanziMesh::Entity_ID_List cells;

  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<double>& bc_value = bc.bc_value();
  const std::vector<double>& bc_mixed = bc.bc_mixed();
  ASSERT(bc_model.size() == nfaces_wghost);
  ASSERT(bc_value.size() == nfaces_wghost);

  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");
    
  for (int f = 0; f != nfaces_owned; ++f) {
    WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
      
    if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      rhs_cell[0][cells[0]] += bc_value[f] * Aface(0, 0);
    }
    else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
      local_op_->matrices_shadow[f] = Aface;
      
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      rhs_cell[0][cells[0]] -= bc_value[f] * mesh_->face_area(f);
      Aface *= 0.0;
    }
    // solve system of two equations in three unknowns
    else if (bc_model[f] == OPERATOR_BC_MIXED) {
      local_op_->matrices_shadow[f] = Aface;
      
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      double area = mesh_->face_area(f);
      double factor = area / (1.0 + bc_mixed[f] * area / Aface(0, 0));
      rhs_cell[0][cells[0]] -= bc_value[f] * factor;
      Aface(0, 0) = bc_mixed[f] * factor;
    }
  }
}


/* ******************************************************************
* Apply BCs on cell operators
****************************************************************** */
void OperatorDiffusionMFD::ApplyBCs_Nodal_(const Teuchos::Ptr<BCs>& bc_f,
                                           const Teuchos::Ptr<BCs>& bc_v,
                                           bool primary)
{
  AmanziMesh::Entity_ID_List faces, nodes, cells;

  global_op_->rhs()->PutScalarGhosted(0.);
  Epetra_MultiVector& rhs_node = *global_op_->rhs()->ViewComponent("node", true);

  int nn(0), nm(0);
  for (int c=0; c!=ncells_owned; ++c) {
    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

    if (bc_f != Teuchos::null) {
      const std::vector<int>& bc_model = bc_f->bc_model();
      const std::vector<double>& bc_value = bc_f->bc_value();
      const std::vector<double>& bc_mixed = bc_f->bc_mixed();

      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      for (int n=0; n!=nfaces; ++n) {
        int f = faces[n];

        if (bc_model[f] == OPERATOR_BC_NEUMANN) {
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
        } else if (bc_model[f] == OPERATOR_BC_MIXED) {
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

      for (int n=0; n!=nnodes; ++n) {
        int v = nodes[n];
        double value = bc_value[v];

        if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
          if (flag) {  // make a copy of cell-based matrix
            local_op_->matrices_shadow[c] = Acell;
            flag = false;
          }
          for (int m = 0; m < nnodes; m++) {
            if (bc_model[nodes[m]] != OPERATOR_BC_DIRICHLET)
              rhs_node[0][nodes[m]] -= Acell(m, n) * value;
            Acell(n, m) = Acell(m, n) = 0.0;
          }

          if (primary) {
            mesh_->node_get_cells(v, AmanziMesh::USED, &cells);
            rhs_node[0][v] += value / cells.size();
            Acell(n,n) = 1.0 / cells.size();
          }
        }
      }
    }
  } 

  global_op_->rhs()->GatherGhostedToMaster("node", Add);
}


/* ******************************************************************
* Modify operator by addition approximation of Newton corection.
* We ignore the right-hand side for the moment.
****************************************************************** */
void OperatorDiffusionMFD::AddNewtonCorrectionCell_(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  // hack: ignore correction if no flux provided.
  if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

  // only works on upwinded methods
  if (upwind_ == OPERATOR_UPWIND_NONE) return;

  const Epetra_MultiVector& kf = *k_->ViewComponent("face");
  const Epetra_MultiVector& dkdp_f = *dkdp_->ViewComponent("face");
  const Epetra_MultiVector& flux_f = *flux->ViewComponent("face");

  // populate the local matrices
  AmanziMesh::Entity_ID_List cells;
  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface.PutScalar(0.0);

    double v = flux_f[0][f];
    double vmod = kf[0][f] > 0. ? fabs(v) * dkdp_f[0][f] / kf[0][f] : 0.;
    if (scalar_rho_mu_) {
      vmod *= rho_;
    } else {
      ASSERT(false);
    }

    // interior face
    int i, dir, c1, c2;
    c1 = cells[0];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c1, &dir);
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
* This method is entirely unclear why it exists and should be documented.
****************************************************************** */
void OperatorDiffusionMFD::ModifyMatrices(const CompositeVector& u)
{
  if (local_op_schema_ != (OPERATOR_SCHEMA_BASE_CELL |
                           OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
    std::cout << "Schema " << global_op_schema_ << " is not supported" << std::endl;
    ASSERT(0);
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
* WARNING: Since diffusive flux is not continuous, we derive it only
* once (using flag) and in exactly the same manner as other routines.
* **************************************************************** */
void OperatorDiffusionMFD::UpdateFlux(const CompositeVector& u, CompositeVector& flux)
{
  // Initialize intensity in ghost faces.
  flux.PutScalar(0.0);
  u.ScatterMasterToGhosted("face");

  const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
  const Epetra_MultiVector& u_face = *u.ViewComponent("face", true);
  Epetra_MultiVector& flux_data = *flux.ViewComponent("face", true);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<int> flag(nfaces_wghost, 0);

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
      if (f < nfaces_owned && !flag[f]) {
        flux_data[0][f] -= av(n) * dirs[n];
        flag[f] = 1;
      }
    }
  }
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.
****************************************************************** */
void OperatorDiffusionMFD::CreateMassMatrices_()
{
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  mfd.ModifyStabilityScalingFactor(factor_);

  bool surface_mesh = (mesh_->cell_dimension() != mesh_->space_dimension());
  AmanziMesh::Entity_ID_List faces;

  Wff_cells_.resize(ncells_owned);

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0,0) = 1.0;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    int ok;
    if (K_.get()) Kc = (*K_)[c];
    WhetStone::DenseMatrix Wff(nfaces, nfaces);

    if (surface_mesh) {
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
          ok = mfd.MassMatrixInverseOptimizedScaled(c, Kc, Wff);
        } else if(method == WhetStone::DIFFUSION_TPFA) {
          ok = mfd.MassMatrixInverseTPFA(c, Kc, Wff);
        } else if(method == WhetStone::DIFFUSION_SUPPORT_OPERATOR) {
          ok = mfd.MassMatrixInverseSO(c, Kc, Wff);
        } else if(method == WhetStone::DIFFUSION_POLYHEDRA_SCALED) {
          ok = mfd.MassMatrixInverseScaled(c, Kc, Wff);
        }
      }
    }

    if (scalar_rho_mu_) {
      Wff *= rho_ / mu_;
    } else {
      const Epetra_MultiVector& rho = *rho_cv_->ViewComponent("cell");
      const Epetra_MultiVector& mu = *mu_cv_->ViewComponent("cell");
      Wff *= rho[0][c] / mu[0][c];
    }

    Wff_cells_[c] = Wff;

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("OperatorDiffusion: unexpected failure in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void OperatorDiffusionMFD::InitDiffusion_(Teuchos::ParameterList& plist)
{
  // Determine discretization
  std::string primary = plist.get<std::string>("discretization primary");
  std::string secondary = plist.get<std::string>("discretization secondary", primary);

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
    msg << "OperatorDiffusion: primary discretization method \"" << primary << "\" is not supported.";
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
    msg << "OperatorDiffusion: secondary discretization method \"" << secondary << "\" is not supported.";
    Exceptions::amanzi_throw(msg);
  }

  // Define stencil for the MFD diffusion method.
  std::vector<std::string> names;
  names = plist.get<Teuchos::Array<std::string> > ("schema").toVector();

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
    msg << "OperatorDiffusion: \"schema\" must be CELL, FACE+CELL, or NODE";
    Exceptions::amanzi_throw(msg);
  }

  // define stencil for the assembled matrix
  int schema_prec_dofs = 0;
  if (plist.isParameter("preconditioner schema")) {
    names = plist.get<Teuchos::Array<std::string> > ("preconditioner schema").toVector();
    for (int i = 0; i < names.size(); i++) {
      if (names[i] == "cell") {
        schema_prec_dofs += OPERATOR_SCHEMA_DOFS_CELL;
      } else if (names[i] == "node") {
        schema_prec_dofs += OPERATOR_SCHEMA_DOFS_NODE;
      } else if (names[i] == "face") {
        schema_prec_dofs += OPERATOR_SCHEMA_DOFS_FACE;
      }
    } 
  } else {
    schema_prec_dofs = schema_dofs;
  }

  // create or check the existing Operator
  int global_op_schema = schema_prec_dofs;  
  if (global_op_ == Teuchos::null) {
    global_op_schema_ = global_op_schema;

    // build the CVS from the global schema
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);

    if (global_op_schema & OPERATOR_SCHEMA_DOFS_CELL)
      cvs->AddComponent("cell", AmanziMesh::CELL, 1);
    if (global_op_schema & OPERATOR_SCHEMA_DOFS_FACE)
      cvs->AddComponent("face", AmanziMesh::FACE, 1);
    if (global_op_schema & OPERATOR_SCHEMA_DOFS_NODE)
      cvs->AddComponent("node", AmanziMesh::NODE, 1);

    // choose the Operator from the prec schema
    Teuchos::ParameterList operator_list = plist.sublist("operator");
    if (schema_prec_dofs == OPERATOR_SCHEMA_DOFS_NODE) {
      global_op_ = Teuchos::rcp(new Operator_Node(cvs, plist));
    } else if (schema_prec_dofs == OPERATOR_SCHEMA_DOFS_CELL) {
      //      cvs->AddComponent("face", AmanziMesh::FACE, 1);
      //      global_op_ = Teuchos::rcp(new Operator_FaceCellScc(cvs, plist));
      global_op_ = Teuchos::rcp(new Operator_Cell(cvs, plist, schema_prec_dofs));
    } else if (schema_prec_dofs == OPERATOR_SCHEMA_DOFS_FACE) {
      cvs->AddComponent("cell", AmanziMesh::CELL, 1);
      global_op_ = Teuchos::rcp(new Operator_FaceCellSff(cvs, plist));
    } else if (schema_prec_dofs == (OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
      global_op_ = Teuchos::rcp(new Operator_FaceCell(cvs, plist));
    } else {
      Errors::Message msg;
      msg << "OperatorDiffusion: \"preconditioner schema\" must be NODE, CELL, FACE, or FACE+CELL";
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
    } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL |
            OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL)) {
      std::string name = "Diffusion: CELL_FACE+CELL";
      local_op_ = Teuchos::rcp(new Op_Cell_FaceCell(name, mesh_));
    } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_FACE |
            OPERATOR_SCHEMA_DOFS_CELL)) {
      if (plist.get<bool>("surface operator", false)) {
        std::string name = "Diffusion: FACE_CELL Surface";
        local_op_ = Teuchos::rcp(new Op_SurfaceFace_SurfaceCell(name, mesh_));
      } else {
        std::string name = "Diffusion: FACE_CELL";
        local_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
      }
    } else {
      ASSERT(0);
    }
    global_op_->OpPushBack(local_op_);
  }
  
  // scaled constraint -- enables zero rel perm
  scaled_constraint_ = plist.get<bool>("scaled constraint equation", false);

  // upwind options
  std::string name = plist.get<std::string>("upwind method", "none");
  if (name == "standard") {
    upwind_ = OPERATOR_UPWIND_FLUX;
  } else if (name == "artificial diffusion") {  
    upwind_ = OPERATOR_UPWIND_AMANZI_ARTIFICIAL_DIFFUSION;
  } else if (name == "divk") {  
    upwind_ = OPERATOR_UPWIND_AMANZI_DIVK;
  } else if (name == "second-order") {  
    upwind_ = OPERATOR_UPWIND_AMANZI_SECOND_ORDER;
  } else if (name == "none") {
    upwind_ = OPERATOR_UPWIND_NONE;  // cell-centered scheme.
  } else {
    ASSERT(false);
  }

  // Do we need to calculate Newton correction terms?
  newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;
  std::string jacobian = plist.get<std::string>("newton correction", "none");
  if (jacobian == "true jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_TRUE;
  } else if (jacobian == "approximate jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE;

    // create a local op
    jac_op_schema_ = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL;
    std::string name("Jacobian FACE_CELL");
    jac_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
    global_op_->OpPushBack(jac_op_);
  }

  // mesh info
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
}


/* ******************************************************************
* Given a set of cell values, update faces using the consistency equations,
* x_f = Aff^-1 * (y_f - Afc * x_c)
****************************************************************** */
void
OperatorDiffusionMFD::UpdateConsistentFaces(CompositeVector& u)
{
  if (consistent_face_op_ == Teuchos::null) {
    // create the op
    Teuchos::RCP<CompositeVectorSpace> cface_cvs = Teuchos::rcp(new CompositeVectorSpace());
    cface_cvs->SetMesh(mesh_)->SetGhosted()
        ->AddComponent("face", AmanziMesh::FACE, 1);

    consistent_face_op_ = Teuchos::rcp(new Operator_ConsistentFace(cface_cvs, plist_.sublist("consistent faces")));
    consistent_face_op_->OpPushBack(local_op_);
    consistent_face_op_->SymbolicAssembleMatrix();
  }

  // calculate the rhs, given by y_f - Afc * x_c
  CompositeVector& y = *consistent_face_op_->rhs();
  Epetra_MultiVector& y_f = *y.ViewComponent("face",true);
  y_f = *global_op_->rhs()->ViewComponent("face",true);
  consistent_face_op_->rhs()->PutScalarGhosted(0.);

  // y_f - Afc * x_c
  const Epetra_MultiVector& x_c = *u.ViewComponent("cell",false);
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
  consistent_face_op_->InitPreconditioner(plist_.sublist("consistent faces").sublist("preconditioner"));

  if (plist_.sublist("consistent faces").isSublist("linear solver")) {
    AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> fac;
    Teuchos::RCP<Operator> lin_solver = fac.Create(
        plist_.sublist("consistent faces").sublist("linear solver"), consistent_face_op_);

    CompositeVector u_f_copy(y);
    int ierr = lin_solver->ApplyInverse(y, u_f_copy);
    *u.ViewComponent("face",false) = *u_f_copy.ViewComponent("face",false);
    ASSERT(!ierr);
  } else {
    CompositeVector u_f_copy(y);
    int ierr = consistent_face_op_->ApplyInverse(y, u);
    ASSERT(!ierr);
    *u.ViewComponent("face",false) = *u_f_copy.ViewComponent("face",false);
  }
}

  
}  // namespace Operators
}  // namespace Amanzi
