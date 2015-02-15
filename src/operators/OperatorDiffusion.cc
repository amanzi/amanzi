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

#include "Op.hh"
#include "Op_Cell_Node.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Face_Cell.hh"

#include "OperatorDefs.hh"
//#include "Operator_Node.hh"
#include "Operator_FaceCell.hh"
#include "Operator_FaceCellScc.hh"
#include "Operator_FaceCellSff.hh"

#include "OperatorDiffusion.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator.                                           
****************************************************************** */
void OperatorDiffusion::Setup(
    std::vector<WhetStone::Tensor>& K, 
    Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
    double rho, double mu)
{
  K_ = &K;
  k_ = k;
  dkdp_ = dkdp;

  rho_ = rho;
  mu_ = mu;
  scalar_rho_mu_ = true;

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

  if (local_op_schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
    ASSERT(K_->size() == ncells_owned);
    CreateMassMatrices_();
  }
}


/* ******************************************************************
* Initialization of the operator.                                           
****************************************************************** */
void OperatorDiffusion::Setup(
    std::vector<WhetStone::Tensor>& K,
    Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
    Teuchos::RCP<const CompositeVector> rho, Teuchos::RCP<const CompositeVector> mu)
{
  K_ = &K;
  k_ = k;
  dkdp_ = dkdp;

  rho_cv_ = rho;
  mu_cv_ = mu;
  scalar_rho_mu_ = false;

  if (local_op_schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
    ASSERT(K_->size() == ncells_owned);
    CreateMassMatrices_();
  }
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void OperatorDiffusion::UpdateMatrices(Teuchos::RCP<const CompositeVector> flux,
                                       Teuchos::RCP<const CompositeVector> u)
{
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


/* ******************************************************************
* Second-order upwind. Mass matrices are recalculated.
****************************************************************** */
void OperatorDiffusion::UpdateMatricesMixedWithGrad_(Teuchos::RCP<const CompositeVector> flux)
{
  // preparing upwind data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_grad = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_twin = Teuchos::null;
  if (k_ != Teuchos::null) {
    k_cell = k_->ViewComponent("cell");
    k_face = k_->ViewComponent("face");
    k_grad = k_->ViewComponent("grad");
    if (k_->HasComponent("twin")) k_twin = k_->ViewComponent("twin", true);
  }

  // update matrix blocks
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);

  AmanziMesh::Entity_ID_List faces, cells;
  std::vector<int> dirs;

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
    WhetStone::Tensor& Kc = (*K_)[c];
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
void OperatorDiffusion::UpdateMatricesMixed_(Teuchos::RCP<const CompositeVector> flux)
{
  // un-rolling upwind data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_twin = Teuchos::null;
  if (k_ != Teuchos::null) {
    k_cell = k_->ViewComponent("cell");
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
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    } else if (upwind_ == OPERATOR_UPWIND_AMANZI_DIVK && k_twin == Teuchos::null) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    } else if (upwind_ == OPERATOR_UPWIND_AMANZI_DIVK && k_twin != Teuchos::null) {
      kc = (*k_cell)[0][c];
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
      double matsum = 0.0;  // elimination of mass matrix
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
    }

    // Amanzi's first upwind: add additional flux 
    if (upwind_ == OPERATOR_UPWIND_AMANZI_ARTIFICIAL_DIFFUSION) {
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
void OperatorDiffusion::UpdateMatricesNodal_()
{
  // update matrix blocks
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  mfd.ModifyStabilityScalingFactor(factor_);

  AmanziMesh::Entity_ID_List nodes;

  nfailed_primary_ = 0;
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    WhetStone::DenseMatrix Acell(nnodes, nnodes);

    int method = mfd_primary_;
    int ok = WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED;

    if (method == WhetStone::DIFFUSION_OPTIMIZED_FOR_MONOTONICITY) {
      ok = mfd.StiffnessMatrixMMatrix(c, (*K_)[c], Acell);
      method = mfd_secondary_;
    } else {
      ok = mfd.StiffnessMatrix(c, (*K_)[c], Acell);
      method = mfd_secondary_;
    }

    if (ok != WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) {
      nfailed_primary_++;
      ok = mfd.StiffnessMatrix(c, (*K_)[c], Acell);
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
void OperatorDiffusion::UpdateMatricesTPFA_()
{
  // populate transmissibilities
  WhetStone::MFD3D_Diffusion mfd(mesh_);

  CompositeVectorSpace cv_space;
  cv_space.SetMesh(mesh_);
  cv_space.SetGhosted(true);
  cv_space.SetComponent("face", AmanziMesh::FACE, 1);

  Teuchos::RCP<CompositeVector> T = Teuchos::RCP<CompositeVector>(new CompositeVector(cv_space, true));
  Epetra_MultiVector& Ttmp = *T->ViewComponent("face", true);

  AmanziMesh::Entity_ID_List cells, faces;
  Ttmp.PutScalar(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    if ((*K_)[c].isZero()) continue;  // We skip zero matrices

    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix Mff(nfaces, nfaces);
    mfd.MassMatrixInverseTPFA(c, (*K_)[c], Mff);
   
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



// /* ******************************************************************
// * Special assemble of elemental face-based matrices. 
// ****************************************************************** */
// void OperatorDiffusion::ModifyMatrices(const CompositeVector& u)
// {
//   if (schema_dofs_ != OPERATOR_SCHEMA_DOFS_CELL + OPERATOR_SCHEMA_DOFS_FACE) {
//     std::cout << "Schema " << schema_dofs_ << " is not supported" << std::endl;
//     ASSERT(0);
//   }

//   // find location of face-based matrices
//   int m = FindMatrixBlock(schema_, OPERATOR_SCHEMA_RULE_EXACT, true);
//   std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];

//   // populate the matrix
//   AmanziMesh::Entity_ID_List faces;
//   const Epetra_MultiVector& u_c = *u.ViewComponent("cell");
//   Epetra_MultiVector& rhs_f = *rhs_->ViewComponent("face", true);

//   for (int f = nfaces_owned; f < nfaces_wghost; f++) rhs_f[0][f] = 0.0;

//   for (int c = 0; c < ncells_owned; c++) {
//     mesh_->cell_get_faces(c, &faces);
//     int nfaces = faces.size();

//     WhetStone::DenseMatrix& Acell = matrix[c];

//     for (int n = 0; n < nfaces; n++) {
//       int f = faces[n];
//       rhs_f[0][f] -= Acell(n, nfaces) * u_c[0][c];
//       Acell(n, nfaces) = 0.0;
//       Acell(nfaces, n) = 0.0;
//     }
//   }

//   // Assemble all right-hand sides
//   rhs_->GatherGhostedToMaster("face", Add);
// }


/* ******************************************************************
* WARNING: Since diffusive flux is not continuous, we derive it only
* once (using flag) and in exactly the same manner as other routines.
* **************************************************************** */
void OperatorDiffusion::UpdateFlux(const CompositeVector& u, CompositeVector& flux)
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
* Modify operator by addition approximation of NEwton corection.
****************************************************************** */
void OperatorDiffusion::AddNewtonCorrection(Teuchos::RCP<const CompositeVector> flux)
{
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.
****************************************************************** */
void OperatorDiffusion::CreateMassMatrices_()
{
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  mfd.ModifyStabilityScalingFactor(factor_);

  bool surface_mesh = (mesh_->cell_dimension() != mesh_->space_dimension());
  AmanziMesh::Entity_ID_List faces;

  Wff_cells_.resize(ncells_owned);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    int ok;
    WhetStone::Tensor& Kc = (*K_)[c];
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
void OperatorDiffusion::InitDiffusion_(Teuchos::ParameterList& plist)
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
      ASSERT(0);
      //      global_op_ = Teuchos::rcp(new Operator_Node(cvs, plist));
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

  // create the local Op and register it with the global Operator
  if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_NODE)) {
    std::string name = "Diffusion: CELL_NODE";
    local_op_ = Teuchos::rcp(new Op_Cell_Node(name, mesh_));
  } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL |
          OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL)) {
    std::string name = "Diffusion: CELL_FACE+CELL";
    local_op_ = Teuchos::rcp(new Op_Cell_FaceCell(name, mesh_));
  } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_FACE |
          OPERATOR_SCHEMA_DOFS_CELL)) {
    std::string name = "Diffusion: FACE_CELL";
    local_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
  } else {
    ASSERT(0);
  }
  global_op_->OpPushBack(local_op_);
  
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

  // mesh info
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
}

}  // namespace Operators
}  // namespace Amanzi
