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
#include "Epetra_FECrsGraph.h"

#include "errors.hh"
#include "WhetStoneDefs.hh"
#include "mfd3d_diffusion.hh"

#include "PreconditionerFactory.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusion.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator.                                           
****************************************************************** */
void OperatorDiffusion::InitOperator(
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
      upwind_ == OPERATOR_UPWIND_AMANZI_MFD) { 
    ASSERT(k->HasComponent("face"));
  }

  if (schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
    CreateMassMatrices_();
  }
}


/* ******************************************************************
* Initialization of the operator.                                           
****************************************************************** */
void OperatorDiffusion::InitOperator(
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

  if (schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
    CreateMassMatrices_();
  }
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void OperatorDiffusion::UpdateMatrices(Teuchos::RCP<const CompositeVector> flux,
                                       Teuchos::RCP<const CompositeVector> u)
{
  if (schema_dofs_ == OPERATOR_SCHEMA_DOFS_NODE) {
    UpdateMatricesNodal_();
  } else if (schema_dofs_ == OPERATOR_SCHEMA_DOFS_CELL + OPERATOR_SCHEMA_DOFS_FACE) {
    UpdateMatricesMixed_(flux);
  } else if (schema_dofs_ == OPERATOR_SCHEMA_DOFS_CELL) {
    UpdateMatricesTPFA_();
  }
}


/* ******************************************************************
* Basic routine of each operator: creation of matrices.
****************************************************************** */
void OperatorDiffusion::UpdateMatricesMixed_(Teuchos::RCP<const CompositeVector> flux)
{
  // find location of matrix blocks
  int schema_dofs = OPERATOR_SCHEMA_DOFS_CELL + OPERATOR_SCHEMA_DOFS_FACE;
  int m = FindMatrixBlock(schema_dofs, OPERATOR_SCHEMA_RULE_EXACT, false);
  bool flag = (m >= 0); 

  if (flag == false) { 
    m = blocks_.size();
    blocks_schema_.push_back(OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL);
    blocks_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
    blocks_shadow_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
  }
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[m];
  WhetStone::DenseMatrix null_matrix;

  // preparing upwind data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) k_cell = k_->ViewComponent("cell");
  if (upwind_ == OPERATOR_UPWIND_FLUX || 
      upwind_ == OPERATOR_UPWIND_AMANZI_ARTIFICIAL_DIFFUSION ||
      upwind_ == OPERATOR_UPWIND_AMANZI_MFD) {
    k_face = k_->ViewComponent("face", true);
  }

  // update matrix blocks
  AmanziMesh::Entity_ID_List faces;
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
    } else if (upwind_ == OPERATOR_UPWIND_AMANZI_MFD) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    } else if (upwind_ == OPERATOR_UPWIND_NONE && k_cell != Teuchos::null) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    } else if (upwind_ == OPERATOR_UPWIND_FLUX) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    }

    if (upwind_ != OPERATOR_UPWIND_AMANZI_MFD) {
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
    if (upwind_ == OPERATOR_UPWIND_AMANZI_MFD) {
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

    if (flag) {
      matrix[c] += Acell;
    } else {
      matrix.push_back(Acell);
      matrix_shadow.push_back(null_matrix);
    }
  }
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.                                           
****************************************************************** */
void OperatorDiffusion::UpdateMatricesNodal_()
{
  // find location of matrix blocks
  int schema_dofs = OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_NODE;
  int m = FindMatrixBlock(schema_dofs, OPERATOR_SCHEMA_RULE_EXACT, false);
  bool flag = (m >= 0);

  if (flag == false) { 
    m = blocks_.size();
    blocks_schema_.push_back(OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_NODE);
    blocks_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
    blocks_shadow_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
  }
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[m];
  WhetStone::DenseMatrix null_matrix;

  // update matrix blocks
  int dim = mesh_->space_dimension();
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

    if (flag) {
      matrix[c] += Acell;
    } else {
      matrix.push_back(Acell);
      matrix_shadow.push_back(null_matrix);
    }
  }
}


/* ******************************************************************
* Calculate and assemble fluxes using the TPFA scheme.
****************************************************************** */
void OperatorDiffusion::UpdateMatricesTPFA_()
{
  // find location of matrix blocks
  int schema_dofs = OPERATOR_SCHEMA_BASE_FACE + OPERATOR_SCHEMA_DOFS_CELL;
  int m = FindMatrixBlock(schema_dofs, OPERATOR_SCHEMA_RULE_EXACT, false);
  bool flag = (m >= 0);

  if (flag == false) { 
    m = blocks_.size();
    blocks_schema_.push_back(OPERATOR_SCHEMA_BASE_FACE + OPERATOR_SCHEMA_DOFS_CELL);
    blocks_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
    blocks_shadow_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
  }
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[m];
  WhetStone::DenseMatrix null_matrix;

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

    if (flag) {
      matrix[f] += Aface;
    } else {
      matrix.push_back(Aface);
      matrix_shadow.push_back(null_matrix);
    }
  }
}


/* ******************************************************************
* A small factory for assembling of matrices for preconditioners.
****************************************************************** */
void OperatorDiffusion::AssembleMatrix(int schema)
{
  if (special_assembling_ != 0) {
    // We do not need it since preconditoner creates a special matrix.
    // AssembleMatrixSpecialSff_();
  } else {
    Operator::AssembleMatrix(schema);
  }
}


/* ******************************************************************
* Special assemble of elemental face-based matrices. 
****************************************************************** */
void OperatorDiffusion::ModifyMatrices(const CompositeVector& u)
{
  if (schema_dofs_ != OPERATOR_SCHEMA_DOFS_CELL + OPERATOR_SCHEMA_DOFS_FACE) {
    std::cout << "Schema " << schema_dofs_ << " is not supported" << std::endl;
    ASSERT(0);
  }

  // find location of face-based matrices
  int m = FindMatrixBlock(schema_, OPERATOR_SCHEMA_RULE_EXACT, true);
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];

  // populate the matrix
  AmanziMesh::Entity_ID_List faces;
  const Epetra_MultiVector& u_c = *u.ViewComponent("cell");
  Epetra_MultiVector& rhs_f = *rhs_->ViewComponent("face", true);

  for (int f = nfaces_owned; f < nfaces_wghost; f++) rhs_f[0][f] = 0.0;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Acell = matrix[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      rhs_f[0][f] -= Acell(n, nfaces) * u_c[0][c];
      Acell(n, nfaces) = 0.0;
      Acell(nfaces, n) = 0.0;
    }
  }

  // Assemble all right-hand sides
  rhs_->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* The cell-based and face-based d.o.f. are packed together into 
* the X and Y vectors.
****************************************************************** */
int OperatorDiffusion::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  int ierr;
  if (special_assembling_ == 1) {
    ierr = ApplyInverseSpecialSff_(X, Y);
  } else if (special_assembling_ == 2) {
    ierr = ApplyInverseSpecialScc_(X, Y);
  } else {
    ierr = Operator::ApplyInverse(X, Y);
  }
  return ierr;
}

 
/* ******************************************************************
* The cell-based and face-based d.o.f. are packed together into 
* the X and Y vectors.
****************************************************************** */
int OperatorDiffusion::ApplyInverseSpecialSff_(const CompositeVector& X, CompositeVector& Y) const
{
  // Y = X;
  // return 0;

  // find the block of matrices
  int schema_dofs = OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL;
  int m = FindMatrixBlock(schema_dofs, OPERATOR_SCHEMA_RULE_SUBSET, true);
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];

  // apply preconditioner inversion
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  Epetra_MultiVector& Yc = *Y.ViewComponent("cell");
  Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);

  // Temporary cell and face vectors.
  CompositeVector T(X);
  Epetra_MultiVector& Tf = *T.ViewComponent("face", true);

  // FORWARD ELIMINATION:  Tf = Xf - Afc inv(Acc) Xc
  AmanziMesh::Entity_ID_List faces;
  Epetra_MultiVector& diag = *diagonal_->ViewComponent("cell");

  for (int f = nfaces_owned; f < nfaces_wghost; f++) {
    Tf[0][f] = 0.0;
  }

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Acell = matrix[c];

    double tmp = Xc[0][c] / (Acell(nfaces, nfaces) + diag[0][c]);
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      Tf[0][f] -= Acell(n, nfaces) * tmp;
    }
  }

  // Solve the Schur complement system Sff * Yf = Tf.
  Epetra_MultiVector& Tf_short = *T.ViewComponent("face", false);
  Epetra_MultiVector& Yf_short = *Y.ViewComponent("face", false);

  T.GatherGhostedToMaster("face", Add);

  preconditioner_->ApplyInverse(Tf_short, Yf_short);

  Y.ScatterMasterToGhosted("face");

  // BACKWARD SUBSTITUTION:  Yc = inv(Acc) (Xc - Acf Yf)
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Acell = matrix[c];

    double tmp = Xc[0][c];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      tmp -= Acell(nfaces, n) * Yf[0][f];
    }
    Yc[0][c] = tmp / (Acell(nfaces, nfaces) + diag[0][c]);
  }

  return 0;
}


/* ******************************************************************
* The cell-based and face-based d.o.f. are packed together into 
* the X and Y vectors.
****************************************************************** */
int OperatorDiffusion::ApplyInverseSpecialScc_(const CompositeVector& X, CompositeVector& Y) const
{
  // Y = X;
  // return 0;

  // find the block of matrices
  int schema_dofs = OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL;
  int m = FindMatrixBlock(schema_dofs, OPERATOR_SCHEMA_RULE_SUBSET, true);
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];

  // apply preconditioner inversion
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  Epetra_MultiVector& Yc = *Y.ViewComponent("cell");
  Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);

  // Temporary cell and face vectors.
  CompositeVector T(X);
  Epetra_MultiVector& Tf = *T.ViewComponent("face", true);
  Epetra_MultiVector& Tc = *T.ViewComponent("cell");

  // populate Aff
  AmanziMesh::Entity_ID_List faces;

  Epetra_MultiVector& diag_face = *diagonal_->ViewComponent("face", true);
  Tf = diag_face;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Acell = matrix[c];
   
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      Tf[0][f] += Acell(n, n);
    }
  }
  T.GatherGhostedToMaster("face");

  // FORWARD ELIMINATION:  Tc = Xc - Acf inv(Aff) Xf
  T.ScatterMasterToGhosted("face");
  X.ScatterMasterToGhosted("face");

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Acell = matrix[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      Tc[0][c] -= Acell(n, nfaces) / Tf[0][f] * Xf[0][f];
    }
  }

  // Solve the Schur complement system Scc * Yc = Tc.
  preconditioner_->ApplyInverse(Tc, Yc);

  // BACKWARD SUBSTITUTION:  Yf = inv(Aff) (Xf - Afc Yc)
  Yf = Xf;
  for (int f = nfaces_owned; f < nfaces_wghost; f++) {
    Yf[0][f] = 0.0;
  }

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Acell = matrix[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      Yf[0][f] -= Acell(nfaces, n) * Yc[0][c];
    }
  }

  Y.GatherGhostedToMaster("face", Add);
  for (int f = 0; f < nfaces_owned; f++) Yf[0][f] /= Tf[0][f];

  return 0;
}


/* ******************************************************************
* Initialization of the preconditioner                                                 
****************************************************************** */
void OperatorDiffusion::InitPreconditioner(
    const std::string& prec_name, const Teuchos::ParameterList& plist)
{
  if (special_assembling_ == 1) { 
    InitPreconditionerSpecialSff_(prec_name, plist);
  } else if (special_assembling_ == 2) { 
    InitPreconditionerSpecialScc_(prec_name, plist);
  } else {
    Operator::InitPreconditioner(prec_name, plist);
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Initializing preconditioner, ||A||=" << A_->NormOne() << std::endl; 
  }
}


/* ******************************************************************
* Routine assembles the Schur complement for face-based degrees 
* of freedom, Sff = Aff - Afc Acc^{-1} Acf. 
****************************************************************** */
void OperatorDiffusion::InitPreconditionerSpecialSff_(
    const std::string& prec_name, const Teuchos::ParameterList& plist)
{
  const std::vector<int>& bc_model = GetBCofType(OPERATOR_BC_TYPE_FACE)->bc_model();
  const std::vector<double>& bc_value = GetBCofType(OPERATOR_BC_TYPE_FACE)->bc_value();

  // find the block of matrices
  int schema_dofs = OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL;
  int m = FindMatrixBlock(schema_dofs, OPERATOR_SCHEMA_RULE_SUBSET, true);

  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];

  // create a face-based stiffness matrix from A.
  A_->PutScalar(0.0);
  A_off_->PutScalar(0.0);
  int ndof = A_->NumMyRows();

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  AmanziMesh::Entity_ID_List faces;

  int lid[OPERATOR_MAX_FACES];
  int gid[OPERATOR_MAX_FACES];
  double values[OPERATOR_MAX_FACES];

  Epetra_MultiVector& diag = *diagonal_->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix Scell(nfaces, nfaces);
    WhetStone::DenseMatrix& Acell = matrix[c];

    double tmp = Acell(nfaces, nfaces) + diag[0][c];
    if (tmp == 0.0 && (*K_)[c].isZero()) continue;  // We skip zero matrices

    for (int n = 0; n < nfaces; n++) {
      for (int m = 0; m < nfaces; m++) {
        Scell(n, m) = Acell(n, m) - Acell(n, nfaces) * Acell(nfaces, m) / tmp;
      }
    }

    for (int n = 0; n < nfaces; n++) {  // Symbolic boundary conditions
      int f = faces[n];
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        for (int m = 0; m < nfaces; m++) Scell(n, m) = Scell(m, n) = 0.0;
        Scell(n, n) = 1.0;
      }
    }

    for (int n = 0; n < nfaces; n++) {
      lid[n] = faces[n];
      gid[n] = fmap_wghost.GID(faces[n]);
    }
    for (int n = 0; n < nfaces; n++) {
      for (int m = 0; m < nfaces; m++) values[m] = Scell(n, m);
      if (lid[n] < ndof) {
        A_->SumIntoMyValues(lid[n], nfaces, values, lid);
      } else {
        A_off_->SumIntoMyValues(lid[n] - ndof, nfaces, values, lid);
      }
    }
  }

  // Scatter off proc to their locations
  A_->Export(*A_off_, *exporter_, Add);
      
  const Epetra_Map& map = A_->RowMap();
  A_->FillComplete(map, map);

  // redefine (if necessary) preconditioner since only 
  // one preconditioner is allowed.
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(prec_name, plist);
  preconditioner_->Update(A_);
}


/* ******************************************************************
* Routine assembles the Schur complement for face-based degrees 
* of freedom, Scc = Acc - Acf Aff^{-1} Afc. 
****************************************************************** */
void OperatorDiffusion::InitPreconditionerSpecialScc_(
    const std::string& prec_name, const Teuchos::ParameterList& plist)
{
  const std::vector<int>& bc_model = GetBCofType(OPERATOR_BC_TYPE_FACE)->bc_model();

  // find location of matrix blocks
  int schema_dofs = OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL;
  int m = FindMatrixBlock(schema_dofs, OPERATOR_SCHEMA_RULE_SUBSET, true);

  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];

  // create a cell-based stiffness matrix from A.
  A_->PutScalar(0.0);
  A_off_->PutScalar(0.0);
  int ndof = A_->NumMyRows();

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  // populate ecofficients that form transmissibility
  CompositeVectorSpace cv_space;
  cv_space.SetMesh(mesh_);
  cv_space.SetGhosted(true);
  cv_space.SetComponent("face", AmanziMesh::FACE, 2);

  Teuchos::RCP<CompositeVector> T = Teuchos::RCP<CompositeVector>(new CompositeVector(cv_space, true));
  Epetra_MultiVector& Ttmp = *T->ViewComponent("face", true);

  AmanziMesh::Entity_ID_List cells, faces;
  Epetra_MultiVector& diag_face = *diagonal_->ViewComponent("face", true);

  Ttmp.PutScalar(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    int c0 = cmap_wghost.GID(c);
    WhetStone::DenseMatrix& Acell = matrix[c];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

      if (cells.size() == 1) { 
        Ttmp[0][f] = Acell(n, nfaces);
      } else {
        int c1 = cmap_wghost.GID(cells[0]);
        int c2 = cmap_wghost.GID(cells[1]);
        int i = (c0 == std::min(c1, c2)) ? 0 : 1; 
        Ttmp[i][f] = Acell(n, nfaces);
      }
    }
  }
  T->GatherGhostedToMaster();
 
  // populate the global matrix
  int lid[2], gid[2];
  double a1, a2, values[2];

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++) {
      lid[n] = cells[n];
      gid[n] = cmap_wghost.GID(cells[n]);
    }

    a1 = Ttmp[0][f]; 
    a2 = Ttmp[1][f];
    if (ncells == 2 && gid[0] > gid[1]) {
      a1 = Ttmp[1][f]; 
      a2 = Ttmp[0][f];
    }

    double coef = fabs(a1 + a2) + diag_face[0][f];
    if (coef == 0.0) continue;

    for (int n = 0; n < ncells; n++) {
      if (n == 0) {
        values[0] = -a1 * a1 / coef;
        values[1] = -a1 * a2 / coef;
      } else {
        values[0] = -a1 * a2 / coef;
        values[1] = -a2 * a2 / coef;
      }

      if (lid[n] < ndof) {
        A_->SumIntoMyValues(lid[n], ncells, values, lid);
      } else {
        A_off_->SumIntoMyValues(lid[n] - ndof, ncells, values, lid);
      }
    }
  }

  // Scatter off proc to their locations
  A_->Export(*A_off_, *exporter_, Add);
      
  const Epetra_Map& map = A_->RowMap();
  A_->FillComplete(map, map);

  // Add diagonal (a hack)
  Epetra_Vector tmp(A_->RowMap());
  A_->ExtractDiagonalCopy(tmp);

  if (diagonal_->HasComponent("cell")) {
    Epetra_MultiVector& diag_cell = *diagonal_->ViewComponent("cell");
    for (int c = 0; c < ncells_owned; c++) {
      WhetStone::DenseMatrix& Acell = matrix[c];
      int n = Acell.NumRows() - 1;
      tmp[c] += Acell(n, n) + diag_cell[0][c];
    }
  }

  A_->ReplaceDiagonalValues(tmp);
  // std::cout << *A_ << std::endl;

  // redefine (if necessary) preconditioner since only 
  // one preconditioner is allowed.
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(prec_name, plist);
  preconditioner_->Update(A_);
}


/* ******************************************************************
* WARNING: Since diffusive flux is not continuous, we derive it only
* once (using flag) and in exactly the same manner as other routines.
* **************************************************************** */
void OperatorDiffusion::UpdateFlux(const CompositeVector& u, CompositeVector& flux)
{
  // find location of face-based matrices
  int schema_dofs = OPERATOR_SCHEMA_DOFS_CELL + OPERATOR_SCHEMA_DOFS_FACE;
  int m = FindMatrixBlock(schema_dofs, OPERATOR_SCHEMA_RULE_SUBSET, true);

  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[m];

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

    if (matrix_shadow[c].NumRows() == 0) { 
      WhetStone::DenseMatrix& Acell = matrix[c];
      Acell.Multiply(v, av, false);
    } else {
      WhetStone::DenseMatrix& Acell = matrix_shadow[c];
      Acell.Multiply(v, av, false);
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
void OperatorDiffusion::CreateMassMatrices_()
{
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  mfd.ModifyStabilityScalingFactor(factor_);

  bool surface_mesh = (mesh_->cell_dimension() != mesh_->space_dimension());
  AmanziMesh::Entity_ID_List faces;

  Wff_cells_.clear();

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

    Wff_cells_.push_back(Wff);

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("OperatorDiffusion: unexpected failure in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void OperatorDiffusion::InitDiffusion_(Teuchos::RCP<BCs> bc, Teuchos::ParameterList& plist)
{
  bc_.push_back(bc);

  // Define stencil for the MFD diffusion method.
  std::vector<std::string> names;
  names = plist.get<Teuchos::Array<std::string> > ("schema").toVector();

  schema_dofs_ = 0;
  for (int i = 0; i < names.size(); i++) {
    if (names[i] == "cell") {
      schema_dofs_ += OPERATOR_SCHEMA_DOFS_CELL;
    } else if (names[i] == "node") {
      schema_dofs_ += OPERATOR_SCHEMA_DOFS_NODE;
    } else if (names[i] == "face") {
      schema_dofs_ += OPERATOR_SCHEMA_DOFS_FACE;
    }
  }

  // define stencil for preconditionre
  if (plist.isParameter("preconditioner schema")) {
    names = plist.get<Teuchos::Array<std::string> > ("preconditioner schema").toVector();

    schema_prec_dofs_ = 0;
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
    schema_prec_dofs_ = schema_dofs_;
  }

  special_assembling_ = 0;
  if (schema_prec_dofs_ != schema_dofs_) {
    if (schema_prec_dofs_ == OPERATOR_SCHEMA_DOFS_FACE) special_assembling_ = 1;
    if (schema_prec_dofs_ == OPERATOR_SCHEMA_DOFS_CELL) special_assembling_ = 2;
  }

  // Define base for assembling.
  std::string primary = plist.get<std::string>("discretization primary");
  std::string secondary = plist.get<std::string>("discretization secondary");

  schema_base_ = OPERATOR_SCHEMA_BASE_CELL;
  if (primary == "fv: default") {
    schema_base_ = OPERATOR_SCHEMA_BASE_FACE;
  }
  if (primary == "mfd: two-point flux approximation" && schema_dofs_ == OPERATOR_SCHEMA_DOFS_CELL) {
    schema_base_ = OPERATOR_SCHEMA_BASE_FACE;
  }

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
  } else if (primary == "fv: default") {
    mfd_primary_ = -1;  // not a mfd scheme
  } else {
    Errors::Message msg("OperatorDiffusion: primary discretization method is not supported.");
    Exceptions::amanzi_throw(msg);
  }

  // Secondary discretization methods
  if (secondary == "mfd: two-point flux approximation") {
    mfd_secondary_ = WhetStone::DIFFUSION_TPFA;
  } else if (secondary == "mfd: optimized for sparsity") {
    mfd_secondary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY;
  } else if (secondary == "mfd: support operator") {
    mfd_secondary_ = WhetStone::DIFFUSION_SUPPORT_OPERATOR;
  } else if (primary == "mfd: default") {
    mfd_primary_ = WhetStone::DIFFUSION_POLYHEDRA_SCALED;
  } else if (primary == "fv: default") {
    mfd_primary_ = -1;  // not a mfd schems
  } else {
    Errors::Message msg("OperatorDiffusion: secondary discretization method is not supported.");
    Exceptions::amanzi_throw(msg);
  }

  // Define other parameters.
  schema_ = schema_base_ + schema_dofs_;
  factor_ = 1.0;

  // upwind options
  std::string name = plist.get<std::string>("upwind method", "none");
  if (name == "standard") {
    upwind_ = OPERATOR_UPWIND_FLUX;
  } else if (name == "amanzi: artificial diffusion") {  
    upwind_ = OPERATOR_UPWIND_AMANZI_ARTIFICIAL_DIFFUSION;
  } else if (name == "amanzi: mfd") {  
    upwind_ = OPERATOR_UPWIND_AMANZI_MFD;
  } else if (name == "none") {
    upwind_ = OPERATOR_UPWIND_NONE;  // cell-centered scheme.
  } else {
    ASSERT(false);
  }

  // experimental options
  nonstandard_symbolic_ = plist.get<int>("nonstandard symbolic assembling", 0);
}

}  // namespace Operators
}  // namespace Amanzi
