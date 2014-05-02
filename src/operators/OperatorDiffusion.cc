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
#include "mfd3d_diffusion.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusion.hh"


namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator.                                           
****************************************************************** */
void OperatorDiffusion::InitOperator(
    std::vector<WhetStone::Tensor>& K, Teuchos::RCP<NonlinearCoefficient> k,
    int schema_base, int schema_dofs, const Teuchos::ParameterList& plist)
{
  plist_ = plist; 
  k_ = k;
  schema_base_ = schema_base;
  schema_dofs_ = schema_dofs;
  schema_ = schema_base_ + schema_dofs_;

  if (schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL) {
    CreateMassMatrices_(K);
  }

  // if upwind is requested, we will need to update nonlinear coefficient 
  std::string str_upwind = plist_.get<std::string>("upwind method", "amanzi");
}


/* ******************************************************************
* Basic routine of each operator: creation of matrices.
****************************************************************** */
void OperatorDiffusion::UpdateMatrices(Teuchos::RCP<const CompositeVector> flux)
{
  // find location of matrix blocks
  int m(0), nblocks = blocks_.size();
  bool flag(false);

  for (int n = 0; n < nblocks; n++) {
    int schema = blocks_schema_[n];
    if ((schema & OPERATOR_SCHEMA_DOFS_CELL) && (schema & OPERATOR_SCHEMA_DOFS_FACE)) {
      m = n;
      flag = true;
      break;
    }
  }

  if (flag == false) { 
    m = nblocks++;
    blocks_schema_.push_back(OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_FACE + OPERATOR_SCHEMA_DOFS_CELL);
    blocks_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
    blocks_shadow_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
  }
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[m];
  WhetStone::DenseMatrix null_matrix;

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
    if (k_ != Teuchos::null) {
      kc = (*k_->cvalues())[c];
    }

    double matsum = 0.0;  // elimination of mass matrix
    for (int n = 0; n < nfaces; n++) {
      double rowsum = 0.0;
      for (int m = 0; m < nfaces; m++) {
        double tmp = Wff(n, m) * kc;
        rowsum += tmp;
        Acell(n, m) = tmp;
      }

      Acell(n, nfaces) = -rowsum;
      Acell(nfaces, n) = -rowsum;
      matsum += rowsum;
    }
    Acell(nfaces, nfaces) = matsum;


    // Update terms due to dependence of k on the solution.
    if (flux !=  Teuchos::null && k_ != Teuchos::null) {
      const Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        double dkf = (*k_->fderivatives())[f];
        double  kf = (*k_->fvalues())[f];
        double alpha = (dkf / kf) * flux_data[0][f] * dirs[n];
        if (alpha > 0) {
          Acell(n, n) += kc * alpha;
        }
      }
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
void OperatorDiffusion::UpdateMatricesStiffness(std::vector<WhetStone::Tensor>& K)
{
  // find location of matrix blocks
  int m(0), nblocks = blocks_.size();
  bool flag(false);

  for (int nb = 0; nb < nblocks; nb++) {
    int schema = blocks_schema_[nb];
    if (schema == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_NODE) {
      m = nb;
      flag = true;
      break;
    }
  }

  if (flag == false) { 
    m = nblocks++;
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
  AmanziMesh::Entity_ID_List nodes;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    WhetStone::DenseMatrix Acell(nnodes, nnodes);
    int ok = mfd.StiffnessMatrix(c, K[c], Acell);

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
* Calculate elemental inverse mass matrices.                                           
****************************************************************** */
void OperatorDiffusion::CreateMassMatrices_(std::vector<WhetStone::Tensor>& K)
{
  WhetStone::MFD3D_Diffusion mfd(mesh_);

  AmanziMesh::Entity_ID_List faces;

  Wff_cells_.clear();

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix Wff(nfaces, nfaces);
    int ok = mfd.MassMatrixInverse(c, K[c], Wff);

    Wff_cells_.push_back(Wff);

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("OperatorDiffusionSurface: unexpected failure in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi

