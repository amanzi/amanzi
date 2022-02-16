/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include "Epetra_Vector.h"

#include "DenseMatrix.hh"
#include "lapack.hh"
#include "Mesh.hh"
#include "NumericalIntegration.hh"
#include "OperatorDefs.hh"
#include "Point.hh"
#include "Polynomial.hh"
#include "ReconstructionCellPoly.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of basic parameters.
* NOTE: we assume that ghost values of field were already populated.
****************************************************************** */
void ReconstructionCellPoly::Init(Teuchos::ParameterList& plist)
{
  d_ = mesh_->space_dimension();

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, d_ * (d_ + 1) / 2 + d_);

  poly_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
  poly_c_ = poly_->ViewComponent("cell");

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ortho_.resize(ncells);

  // process other parameters
  degree_ = plist.get<int>("polynomial order", 0);
}


/* ******************************************************************
* Gradient of linear reconstruction is based on stabilized 
* least-square fit.
****************************************************************** */
void ReconstructionCellPoly::Compute(
    const AmanziMesh::Entity_ID_List& ids,
    const Teuchos::RCP<const Epetra_MultiVector>& field, int component,
    const Teuchos::RCP<const BCs>& bc)
{
  field_ = field;
  component_ = component;

  auto& poly = *poly_->ViewComponent("cell");
  std::set<AmanziMesh::Entity_ID> cells, faces;
  AmanziGeometry::Point xcc(d_);

  int npoly = poly.NumVectors();
  WhetStone::DenseMatrix matrix(npoly, npoly);
  WhetStone::DenseVector rhs(npoly), coef(npoly);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  WhetStone::Polynomial quad(d_, 2);
  std::vector<const WhetStone::WhetStoneFunction*> funcs(1);
  funcs[0] = &quad;
  WhetStone::NumericalIntegration numi(mesh_);

  for (int c = 0; c < ncells_owned; c++) {
    double vol = mesh_->cell_volume(c);
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    quad.set_origin(xc);

    CellAllAdjCells_(c, AmanziMesh::Parallel_type::ALL, cells);

    // compute othogonality coefficients
    ortho_[c].clear();
    int n(0);
    for (int i = 0; i < d_; ++i) {
      for (int j = i; j < d_; ++j, ++n) {
        quad(2, n) = 1.0;
        ortho_[c].push_back(numi.IntegrateFunctionsTriangulatedCell(c, funcs, 2) / vol);
        quad(2, n) = 0.0;
      }
    }

    // populate least-squate matrix
    matrix.PutScalar(0.0);
    rhs.PutScalar(0.0);

    for (int c1 : cells) {
      double vol1 = mesh_->cell_volume(c1);
      const AmanziGeometry::Point& xc1 = mesh_->cell_centroid(c1);
      xcc = xc1 - xc;

      // -- gradient terms
      for (int i = 0; i < d_; ++i) coef(i) = xcc[i];

      // -- Hessian terms
      int n(0);
      for (int i = 0; i < d_; ++i) {
        for (int j = i; j < d_; ++j, ++n) {
          quad(0) = -ortho_[c][n];
          quad(2, n) = 1.0;
          coef(d_ + n) = numi.IntegrateFunctionsTriangulatedCell(c1, funcs, 2) / vol1;
          quad(2, n) = 0.0;
        }
      }
      quad(0) = 0.0;

      double value = (*field_)[component_][c1] - (*field_)[component_][c];
      PopulateLeastSquareSystem_(coef, value, matrix, rhs);
    }

    // -- add boundary faces
    if (bc.get()) {
      const auto& bc_model = bc->bc_model(); 
      const auto& bc_value = bc->bc_value(); 

      cells.insert(c);
      CellAllAdjFaces_(c, cells, faces);
      for (int f : faces) {
        if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
          const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
          xcc = xf - xc;

          // -- gradient
          for (int i = 0; i < d_; ++i) coef(i) = xcc[i];

          // -- Hessian terms
          int n(0);
          for (int i = 0; i < d_; ++i) {
            for (int j = i; j < d_; ++j, ++n) coef(d_ + n) = xcc[i] * xcc[j] - ortho_[c][n];
          }

          double value = bc_value[f] - (*field_)[component_][c];
          PopulateLeastSquareSystem_(coef, value, matrix, rhs);
        }
      }
    }

    // improve robustness w.r.t degenerate matrices
    WhetStone::DenseVector rhs_copy(rhs);
    WhetStone::DenseMatrix matrix_copy(matrix);

    int info, nrhs = 1;
    WhetStone::DPOSV_F77("U", &npoly, &nrhs, matrix.Values(), &npoly, rhs.Values(), &npoly, &info);
    if (info) {
      // try regularized matrix 
      double norm = matrix_copy.NormInf() * OPERATOR_RECONSTRUCTION_MATRIX_CORRECTION;
      for (int i = 0; i < npoly; i++) matrix_copy(i, i) += norm;

      WhetStone::DPOSV_F77("U", &npoly, &nrhs, matrix_copy.Values(), &npoly, rhs_copy.Values(), &npoly, &info);
      if (info) {  // reduce reconstruction order
        for (int i = 0; i < npoly; i++) rhs(i) = 0.0;
      }
    }

    for (int i = 0; i < npoly; i++) poly[i][c] = rhs(i);
  }

  poly_->ScatterMasterToGhosted("cell");
}


/* ******************************************************************
* Assemble a SPD least square matrix
****************************************************************** */
void ReconstructionCellPoly::PopulateLeastSquareSystem_(
    WhetStone::DenseVector& coef, double field_value,
    WhetStone::DenseMatrix& matrix, WhetStone::DenseVector& rhs)
{
  int npoly = matrix.NumRows();
  for (int i = 0; i < npoly; i++) {
    for (int j = i; j < npoly; j++) matrix(j, i) = matrix(i, j) += coef(i) * coef(j);

    rhs(i) += coef(i) * field_value;
  }
}


/* ******************************************************************
* On intersecting manifolds, we extract neighboors living in the same 
* manifold using a smoothness criterion.
****************************************************************** */
void ReconstructionCellPoly::CellAllAdjCells_(
    AmanziMesh::Entity_ID c, AmanziMesh::Parallel_type ptype,
    std::set<AmanziMesh::Entity_ID>& cells) const
{
  std::set<int> ids;
  AmanziMesh::Entity_ID_List nodes, vcells;

  cells.clear();

  mesh_->cell_get_nodes(c, &nodes);
  for (int i = 0; i < nodes.size(); i++) {
    int v = nodes[i];
    mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &vcells);

    for (int k = 0; k < vcells.size(); ++k) {
      int c1 = vcells[k];
      if (cells.find(c1) == cells.end() && c != c1) cells.insert(c1);
    }
  }
}


void ReconstructionCellPoly::CellAllAdjFaces_(
    AmanziMesh::Entity_ID c,
    const std::set<AmanziMesh::Entity_ID>& cells,
    std::set<AmanziMesh::Entity_ID>& faces) const
{
  AmanziMesh::Entity_ID_List cfaces, fcells;

  faces.clear();
  for (int c : cells) {
    mesh_->cell_get_faces(c, &cfaces);
    for (int f : cfaces) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &fcells);
      if (fcells.size() == 1 && faces.find(f) == faces.end()) faces.insert(f);
    }
  }
}


/* ******************************************************************
* Calculates reconstructed value at point p.
****************************************************************** */
double ReconstructionCellPoly::getValue(int c, const AmanziGeometry::Point& p)
{
  AmanziGeometry::Point xcc = p - mesh_->cell_centroid(c);

  double value = (*field_)[0][c];
  for (int i = 0; i < d_; i++) value += (*poly_c_)[i][c] * xcc[i];

  int n(0);
  for (int i = 0; i < d_; ++i) {
    for (int j = i; j < d_; ++j, ++n) {
      value += (*poly_c_)[d_ + n][c] * (xcc[i] * xcc[j] - ortho_[c][n]);
    }
  }

  return value;
}


/* ******************************************************************
* Calculates deviation from mean value at point p.
****************************************************************** */
double ReconstructionCellPoly::getValueSlope(int c, const AmanziGeometry::Point& p)
{
  AmanziGeometry::Point xcc = p - mesh_->cell_centroid(c);

  double value(0.0);
  for (int i = 0; i < d_; i++) value += (*poly_c_)[i][c] * xcc[i];

  int n(0);
  for (int i = 0; i < d_; ++i) {
    for (int j = i; j < d_; ++j, ++n) {
      value += (*poly_c_)[d_ + n][c] * (xcc[i] * xcc[j] - ortho_[c][n]);
    }
  }

  return value;
}


/* ******************************************************************
* Returns full polynomial
****************************************************************** */
WhetStone::Polynomial ReconstructionCellPoly::getPolynomial(int c) 
{
  WhetStone::Polynomial tmp(d_, 2);
  tmp(0) = (*field_)[0][c];
  for (int i = 0; i < d_; i++) tmp(i + 1) = (*poly_c_)[i][c];

  int n(0);
  for (int i = 0; i < d_; ++i) {
    for (int j = i; j < d_; ++j, ++n) {
      tmp(d_ + n + 1) = (*poly_c_)[d_ + n][c];
      tmp(0) -= (*poly_c_)[d_ + n][c] * ortho_[c][n];
    }
  }

  tmp.set_origin(mesh_->cell_centroid(c));
  return tmp;
}

}  // namespace Operator
}  // namespace Amanzi

