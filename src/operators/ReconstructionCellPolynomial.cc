/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

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
#include "ReconstructionCellPolynomial.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of basic parameters.
* NOTE: we assume that ghost values of field were already populated.
****************************************************************** */
void
ReconstructionCellPolynomial::Init(Teuchos::ParameterList& plist)
{
  d_ = mesh_->getSpaceDimension();

  CompositeVectorSpace cvs, cvs2;
  cvs.SetMesh(mesh_)->SetGhosted(true)->SetComponent(
    "cell", AmanziMesh::Entity_kind::CELL, d_ * (d_ + 1) / 2 + d_);

  poly_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
  poly_c_ = poly_->ViewComponent("cell", true);

  cvs2.SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::Entity_kind::CELL, d_ * (d_ + 1) / 2);
  ortho_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs2, true));
  ortho_c_ = ortho_->ViewComponent("cell", true);

  // process other parameters
  degree_ = plist.get<int>("polynomial order", 2);
  AMANZI_ASSERT(degree_ == 2);
}


/* ******************************************************************
* Gradient of linear reconstruction is based on stabilized
* least-square fit.
****************************************************************** */
void
ReconstructionCellPolynomial::Compute(const AmanziMesh::Entity_ID_View& ids,
                                      const Teuchos::RCP<const Epetra_MultiVector>& field,
                                      int component,
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

  int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  WhetStone::Polynomial quad(d_, 2);
  std::vector<const WhetStone::WhetStoneFunction*> funcs(1);
  funcs[0] = &quad;
  WhetStone::NumericalIntegration numi(mesh_);

  for (int c = 0; c < ncells_owned; c++) {
    double vol = mesh_->getCellVolume(c);
    const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
    quad.set_origin(xc);

    CellAllAdjCells_(c, AmanziMesh::Parallel_kind::ALL, cells);

    // compute othogonality coefficients
    int n(0);
    for (int i = 0; i < d_; ++i) {
      for (int j = i; j < d_; ++j, ++n) {
        quad(2, n) = 1.0;
        (*ortho_c_)[n][c] = numi.IntegrateFunctionsTriangulatedCell(c, funcs, 2) / vol;
        quad(2, n) = 0.0;
      }
    }

    // populate least-squate matrix
    matrix.PutScalar(0.0);
    rhs.PutScalar(0.0);

    for (int c1 : cells) {
      double vol1 = mesh_->getCellVolume(c1);
      const AmanziGeometry::Point& xc1 = mesh_->getCellCentroid(c1);
      xcc = xc1 - xc;

      // -- gradient terms
      for (int i = 0; i < d_; ++i) coef(i) = xcc[i];

      // -- Hessian terms
      n = 0;
      for (int i = 0; i < d_; ++i) {
        for (int j = i; j < d_; ++j, ++n) {
          quad(0) = -(*ortho_c_)[n][c];
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
          const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
          xcc = xf - xc;

          // -- gradient
          for (int i = 0; i < d_; ++i) coef(i) = xcc[i];

          // -- Hessian terms
          n = 0;
          for (int i = 0; i < d_; ++i) {
            for (int j = i; j < d_; ++j, ++n) coef(d_ + n) = xcc[i] * xcc[j] - (*ortho_c_)[n][c];
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

      WhetStone::DPOSV_F77(
        "U", &npoly, &nrhs, matrix_copy.Values(), &npoly, rhs_copy.Values(), &npoly, &info);
      if (info) { // reduce reconstruction order
        for (int i = 0; i < npoly; i++) rhs(i) = 0.0;
      }
    }

    for (int i = 0; i < npoly; i++) poly[i][c] = rhs(i);
  }

  poly_->ScatterMasterToGhosted();
  ortho_->ScatterMasterToGhosted();
}


/* ******************************************************************
* Assemble a SPD least square matrix
****************************************************************** */
void
ReconstructionCellPolynomial::PopulateLeastSquareSystem_(WhetStone::DenseVector& coef,
                                                         double field_value,
                                                         WhetStone::DenseMatrix& matrix,
                                                         WhetStone::DenseVector& rhs)
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
void
ReconstructionCellPolynomial::CellAllAdjCells_(AmanziMesh::Entity_ID c,
                                               AmanziMesh::Parallel_kind ptype,
                                               std::set<AmanziMesh::Entity_ID>& cells) const
{
  std::set<int> ids;
  AmanziMesh::Entity_ID_View nodes, vcells;

  cells.clear();

  nodes = mesh_->getCellNodes(c);
  for (int i = 0; i < nodes.size(); i++) {
    int v = nodes[i];
    vcells = mesh_->getNodeCells(v, AmanziMesh::Parallel_kind::ALL);

    for (int k = 0; k < vcells.size(); ++k) {
      int c1 = vcells[k];
      if (cells.find(c1) == cells.end() && c != c1) cells.insert(c1);
    }
  }
}


void
ReconstructionCellPolynomial::CellAllAdjFaces_(AmanziMesh::Entity_ID c,
                                               const std::set<AmanziMesh::Entity_ID>& cells,
                                               std::set<AmanziMesh::Entity_ID>& faces) const
{
  AmanziMesh::Entity_ID_View cfaces, fcells;

  faces.clear();
  for (int c1 : cells) {
    cfaces = mesh_->getCellFaces(c1);
    for (int f : cfaces) {
      fcells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      if (fcells.size() == 1 && faces.find(f) == faces.end()) faces.insert(f);
    }
  }
}


/* ******************************************************************
* Calculates reconstructed value at point p.
****************************************************************** */
double
ReconstructionCellPolynomial::getValue(int c, const AmanziGeometry::Point& p)
{
  AmanziGeometry::Point xcc = p - mesh_->getCellCentroid(c);

  double value = (*field_)[0][c];
  for (int i = 0; i < d_; i++) value += (*poly_c_)[i][c] * xcc[i];

  int n(0);
  for (int i = 0; i < d_; ++i) {
    for (int j = i; j < d_; ++j, ++n) {
      value += (*poly_c_)[d_ + n][c] * (xcc[i] * xcc[j] - (*ortho_c_)[n][c]);
    }
  }

  return value;
}


/* ******************************************************************
* Calculates deviation from mean value at point p.
****************************************************************** */
double
ReconstructionCellPolynomial::getValueSlope(int c, const AmanziGeometry::Point& p)
{
  AmanziGeometry::Point xcc = p - mesh_->getCellCentroid(c);

  double value(0.0);
  for (int i = 0; i < d_; i++) value += (*poly_c_)[i][c] * xcc[i];

  int n(0);
  for (int i = 0; i < d_; ++i) {
    for (int j = i; j < d_; ++j, ++n) {
      value += (*poly_c_)[d_ + n][c] * (xcc[i] * xcc[j] - (*ortho_c_)[n][c]);
    }
  }

  return value;
}


/* ******************************************************************
* Returns full polynomial
****************************************************************** */
WhetStone::Polynomial
ReconstructionCellPolynomial::getPolynomial(int c) const
{
  WhetStone::Polynomial tmp(d_, 2);
  tmp(0) = (*field_)[0][c];
  for (int i = 0; i < d_; i++) tmp(i + 1) = (*poly_c_)[i][c];

  int n(0);
  for (int i = 0; i < d_; ++i) {
    for (int j = i; j < d_; ++j, ++n) {
      tmp(d_ + n + 1) = (*poly_c_)[d_ + n][c];
      tmp(0) -= (*poly_c_)[d_ + n][c] * (*ortho_c_)[n][c];
    }
  }

  tmp.set_origin(mesh_->getCellCentroid(c));
  return tmp;
}

} // namespace Operators
} // namespace Amanzi
