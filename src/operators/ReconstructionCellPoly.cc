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
#include "OperatorDefs.hh"
#include "Point.hh"
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

  // process other parameters
  degree_ = plist.get<int>("polynomial order", 0);
}


/* ******************************************************************
* Gradient of linear reconstruction is based on stabilized 
* least-square fit.
****************************************************************** */
void ReconstructionCellPoly::ComputePoly(
    const AmanziMesh::Entity_ID_List& ids,
    const Teuchos::RCP<const Epetra_MultiVector>& field, int component,
    const Teuchos::RCP<const BCs>& bc)
{
  field_ = field;
  component_ = component;

  Epetra_MultiVector& poly = *poly_->ViewComponent("cell");
  std::set<AmanziMesh::Entity_ID> cells, faces;
  AmanziGeometry::Point xcc(d_);

  int npoly = poly.NumVectors();
  WhetStone::DenseMatrix matrix(npoly, npoly);
  WhetStone::DenseVector rhs(npoly), coef(npoly);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    // mesh_->cell_get_face_adj_cells(c, AmanziMesh::Parallel_type::ALL, &cells);
    CellAllAdjCells_(c, AmanziMesh::Parallel_type::ALL, cells);

    matrix.PutScalar(0.0);
    rhs.PutScalar(0.0);

    // -- add neighboring cells
    for (int c1 : cells) {
      const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c1);
      int k(d_);
      xcc = xc2 - xc;

      for (int i = 0; i < d_; ++i) coef(i) = xcc[i];  // gradient
      for (int i = 0; i < d_; ++i)  // Hessian
        for (int j = i; j < d_; ++j) coef(k++) = xcc[i] * xcc[j];

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
          int k(d_);
          xcc = xf - xc;

          for (int i = 0; i < d_; ++i) coef(i) = xcc[i];  // gradient
          for (int i = 0; i < d_; ++i)  // Hessian
            for (int j = i; j < d_; ++j) coef(k++) = xcc[i] * xcc[j];

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

    // rhs[0] = rhs[1] = rhs[2] = 0.0;  // TESTING COMPATABILITY
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
  Epetra_MultiVector poly = *poly_->ViewComponent("cell");
  AmanziGeometry::Point xcc = p - mesh_->cell_centroid(c);

  int k(d_);
  double value = (*field_)[component_][c];
  for (int i = 0; i < d_; i++) value += poly[i][c] * xcc[i];
  for (int i = 0; i < d_; ++i)
    for (int j = i; j < d_; ++j) value += poly[k++][c] * xcc[i] * xcc[j];
  return value;
}

}  // namespace Operator
}  // namespace Amanzi

