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
#include "ReconstructionCell.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of basic parameters.
* NOTE: we assume that ghost values of field were already populated.
****************************************************************** */
void ReconstructionCell::Init(Teuchos::ParameterList& plist)
{
  dim = mesh_->space_dimension();

  CompositeVectorSpace cv_space;
  cv_space.SetMesh(mesh_);
  cv_space.SetGhosted(true);
  cv_space.SetComponent("cell", AmanziMesh::CELL, dim);

  gradient_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cv_space, true));

  // process other parameters
  poly_order_ = plist.get<int>("polynomial order", 0);
}


/* ******************************************************************
* Gradient of linear reconstruction is based on stabilized 
* least-square fit.
****************************************************************** */
void ReconstructionCell::ComputeGradient(
    const AmanziMesh::Entity_ID_List& ids,
    const Teuchos::RCP<const Epetra_MultiVector>& field, int component)
{
  field_ = field;
  component_ = component;

  Epetra_MultiVector& grad = *gradient_->ViewComponent("cell", false);
  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point xcc(dim);

  WhetStone::DenseMatrix matrix(dim, dim);
  WhetStone::DenseVector rhs(dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    // mesh_->cell_get_face_adj_cells(c, AmanziMesh::Parallel_type::ALL, &cells);
    CellFaceAdjCellsNonManifold_(c, AmanziMesh::Parallel_type::ALL, cells);
    int ncells = cells.size();

    matrix.PutScalar(0.0);
    rhs.PutScalar(0.0);

    for (int n = 0; n < ncells; n++) {
      const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(cells[n]);
      for (int i = 0; i < dim; i++) xcc[i] = xc2[i] - xc[i];

      double value = (*field_)[component_][cells[n]] - (*field_)[component_][c];
      PopulateLeastSquareSystem_(xcc, value, matrix, rhs);
    }

    // improve robustness w.r.t degenerate matrices
    double det = matrix.Det();
    double norm = matrix.NormInf();

    if (det < pow(norm, 1.0/dim)) {
      norm *= OPERATOR_RECONSTRUCTION_MATRIX_CORRECTION;
      for (int i = 0; i < dim; i++) matrix(i, i) += norm;
    }

    int info, nrhs = 1;
    WhetStone::DPOSV_F77("U", &dim, &nrhs, matrix.Values(), &dim, rhs.Values(), &dim, &info);
    if (info) {  // reduce reconstruction order
      for (int i = 0; i < dim; i++) rhs(i) = 0.0;
    }

    // rhs[0] = rhs[1] = rhs[2] = 0.0;  // TESTING COMPATABILITY
    for (int i = 0; i < dim; i++) grad[i][c] = rhs(i);
  }

  gradient_->ScatterMasterToGhosted("cell");
}


/* ******************************************************************
* Assemble a SPD least square matrix
****************************************************************** */
void ReconstructionCell::PopulateLeastSquareSystem_(
    AmanziGeometry::Point& centroid, double field_value,
    WhetStone::DenseMatrix& matrix, WhetStone::DenseVector& rhs)
{
  for (int i = 0; i < dim; i++) {
    double xyz = centroid[i];

    matrix(i, i) += xyz * xyz;
    for (int j = i + 1; j < dim; j++) matrix(j, i) = matrix(i, j) += xyz * centroid[j];

    rhs(i) += xyz * field_value;
  }
}


/* ******************************************************************
* On intersecting manifolds, we extract neighboors living in the same 
* manifold using a smoothness criterion.
****************************************************************** */
void ReconstructionCell::CellFaceAdjCellsNonManifold_(
    AmanziMesh::Entity_ID c, AmanziMesh::Parallel_type ptype,
    std::vector<AmanziMesh::Entity_ID>& cells) const
{
  AmanziMesh::Entity_ID_List faces, fcells;

  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  cells.clear();

  for (int n = 0; n < nfaces; ++n) {
    AmanziMesh::Entity_ID f = faces[n];
    mesh_->face_get_cells(f, ptype, &fcells);
    int ncells = fcells.size();

    if (ncells == 2) {
      cells.push_back(fcells[0] + fcells[1] - c);
    } else if (ncells > 2) {
      int dir;
      double dmax(0.0);
      AmanziMesh::Entity_ID cmax;
      const AmanziGeometry::Point& normal0 = mesh_->face_normal(f, false, c, &dir);

      for (int i = 0; i < ncells; ++i) {
        AmanziMesh::Entity_ID c1 = fcells[i];
        if (c1 != c) {
          const AmanziGeometry::Point& normal1 = mesh_->face_normal(f, false, c1, &dir);
          double d = fabs(normal0 * normal1) / norm(normal1);
          if (d > dmax) {
            dmax = d;
            cmax = c1; 
          }
        }
      } 

      cells.push_back(cmax);
    }
  }
}


/* ******************************************************************
* Calculates reconstructed value at point p.
****************************************************************** */
double ReconstructionCell::getValue(int c, const AmanziGeometry::Point& p)
{
  Teuchos::RCP<Epetra_MultiVector> grad = gradient_->ViewComponent("cell", false);
  const auto& xc = mesh_->cell_centroid(c);

  double value = (*field_)[component_][c];
  for (int i = 0; i < dim; i++) value += (*grad)[i][c] * (p[i] - xc[i]);
  return value;
}


/* ******************************************************************
* Calculates reconstructed value at point p using external gradient.
****************************************************************** */
double ReconstructionCell::getValue(
    const AmanziGeometry::Point& gradient, int c, const AmanziGeometry::Point& p)
{
  const auto& xc = mesh_->cell_centroid(c);
  return (*field_)[component_][c] + gradient * (p - xc);
}

}  // namespace Operator
}  // namespace Amanzi

