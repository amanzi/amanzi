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
#include "ReconstructionCellLinear.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of basic parameters.
* NOTE: we assume that ghost values of field were already populated.
****************************************************************** */
void ReconstructionCellLinear::Init(Teuchos::ParameterList& plist)
{
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, dim);

  gradient_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
  gradient_c_ = gradient_->ViewComponent("cell", true);
  if (plist.get<std::string>("weight", "constant") == "inverse distance") weight_ = WeightType::WT_INVERSE_DISTANCE;
}


/* ******************************************************************
* Gradient of linear reconstruction is based on stabilized 
* least-square fit.
****************************************************************** */
void ReconstructionCellLinear::Compute(
    const AmanziMesh::Entity_ID_List& ids,
    const Teuchos::RCP<const Epetra_MultiVector>& field, int component,
    const Teuchos::RCP<const BCs>& bc)
{
  field_ = field;
  component_ = component;

  Epetra_MultiVector& grad = *gradient_->ViewComponent("cell");
  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point xcc(dim);

  WhetStone::DenseMatrix matrix(dim, dim);
  WhetStone::DenseVector rhs(dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    // mesh_->cell_get_face_adj_cells(c, AmanziMesh::Parallel_type::ALL, &cells);
    CellFaceAdjCellsManifold_(c, AmanziMesh::Parallel_type::ALL, cells);
    int ncells = cells.size();

    matrix.PutScalar(0.0);
    rhs.PutScalar(0.0);

    for (int n = 0; n < ncells; n++) {
      const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(cells[n]);
      xcc = xc2 - xc;

      double value = (*field_)[component_][cells[n]] - (*field_)[component_][c];
      PopulateLeastSquareSystem_(xcc, value, matrix, rhs);
    }

    // improve robustness w.r.t degenerate matrices
    double det = matrix.Det();
    double norm = matrix.NormInf();

    if (det < pow(norm, 1.0/dim)) {
      norm *= OPERATOR_RECONSTRUCTION_MATRIX_CORRECTION;  // relative 
      norm += OPERATOR_RECONSTRUCTION_MATRIX_CORRECTION;  // absolute
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
void ReconstructionCellLinear::PopulateLeastSquareSystem_(
    AmanziGeometry::Point& centroid, double field_value,
    WhetStone::DenseMatrix& matrix, WhetStone::DenseVector& rhs)
{
  double w = (weight_ == WeightType::WT_INVERSE_DISTANCE) ? 1.0 / norm(centroid) : 1.0; 

  for (int i = 0; i < dim; i++) {
    double xyz = centroid[i];
    double xyz_w = xyz * w;

    matrix(i, i) += xyz * xyz_w;
    for (int j = i + 1; j < dim; j++) matrix(j, i) = matrix(i, j) += xyz_w * centroid[j];

    rhs(i) += xyz_w * field_value;
  }
}


/* ******************************************************************
* On intersecting manifolds, we extract neighboors living in the same 
* manifold using a smoothness criterion.
****************************************************************** */
void ReconstructionCellLinear::CellFaceAdjCellsManifold_(
    AmanziMesh::Entity_ID c, AmanziMesh::Parallel_type ptype,
    std::vector<AmanziMesh::Entity_ID>& cells) const
{
  AmanziMesh::Entity_ID_List fcells;

  const auto& faces = mesh_->cell_get_faces(c);
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
double ReconstructionCellLinear::getValue(
    int c, const AmanziGeometry::Point& p)
{
  const auto& xc = mesh_->cell_centroid(c);

  double value = (*field_)[component_][c];
  for (int i = 0; i < dim; i++) value += (*gradient_c_)[i][c] * (p[i] - xc[i]);
  return value;
}


/* ******************************************************************
* Calculates deviation from a mean value at point p.
****************************************************************** */
double ReconstructionCellLinear::getValueSlope(
    int c, const AmanziGeometry::Point& p)
{
  const auto& xc = mesh_->cell_centroid(c);

  double value(0.0);
  for (int i = 0; i < dim; i++) value += (*gradient_c_)[i][c] * (p[i] - xc[i]);
  return value;
}


/* ******************************************************************
* Returns full polynomial
****************************************************************** */
WhetStone::Polynomial ReconstructionCellLinear::getPolynomial(int c) const
{
  WhetStone::Polynomial tmp(dim, 2);
  tmp(0) = (*field_)[0][c];
  for (int i = 0; i < dim; i++) tmp(i + 1) = (*gradient_c_)[i][c];

  tmp.set_origin(mesh_->cell_centroid(c));
  return tmp;
}

}  // namespace Operator
}  // namespace Amanzi

