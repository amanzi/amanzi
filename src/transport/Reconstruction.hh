/*
  This class operates with a scalar cell-centered field.

  Lisence: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Usage: 
    Reconstruction scalar_cell_function(u);
    scalar_cell_function.calculateCellGradient();
*/

#ifndef AMANZI_RECONSTRACTION_HH_
#define AMANZI_RECONSTRACTION_HH_

#include <vector>

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Mesh.hh"
#include "Point.hh"


namespace Amanzi {
namespace AmanziTransport {

const int RECONSTRUCTION_NULL = 0;
const int RECONSTRUCTION_INIT = 1;

const double RECONSTRUCTION_MATRIX_CORRECTION = 1e-15;

class Reconstruction {  
 public:
  Reconstruction() { status = RECONSTRUCTION_NULL; }
  Reconstruction(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, 
                 Teuchos::RCP<Epetra_Vector> scalar_field) :
    mesh_(mesh), scalar_field_(scalar_field) {}

  ~Reconstruction() {};

  void Init();
  void ResetField(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                  Teuchos::RCP<Epetra_Vector> scalar_field) { 
    mesh_ = mesh;
    scalar_field_ = scalar_field;
    if (status == RECONSTRUCTION_NULL) Init();
  }

  void CalculateCellGradient();
  void applyLimiter(Teuchos::RCP<Epetra_Vector>& limiter);
 
  // access
  double getValue(const int cell, const AmanziGeometry::Point& p);
  double getValue(AmanziGeometry::Point& gradient, const int cell, const AmanziGeometry::Point& p);

  Teuchos::RCP<CompositeVector> gradient() { return gradient_; }

 private:
  void populateLeastSquareSystem(AmanziGeometry::Point& centroid, 
                                 double field_value, 
                                 Teuchos::SerialDenseMatrix<int, double>& matrix, 
                                 double* rhs);
  double calculateMatrixNorm(Teuchos::SerialDenseMatrix<int, double>& matrix);
  double calculateMatrixDeterminant(Teuchos::SerialDenseMatrix<int, double>& matrix);
  int findMinimalDiagonalEntry(Teuchos::SerialDenseMatrix<int, double>& matrix);
  void printLeastSquareSystem(Teuchos::SerialDenseMatrix<int, double>matrix, double* rhs);

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  Teuchos::RCP<Epetra_Vector> scalar_field_;  // scalar cell-centered field
  Teuchos::RCP<CompositeVector> gradient_;

  int cmax_owned, cmax, number_owned_cells, number_wghost_cells;
  int fmin, fmax_owned, fmax, number_owned_faces, number_wghost_faces;

  Teuchos::SerialDenseMatrix<int, double> matrix;

  int status, dim;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif
