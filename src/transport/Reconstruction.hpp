/*
This class operates with a scalar cell-centered field.
Lisence: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
  Reconstruction scalar_cell_function(u);
  scalar_cell_function.calculateCellGradient();
*/

#ifndef __Reconstruction_hpp__
#define __Reconstruction_hpp__

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

class Reconstruction {  
 public:
  Reconstruction() { status = RECONSTRUCTION_NULL; }
  Reconstruction(Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh, 
                 Teuchos::RCP<Epetra_Vector> scalar_field) { 
    mesh_ = mesh;
    scalar_field_ = scalar_field;
  }
  ~Reconstruction() {};

  void Init();
  void reset_field(Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh,
                   Teuchos::RCP<Epetra_Vector> scalar_field) { 
    mesh_ = mesh;
    scalar_field_ = scalar_field;
    if (status == RECONSTRUCTION_NULL) Init();
  }

  void calculateCellGradient();
  void applyLimiter(Teuchos::RCP<Epetra_Vector>& limiter);
 
  // access
  double getValue(const int cell, const AmanziGeometry::Point& p);
  double getValue(AmanziGeometry::Point& gradient, const int cell, const AmanziGeometry::Point& p);

  Teuchos::RCP<Epetra_MultiVector> get_gradient() { return gradient_; }
  std::vector<double>& get_field_local_min() { return field_local_min_; }
  std::vector<double>& get_field_local_max() { return field_local_max_; }

 private:
  void populateLeastSquareSystem(AmanziGeometry::Point& centroid, 
                                 double field_value, 
                                 Teuchos::SerialDenseMatrix<int, double>& matrix, 
                                 double* rhs);
  void printLeastSquareSystem(Teuchos::SerialDenseMatrix<int, double>matrix, double* rhs);

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;

  Teuchos::RCP<Epetra_Vector> scalar_field_;  // scalar cell-centered field
  Teuchos::RCP<Epetra_MultiVector> gradient_;

  std::vector<double> field_local_min_;
  std::vector<double> field_local_max_;

  int cmin, cmax_owned, cmax, number_owned_cells, number_wghost_cells;
  int fmin, fmax_owned, fmax, number_owned_faces, number_wghost_faces;

  Teuchos::SerialDenseMatrix<int, double> matrix;

  int status, dim;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif
