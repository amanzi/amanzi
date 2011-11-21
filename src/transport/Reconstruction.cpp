/*
The transport component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "dbc.hh"

#include "Transport_PK.hpp"
#include "Reconstruction.hpp"

namespace Amanzi {
namespace AmanziTransport {

/* We set up most popular parameters here. */
void Reconstruction::Init()
{
  status = RECONSTRUCTION_NULL;
  const Epetra_Map& cmap = mesh_->cell_map(true);
  const Epetra_Map& fmap = mesh_->face_map(true);

  cmin = cmap.MinLID();
  cmax = cmap.MaxLID();

  number_owned_cells = mesh_->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  cmax_owned = cmin + number_owned_cells - 1;

  fmin = fmap.MinLID();
  fmax = fmap.MaxLID(); 

  number_owned_faces = mesh_->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  fmax_owned = fmin + number_owned_faces - 1;

  double *memory;  // we need to allocate sufficient memory for LAPACK routines
  memory = new double(TRANSPORT_MAX_FACES * TRANSPORT_MAX_FACES);
  Teuchos::SerialDenseMatrix<int, double> matrix(Teuchos::View, 
                                                 memory, 
                                                 TRANSPORT_MAX_FACES,
                                                 TRANSPORT_MAX_FACES,
                                                 TRANSPORT_MAX_FACES);

  dim = mesh_->space_dimension();
  gradient_ = Teuchos::rcp(new Epetra_MultiVector(cmap, 3));

  field_local_min_.resize(cmax_owned+1);
  field_local_max_.resize(cmax_owned+1);

  status = RECONSTRUCTION_INIT;
}


/* ******************************************************************
* Implementation is tuned up for gradient (first-order reconstruction).
* It can be extended easily if needed in the future.
****************************************************************** */
void Reconstruction::calculateCellGradient()
{
  Epetra_Vector& u = *scalar_field_;  // a few aliases
  Teuchos::LAPACK<int, double> lapack;

  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point xc2(dim);
 
  double *rhs; 
  rhs = new double[dim]; 

  for (int c=cmin; c<=cmax_owned; c++) {
    matrix.shape(dim, dim);  // Teuchos will initilize this matrix by zeros
    for (int i=0; i<dim; i++) rhs[i] = 0.0;

    mesh_->cell_get_face_adj_cells(c, AmanziMesh::USED, &cells); 

    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    field_local_min_[c] = field_local_max_[c] = u[c];

    for (int i=0; i<cells.size(); i++) {
      xc2  = mesh_->cell_centroid(cells[i]);
      xc2 -= xc;

      double value = u[cells[i]] - u[c];
      populateLeastSquareSystem(xc2, value, matrix, rhs);

      value = u[cells[i]];   
      field_local_min_[c] = std::min(field_local_min_[c], value);
      field_local_max_[c] = std::max(field_local_max_[c], value);
    }

    // improve robustness w.r.t degenerate matrices
    double det = calculate_matrix_determinant(matrix);
    double norm = calculate_matrix_norm(matrix);

    if (det < pow(norm, 1.0/dim)) {
      norm *= RECONSTRUCTION_MATRIX_CORRECTION;
      for (int i=0; i<dim; i++) matrix(i,i) += norm;
    }
    //printLeastSquareSystem(matrix, rhs);
 
    int info;
    lapack.POSV('U', dim, 1, matrix.values(), dim, rhs, dim, &info); 
    if (info) {  // reduce reconstruction order
      for (int i=0; i<dim; i++) rhs[i] = 0.0; 
    }

    //rhs[0] = rhs[1] = rhs[2] = 0.0;  // TESTING COMPATABILITY 
    for (int i=0; i<dim; i++) (*gradient_)[i][c] = rhs[i];
  }

  delete [] rhs;

#ifdef HAVE_MPI
  const Epetra_BlockMap& source_fmap = (*gradient_).Map();
  const Epetra_BlockMap& target_fmap = (*gradient_).Map();

  Epetra_Import importer(target_fmap, source_fmap);
  (*gradient_).Import(*gradient_, importer, Insert);
#endif
}


/* ******************************************************************
 * The limiter must be between 0 and 1
****************************************************************** */
void Reconstruction::applyLimiter(Teuchos::RCP<Epetra_Vector>& limiter)
{
  for (int c=cmin; c<=cmax; c++) { 
    for (int i=0; i<dim; i++) (*gradient_)[i][c] *= (*limiter)[c];
  }
}


/* ******************************************************************
 * calculates a value at point p using gradinet and centroid
****************************************************************** */
double Reconstruction::getValue(const int cell, const AmanziGeometry::Point& p)
{
  AmanziGeometry::Point xc(dim);
  xc = p;
  xc -= mesh_->cell_centroid(cell); 

  double value = (*scalar_field_)[cell];
  for (int i=0; i<dim; i++) value += (*gradient_)[i][cell] * xc[i];

  return value;
}


/* ******************************************************************
 * calculates a value at point p using gradinet and centroid
****************************************************************** */
double Reconstruction::getValue(
    AmanziGeometry::Point& gradient, const int cell, const AmanziGeometry::Point& p)
{
  AmanziGeometry::Point xc(dim);
  xc = p;
  xc -= mesh_->cell_centroid(cell); 

  double value = (*scalar_field_)[cell] + (gradient * xc);

  return value;
}


/* ******************************************************************
 * Assemble a SPD least square matrix
****************************************************************** */
void Reconstruction::populateLeastSquareSystem(AmanziGeometry::Point& centroid, 
                                               double field_value, 
                                               Teuchos::SerialDenseMatrix<int, double>& matrix, 
                                               double* rhs)
{
  for (int i=0; i<dim; i++) {
    double xyz = centroid[i];

    matrix(i,i) += xyz * xyz; 
    for (int j=i+1; j<dim; j++) matrix(j,i) = matrix(i,j) += xyz * centroid[j];

    rhs[i] += xyz * field_value;
  }
}


/* ******************************************************************
 * Optimized linear algebra: norm, determinant, etc
****************************************************************** */
double Reconstruction::calculate_matrix_norm(Teuchos::SerialDenseMatrix<int, double>& matrix)
{
  double a = 0.0;
  for (int i=0; i<dim; i++) {
    for (int j=i; j<dim; j++) a = std::max(a, matrix(i,j));
  }
  return a;
}

double Reconstruction::calculate_matrix_determinant(Teuchos::SerialDenseMatrix<int, double>& matrix)
{
  double a = 0.0;
  if (dim == 2) {
    a = matrix(0,0) * matrix(1,1) - matrix(0,1) * matrix(0,1);    
  } else if (dim == 3) {
    a = matrix(0,0) * matrix(1,1) * matrix(2,2)
      + matrix(0,1) * matrix(1,2) * matrix(2,0) * 2
      - matrix(0,2) * matrix(1,1) * matrix(2,0)
      - matrix(0,1) * matrix(1,0) * matrix(2,2)
      - matrix(0,0) * matrix(1,2) * matrix(2,1);
  } else {
    a = matrix(0,0);
  }
  return a;
}

int Reconstruction::find_minimal_diagonal_entry(Teuchos::SerialDenseMatrix<int, double>& matrix)
{
  double a = matrix(0,0);  // We assume that matrix is SPD.
  int k = 0;

  for (int i=1; i<dim; i++) {
    double b = matrix(i,i);
    if (b < a) { 
      a = b;
      k = i;
    }
  }
  return k;
}


/* ******************************************************************
 * I/O routines
****************************************************************** */
void Reconstruction::printLeastSquareSystem(Teuchos::SerialDenseMatrix<int, double>matrix, double* rhs)
{
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) std::printf("%6.3f ", matrix(i,j));
    std::printf("  f[%1d] =%8.5f\n", i, rhs[i]);
  }
  std::printf("\n");
}


}  // namespace AmanziTransport
}  // namespace Amanzi

