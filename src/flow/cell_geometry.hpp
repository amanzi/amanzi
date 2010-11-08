#ifndef __cell_geometry__
#define __cell_geometry__

#include "Epetra_SerialDenseMatrix.h"

namespace cell_geometry {
  // Dot product of n-vectors a and b.
  double dot_product(double a[], double b[], int n);
  // Length of the n-vector x.
  double vector_length(double x[], int n);
  // Cross product of two 3-vectors.
  void cross_product(double result[], double a[], double b[]);
  // Triple product of three 3-vectors.
  double triple_product(double a[], double  b[], double c[]);
  // Area-weighted normal (oriented area) of the quadrilateral face {x1, x2, x3, x4} in R^3.
  void quad_face_normal(double result[], double x1[], double x2[], double x3[], double x4[]);
  
  // Area-weighted normal (oriented area) of the quadrilateral face {x[0], x[1], x[2], x[3]} in R^3.
  void quad_face_normal(double result[], double x[4][3]);
  
  void quad_face_normal(double result[], Epetra_SerialDenseMatrix &x);

    // Area of the quadrilateral face {x1, x2, x3, x4} in R^3.
  double quad_face_area(double x1[], double x2[], double x3[], double x4[]);
  // Volume of the tetrahedron {x1, x2, x3, x4}.
  double tet_volume(double x1[], double x2[], double x3[], double x4[]);
  // Volume of the hexahedron {x[0], x[1], ..., x[7]}.
  double hex_volume(Epetra_SerialDenseMatrix &x);
  // Volumes of the hexahedron {x[0], x[1], ..., x[7]} and eight corner tetrahedra volumes.
  void compute_hex_volumes(Epetra_SerialDenseMatrix &x, double &hvol, double cvol[]);
  // Area-weighted outward normals (oriented area) of the quadrilateral faces of hexahedron.
  void compute_hex_face_normals(Epetra_SerialDenseMatrix &x, Epetra_SerialDenseMatrix &a);
}

#endif
