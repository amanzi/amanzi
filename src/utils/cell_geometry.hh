#ifndef __cell_geometry__
#define __cell_geometry__

#include "Epetra_SerialDenseMatrix.h"

namespace cell_geometry {
  // Dot product of n-vectors a and b.
  double dot_product(const double a[], const double b[], int n);
  // Length of the n-vector x.
  double vector_length(const double x[], int n);
  // Cross product of two 3-vectors.
  void cross_product(const double a[], const double b[], double result[]);
  // Triple product of three 3-vectors.
  double triple_product(const double a[], const double  b[], const double c[]);
  // Area-weighted normal (oriented area) of the quadrilateral face {x1, x2, x3, x4} in R^3.
  void quad_face_normal(const double x1[], const double x2[], const double x3[], const double x4[], double result[]);
  
  // Area-weighted normal (oriented area) of the quadrilateral face {x[0], x[1], x[2], x[3]} in R^3.
  void quad_face_normal(const double x[4][3], double result[]);
  void quad_face_normal(const Epetra_SerialDenseMatrix &x, double result[]);
  void quad_face_normal(double *x, double result[]);

    // Area of the quadrilateral face {x1, x2, x3, x4} in R^3.
  double quad_face_area(const double x1[], const double x2[], const double x3[], const double x4[]);
  // Volume of the tetrahedron {x1, x2, x3, x4}.
  double tet_volume(const double x1[], const double x2[], const double x3[], const double x4[]);
  // Volume of the hexahedron {x[0], x[1], ..., x[7]}.
  double hex_volume(const Epetra_SerialDenseMatrix &x);
  // Volumes of the hexahedron {x[0], x[1], ..., x[7]} and eight corner tetrahedra volumes.
  void compute_hex_volumes(const Epetra_SerialDenseMatrix &x, double &hvol, double cvol[]);
  // Area-weighted outward normals (oriented area) of the quadrilateral faces of hexahedron.
  void compute_hex_face_normals(const Epetra_SerialDenseMatrix &x, Epetra_SerialDenseMatrix &a);
  
  // Centroid of the hexahedron {x[0], x[1], ..., x[y]}.  Correct for planar faces.
  void hex_centroid (const Epetra_SerialDenseMatrix &x, double c[]);
  
  // Centroid of the planar quadrilateral face {x[0], ..., x[3]}.
  void quad_face_centroid(const Epetra_SerialDenseMatrix &x, double c[]);
  
  // Area of the triangular face {x0, x1, x2} in R^3.
  double tri_face_area(const double x0[], const double x1[], const double x2[]);
  
  // Centroid of the triangular face {x0, x1, x2} in R^3.
  void tri_face_centroid(const double x0[], const double x1[], const double x2[], double c[]);
}

#endif
