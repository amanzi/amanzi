#ifndef __cell_geometry__
#define __cell_geometry__

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
  // Area of the quadrilateral face {x1, x2, x3, x4} in R^3.
  double quad_face_area(double x1[], double x2[], double x3[], double x4[]);
  // Volume of the tetrahedron {x1, x2, x3, x4}.
  double tet_volume(double x1[], double x2[], double x3[], double x4[]);
}

#endif
