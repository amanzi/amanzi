#ifndef __mimetic_cell_hpp__
#define __mimetic_cell_hpp__


void invert_sym_3x3 (double minv[][3], double mat[][3] );
void transpose_3x3 (double matout[][3], double matin[][3] );
void matmul (double matout[][3], double mat1[][3], double mat2[][3], int m );

double dot_product(double a[], double b[]);
void cross_product(double result[], double a[], double b[]);
double triple_product(double a[], double b[], double c[]);
double tet_volume(double x1[], double x2[], double x3[], double x4[]);




class mimetic_hex {

public:
  mimetic_hex(double x_[][3]);
  ~mimetic_hex() {};

  void mass_matrix(double* matrix);
  void mass_matrix(double* matrix, double);
  void mass_matrix(double* matrix, double[][3]);
  
  void diff_op(double, const double&, const double[], double&, double[]);

private:

  
  void compute_hex_volumes();
  void compute_hex_face_normals();

private:
  double x[8][3];
  double cvol[8];
  double hvol;
  double face_normal[6][3];
};

#endif
