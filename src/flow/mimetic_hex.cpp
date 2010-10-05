#include <iostream>
#include "mimetic_hex.hpp"
#include "upper_packed_matrix.hpp"
#include <math.h>

mimetic_hex::mimetic_hex(double x_[8][3])
{
  for (int j=0; j<8; j++)
    for (int i=0; i<3; i++)
      x[j][i] = x_[j][i];
  
  compute_hex_volumes();
  
  compute_hex_face_normals();
}



void mimetic_hex::mass_matrix(double* matrix, double scalar_coeff)
{
  // the scalar coefficient is assumed to be constant
  // over the cell


  int  cfaceloc[8][3] = { { 1, 2, 4 }, { 1, 4, 3 }, { 0, 3, 4 }, { 0, 4, 2 }, 
				{ 1, 5, 2 }, { 1, 3, 5 }, { 0, 5, 3 }, { 0, 2, 5 } };  
   
  // set the mass matrix to zero
  // matrix is an upper packed matrix, of size 6x6
  // hence, its dimension is 21
  for (int i=0; i<21; i++)  matrix[i] = 0.0;

  // compute the sum of all corner tet volumes in the cell
  double sumvolnode = 0.0;
  for (int i=0; i<8; i++) sumvolnode += cvol[i];

  
  for (int corner=0; corner<8; corner++) {

    // for the current corner initialize Jac with the 
    // there normal vectors of its three faces
    double Jac[3][3] = { 0.0 };
    
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++) 
	Jac[i][j] = face_normal[ cfaceloc[corner][j] ] [i];
    
    double EMM[3][3] = { 0.0 };
    double aux1[3][3] = { 0.0 };
    double aux2[3][3] = { 0.0 };

    transpose_3x3( aux1, Jac );
    matmul( aux2, Jac, aux1, 3 ); 
      
    invert_sym_3x3( EMM, aux2 );
    
    // scaling factor for the corner matrix (conductivity divided by corner volume)
    double scale = scalar_coeff * cvol[corner]*hvol/sumvolnode;
    
    // loop over the faces adjacent to corner cnr
    for (int jface = 0; jface<3; jface++) {
      
      // retrieve the index of the face in terms of the cell
      int jcfl = cfaceloc[jface][corner];
      
      // loop over the faces adjacent to the corner cnr (only do the upper triangular part)
      for (int iface = 0; iface<3; iface++) {
	
	// retrieve the index of this face in term of the cell
	int icfl = cfaceloc[corner][iface];
	
	// compute the index in the cell mass matrix
	// packed upper storage
	if (icfl <= jcfl) {
	  int ijcfl = icfl + jcfl*(jcfl-1)/2;
	  matrix[ijcfl] += scale*EMM[iface][jface];
	}
      }
    }
  }
  
}

void mimetic_hex::mass_matrix(double* matrix, double tensor_coeff[][3])
{
  // the tensor coefficient is assumed to be a 3x3 matrix
  // that is constant over the cell.


  int  cfaceloc[8][3] = { { 1, 2, 4 }, { 1, 4, 3 }, { 0, 3, 4 }, { 0, 4, 2 }, 
				{ 1, 5, 2 }, { 1, 3, 5 }, { 0, 5, 3 }, { 0, 2, 5 } };  
   
  // set the mass matrix to zero
  // matrix is an upper packed matrix, of size 6x6
  // hence, its dimension is 21
  for (int i=0; i<21; i++)  matrix[i] = 0.0;

  // compute the sum of all corner tet volumes in the cell
  double sumvolnode = 0.0;
  for (int i=0; i<8; i++) sumvolnode += cvol[i];

  // compute the transpose of the tensor coefficient, internal to this
  // mimetic_hex class we use transposes of matrices, since we want 
  // to use contiguous storage for vectors 
  double trans_tensor_coeff[3][3];
  transpose_3x3(trans_tensor_coeff, tensor_coeff);

  for (int corner=0; corner<8; corner++) {

    // for the current corner initialize Jac with the 
    // there normal vectors of its three faces
    double Jac[3][3] = { 0.0 };
    
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++) 
	Jac[i][j] = face_normal[ cfaceloc[corner][j] ] [i];
    
    double EMM[3][3] = { 0.0 };
    double trans_Jac[3][3] = { 0.0 };
    double aux1[3][3] = { 0.0 };
    double aux2[3][3] = { 0.0 };

    transpose_3x3( trans_Jac, Jac );
    matmul( aux1, trans_Jac, trans_tensor_coeff, 3 ); 

    matmul( aux2, Jac, aux1, 3 ); 
      
    invert_sym_3x3( EMM, aux2 );
    
    // scaling factor for the corner matrix (conductivity divided by corner volume)
    double scale = cvol[corner]*hvol/sumvolnode;
    
    // loop over the faces adjacent to corner cnr
    for (int jface = 0; jface<3; jface++) {
      
      // retrieve the index of the face in terms of the cell
      int jcfl = cfaceloc[jface][corner];
      
      // loop over the faces adjacent to the corner cnr (only do the upper triangular part)
      for (int iface = 0; iface<3; iface++) {
	
	// retrieve the index of this face in term of the cell
	int icfl = cfaceloc[corner][iface];
	
	// compute the index in the cell mass matrix
	// packed upper storage
	if (icfl <= jcfl) {
	  int ijcfl = icfl + jcfl*(jcfl-1)/2;
	  matrix[ijcfl] += scale*EMM[iface][jface];
	}
      }
    }
  }
  
}

void mimetic_hex::mass_matrix(double* matrix)
{

  int  cfaceloc[8][3] = { { 1, 2, 4 }, { 1, 4, 3 }, { 0, 3, 4 }, { 0, 4, 2 }, 
				{ 1, 5, 2 }, { 1, 3, 5 }, { 0, 5, 3 }, { 0, 2, 5 } };  
   
  // set the mass matrix to zero
  // matrix is an upper packed matrix, of size 6x6
  // hence, its dimension is 21
  for (int i=0; i<21; i++)  matrix[i] = 0.0;

  // compute the sum of all corner tet volumes in the cell
  double sumvolnode = 0.0;
  for (int i=0; i<8; i++) sumvolnode += cvol[i];

  
  for (int corner=0; corner<8; corner++) {

    // for the current corner initialize Jac with the 
    // there normal vectors of its three faces
    double Jac[3][3] = { 0.0 };
    
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++) 
	Jac[i][j] = face_normal[ cfaceloc[corner][j] ] [i];
    
    double EMM[3][3] = { 0.0 };
    double aux1[3][3] = { 0.0 };
    double aux2[3][3] = { 0.0 };

    transpose_3x3( aux1, Jac );
    matmul( aux2, Jac, aux1, 3 ); 
      
    invert_sym_3x3( EMM, aux2 );
    
    // scaling factor for the corner matrix
    // without a coefficient (scalar or tensor) this scaling
    // is purely geometric: 1/(corner volume) 
    double scale = cvol[corner]*hvol/sumvolnode;
    
    // loop over the faces adjacent to corner cnr
    for (int jface = 0; jface<3; jface++) {
      
      // retrieve the index of the face in terms of the cell
      int jcfl = cfaceloc[jface][corner];
      
      // loop over the faces adjacent to the corner cnr (only do the upper triangular part)
      for (int iface = 0; iface<3; iface++) {
	
	// retrieve the index of this face in term of the cell
	int icfl = cfaceloc[corner][iface];
	
	// compute the index in the cell mass matrix
	// packed upper storage
	if (icfl <= jcfl) {
	  int ijcfl = icfl + jcfl*(jcfl-1)/2;
	  matrix[ijcfl] += scale*EMM[iface][jface];
	}
      }
    }
  }
  
}


void mimetic_hex::diff_op(double coef,
    const double &pcell, const double pface[], double &rcell, double rface[])
{
  double matrix[21], aux[6];
  
  mass_matrix(matrix, 1.0/coef);
  upper_packed_matrix minv(matrix, 6);
  minv.invert();
  
  for (int i = 0; i < 6; ++i) aux[i] = pface[i] - pcell;
    
  minv.sym_matmul(aux, rface);
  
  rcell = 0.0;
  for (int i = 0; i < 6; ++i) rcell -= rface[i];
}
  

void  transpose_3x3(double matout[][3], double matin[][3]) 
{
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) 
      matout[i][j] = matin[j][i];
}


void matmul(double matout[][3], double mat1[][3], double mat2[3][3], int m) 
{
  // m is the leading dimension of mat2 and matout

  for (int rows=0; rows<m; rows++)

    for (int j=0; j<3; j++) {
      matout[rows][j] = 0.0;
      for (int k=0; k<3; k++)
	matout[rows][j] += mat1[rows][k] * mat2[k][j];
    }
}


void invert_sym_3x3(double minv[][3], double mat[][3])
{

   double invdet =  1.0 /
     (mat[0][0]*(mat[1][1]*mat[2][2] - pow(mat[2][1],2)) 
      - mat[0][1]*(mat[0][1]*mat[2][2] - 2.0*mat[0][2]*mat[1][2]) 
      - mat[1][1]*pow(mat[0][2],2));


   minv[0][0] =  (mat[1][1]*mat[2][2]-pow(mat[2][1],2))*invdet;
   minv[0][1] = -(mat[1][0]*mat[2][2]-mat[2][0]*mat[1][2])*invdet;
   minv[0][2] =  (mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1])*invdet;
   minv[1][0] = minv[0][1];
   minv[1][1] =  (mat[0][0]*mat[2][2]-pow(mat[2][0],2))*invdet;
   minv[1][2] = -(mat[0][0]*mat[2][1]-mat[2][0]*mat[0][1])*invdet;
   minv[2][0] = minv[0][2];
   minv[2][1] = minv[1][2];
   minv[2][2] =  (mat[0][0]*mat[1][1]-pow(mat[1][0],2))*invdet;
   
}


double dot_product(double a[], double b[])
{
  double dp = a[0]*b[0];
  for (int i=1; i<3; i++)
    dp += a[i]*b[i];
 
  return dp;
}


void cross_product(double result[], double a[], double b[])
{
  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];
};



double triple_product(double a[], double b[], double c[])
{
  return a[0]*(b[1]*c[2] - b[2]*c[1]) +
    a[1]*(b[2]*c[0] - b[0]*c[2]) + 
    a[2]*(b[0]*c[1] - b[1]*c[0]);
};



double tet_volume(double x1[], double x2[], double x3[], double x4[] )
{
  double v1[3], v2[3], v3[3];
  
  for (int i=0; i<3; i++) {
    v1[i] = x2[i]-x1[i];
    v2[i] = x3[i]-x1[i];
    v3[i] = x4[i]-x1[i];
  }

  return triple_product(v1, v2, v3) / 6.0;

};

void mimetic_hex::compute_hex_volumes()
{

  cvol[0] = tet_volume(x[0], x[1], x[3], x[4]);
  cvol[1] = tet_volume(x[1], x[2], x[0], x[5]); 
  cvol[2] = tet_volume(x[2], x[3], x[1], x[6]); 
  cvol[3] = tet_volume(x[3], x[0], x[2], x[7]); 
  cvol[4] = tet_volume(x[4], x[7], x[5], x[0]); 
  cvol[5] = tet_volume(x[5], x[4], x[6], x[1]);
  cvol[6] = tet_volume(x[6], x[5], x[7], x[2]);
  cvol[7] = tet_volume(x[7], x[6], x[4], x[3]); 

  
  hvol = 0.0;
  for (int i=0; i<8; i++) hvol += cvol[i];

  hvol += tet_volume(x[0],x[2],x[7],x[5]) + tet_volume(x[1],x[3],x[4],x[6]);
  
  hvol *= 0.5;
}


void mimetic_hex::compute_hex_face_normals()
{
  double v1[3], v2[3];
  
  for (int i=0; i<3; i++) {
    v1[i] = x[7][i] - x[2][i];
    v2[i] = x[6][i] - x[3][i];
  }
  cross_product(face_normal[0], v1, v2);

  for (int i=0; i<3; i++) {
    v1[i] = x[5][i] - x[0][i];
    v2[i] = x[4][i] - x[1][i];
  }
  cross_product(face_normal[1], v1, v2);  

  for (int i=0; i<3; i++) {
    v1[i] = x[4][i] - x[3][i];
    v2[i] = x[7][i] - x[0][i];
  }
  cross_product(face_normal[2], v1, v2);  
  
  for (int i=0; i<3; i++) {
    v1[i] = x[6][i] - x[1][i];
    v2[i] = x[5][i] - x[2][i];
  }
  cross_product(face_normal[3], v1, v2);  

  for (int i=0; i<3; i++) {
    v1[i] = x[2][i] - x[0][i];
    v2[i] = x[1][i] - x[3][i];
  }
  cross_product(face_normal[4], v1, v2);  

  for (int i=0; i<3; i++) {
    v1[i] = x[6][i] - x[4][i];
    v2[i] = x[7][i] - x[5][i];
  }
  cross_product(face_normal[5], v1, v2);  
  
  for (int j=0; j<6; j++)
    for (int i=0; i<3; i++)
      face_normal[j][i] *= 0.5;
  
}
