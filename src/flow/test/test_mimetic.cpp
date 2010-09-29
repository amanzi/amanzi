#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../upper_packed_matrix.hpp"
#include "../mimetic_hex.hpp"


TEST(CROSS_PRODUCT) {

  

  // unit vectors
  double e[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  // Check result for all combinations of the standard unit vectors; bit-exact.
  for (int k=0; k<3; k++) {
    
    int kp1 = k % 3;
    int kp2 = (k+1) % 3;

    double aux[3];
    double zero_array[3] = { 0.0, 0.0, 0.0 };

    cross_product(aux, e[k], e[kp1]);
    CHECK_ARRAY_EQUAL(aux, zero_array, 3);

    cross_product(aux, e[kp1], e[k]);
    CHECK_ARRAY_EQUAL(aux, zero_array, 3);
    
    cross_product(aux, e[k], e[k]);
    CHECK_ARRAY_EQUAL(aux, zero_array, 3);
  }
    


  // Check result for random vectors against the bilinear reconstruction.
  // Floating point values have at most 8 significant bits so the result is bit-exact.
  
  srand((unsigned)(time(0))); 

  for (int n=0; n<10; n++) {
    
    double a[3], b[3];
    
    for (int i=0; i<3; i++) {
      a[i] = rand()/(double(RAND_MAX)+1); 
      a[i] = int(512*(a[i] - 0.5));

      b[i] = rand()/(double(RAND_MAX)+1);
      b[i] = int(512*(b[i] - 0.5));
    }
      
    double c[3] = { 0.0, 0.0, 0.0 };
    double aux[3];

    for (int j = 0; j<3; j++) {
      for (int k = 0; k<3; k++) {
	cross_product(aux, e[j], e[k]);
	for (int i=0; i<3; i++) 
	  c[i] += a[j]*b[k]* aux[i];
      }
    }
    
    cross_product(aux, a, b);
    
    CHECK_ARRAY_EQUAL(c, aux, 3);
  }


}


TEST(TRIPLE_PRODUCT) {

  srand((unsigned)(time(0))); 

  for (int n=0; n<10; n++) {

    double a[3], b[3], c[3];

    for (int i=0; i<3; i++) {
      a[i] = rand()/(double(RAND_MAX)+1); 
      b[i] = rand()/(double(RAND_MAX)+1); 
      c[i] = rand()/(double(RAND_MAX)+1); 

      a[i] = int(512*(b[i] - 0.5));
      b[i] = int(512*(b[i] - 0.5));
      c[i] = int(512*(b[i] - 0.5));     
    }
    
    double tp = triple_product(a, b, c);

    double aux[3];
    cross_product(aux, b, c);
    
    double dpcp = dot_product(a, aux);

    CHECK_EQUAL(tp, dpcp);

  }

}


TEST(TET_VOLUME) {

  // 125 * orthonormal rotation matrix (based on 3^2 + 4^2 = 5^2).
  double q[3][3] = { { 45, -12, 116 }, { 60, 109, -12 }, { -100, 60, 45 } };

  // Reference tet; origin corner of 1x2x3 orthogonal, axis-aligned hex.
  double ref_tet[4][3] = { { 0,0,0 }, { 1,0,0 }, { 0,2,0 }, { 0,0,3 } };

  CHECK_EQUAL(tet_volume(ref_tet[0], ref_tet[1], ref_tet[2], ref_tet[3]), 1.0);
  

  double a[3];
  
  srand((unsigned)(time(0))); 
  for (int i=0; i<3; i++) {
    a[i] =  rand()/(double(RAND_MAX)+1);
    a[i] =  int(512*(a[i] - 0.5));
  }
  
  double shifted[4][3];
  for (int j=0; j<4; j++)
    for (int i=0; i<3; i++)  
      shifted[j][i] = ref_tet[j][i] - a[i];

  // make sure that translation does not change the tet volume
  CHECK_EQUAL(tet_volume(shifted[0], shifted[1], shifted[2], shifted[3]), 1.0);
  
  double aux[4][3];
  matmul(aux, shifted, q, 4);
  
  CHECK_EQUAL(tet_volume(aux[0], aux[1], aux[2], aux[3]), pow(125,3));

  double s[3][3];
  for (int j=0; j<3; j++)
    for (int k=0; k<3; k++)
      s[k][j] = q[1][j] * q[2][k];
  
  for (int j=0; j<3; j++)
    s[j][j] += 1.0;

  matmul(aux, shifted, s, 4);
  
  CHECK_EQUAL(tet_volume(aux[0], aux[1], aux[2], aux[3]), 1.0);
  
}







TEST(MIMETIC) {

  double x1[8][3] = { { 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 },
		      { 0,0,1 }, { 1,0,1 }, { 1,1,1 }, { 0,1,1 } };


  // 125 * orthonormal rotation matrix (based on 3^2 + 4^2 = 5^2).
  double q[3][3] = { { 45, -12, 116 }, { 60, 109, -12 }, { -100, 60, 45 } };

  // remove the scaling of q
  for (int j=0; j<3; j++) 
    for (int i=0; i<3; i++)  q[j][i] /= 125;

  
  mimetic_hex cell1 ( x1 );
  
  double m1[21];
  
  cell1.mass_matrix(m1);

  
  double disp[3];
  for (int i=0; i<3; i++) {
    disp[i] = rand()/(double(RAND_MAX)+1); 
    disp[i] = int(512*(disp[i] - 0.5));
  }

  double aux[8][3];
  
  for (int j=0; j<8; j++)
    for (int i=0; i<3; i++) 
      aux[j][i] = x1[j][i] - disp[i];

  
  double x2[8][3];
  matmul(x2, aux, q, 8);
  
  mimetic_hex cell2( x2 );
  
  double m2[21];
  
  cell2.mass_matrix( m2 );

  CHECK_ARRAY_CLOSE( m1, m2, 21, 0.000001 );
      

}
