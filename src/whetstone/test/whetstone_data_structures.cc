/*
  WhetStone

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

// TPLs
#include "UnitTest++.h"

// Amanzi::WhetStone
#include <DenseMatrix.hh>

#include <Kokkos_Vector.hpp>

struct test{}; 

/* ****************************************************************
* Test DenseMatrix 
 **************************************************************** */
TEST_FIXTURE(test, DENSEMATRIX)
{
  std::cout<<"Default execution space: "<<typeid(Kokkos::DefaultExecutionSpace).name()<<std::endl;
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;
  Kokkos::vector<DenseMatrix> m2;
  //Kokkos::vector<Mat> m;
  const int nb = 4;
  for(int i = 0 ; i < nb ; ++i){
      //m.push_back(Mat(i+2,i+3));
      //std::cout<<m[i]<<std::endl;
      m2.push_back(DenseMatrix(i+2,i+3)); 
      std::cout<<m2[i]<<std::endl; 
  }
  
  Kokkos::parallel_for(nb, KOKKOS_LAMBDA(const int& i){
      Kokkos::vector<DenseMatrix> lm = m2; 
      DenseMatrix mm = lm[i]; 
      mm.putScalar(i); 
      mm *= 2; 
      //mm.mult(2);
  });

  for(int i = 0 ; i < nb ; ++i){
    std::cout<<m2[i]<<std::endl;
  }
}
