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
  const int nb = 4;
  Kokkos::vector<DenseMatrix> m1;
  Kokkos::vector<DenseMatrix> m2;
  Kokkos::View<double*> det("Determinant",nb); 

  for(int i = 0 ; i < nb ; ++i){
      m1.push_back(DenseMatrix(i+2,i+2));
      std::cout<<m1[i]<<std::endl;
      m2.push_back(DenseMatrix(i+2,i+2)); 
      std::cout<<m2[i]<<std::endl; 
  }

  Kokkos::parallel_for(nb, KOKKOS_LAMBDA(const int& i){
      Kokkos::vector<DenseMatrix> lm1 = m1; 
      Kokkos::vector<DenseMatrix> lm2 = m2; 
      DenseMatrix mm1 = lm1[i]; 
      DenseMatrix mm2 = lm2[i];
      mm1.putScalar(i+1); 
      mm1 *= 2;
      mm2.putScalar(i+1);
      mm1 += mm2;  
      assert(mm1.Multiply(mm1,mm2,true) == 0); 
      det(i) = mm1.Det(); 
  });

  for(int i = 0 ; i < nb ; ++i){
    std::cout<<m1[i]<<std::endl;
    std::cout<<"Det = "<<det(i)<<std::endl;
    std::cout<<m2[i]<<std::endl;
  }
}
