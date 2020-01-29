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
#include <DenseVector.hh>
#include <Tensor.hh>

#include <Kokkos_Vector.hpp>

struct test{}; 

/* ****************************************************************
* Test DenseMatrix 
 **************************************************************** */
// TEST_FIXTURE(test, DENSEMATRIX)
// {
//   std::cout<<"Default execution space: "<<
//     typeid(Kokkos::DefaultExecutionSpace).name()<<std::endl;
//   using namespace Amanzi;
//   using namespace Amanzi::WhetStone;
//   const int nb = 4;
//   Kokkos::vector<DenseMatrix> m1;
//   Kokkos::vector<DenseMatrix> m2;
//   Kokkos::View<double*> det("Determinant",nb); 

//   for(int i = 0 ; i < nb ; ++i){
//       m1.push_back(DenseMatrix(i+2,i+2));
//       std::cout<<m1[i]<<std::endl;
//       m2.push_back(DenseMatrix(i+2,i+2)); 
//       std::cout<<m2[i]<<std::endl; 
//   }

//   Kokkos::parallel_for(nb, KOKKOS_LAMBDA(const int& i){
//       DenseMatrix mm1 = m1[i];
//       DenseMatrix mm2 = m2[i];
//       mm1.putScalar(i+1); 
//       mm1 *= 2;
//       mm2.putScalar(i+1);
//       mm1 += mm2;  
//       assert(mm1.Multiply(mm1,mm2,true) == 0); 
//       det(i) = mm1.Det(); 
//   });

//   for(int i = 0 ; i < nb ; ++i){
//     std::cout<<m1[i]<<std::endl;
//     std::cout<<"Det = "<<det(i)<<std::endl;
//     std::cout<<m2[i]<<std::endl;
//   }
// }

// /* ****************************************************************
// * Test DenseVector 
//  **************************************************************** */
// TEST_FIXTURE(test, DENSEVECTOR)
// {
//   std::cout<<"Default execution space: "<<
//     typeid(Kokkos::DefaultExecutionSpace).name()<<std::endl;
//   using namespace Amanzi;
//   using namespace Amanzi::WhetStone;
//   const int nb = 4;
//   Kokkos::vector<DenseVector> m1;

//   for(int i = 0 ; i < nb ; ++i){
//       m1.push_back(DenseVector(i+2));
//       std::cout<<m1[i]<<std::endl;
//   }

//   Kokkos::parallel_for(nb, KOKKOS_LAMBDA(const int& i){
//       m1[i].putScalar(i+1); 
//   });

//   for(int i = 0 ; i < nb ; ++i){
//     std::cout << m1[i] << std::endl;
//   }
// }

/* ****************************************************************
* Test Tensor 
 **************************************************************** */
TEST_FIXTURE(test, TENSOR)
{
  std::cout<<"Default execution space: "<<
    typeid(Kokkos::DefaultExecutionSpace).name()<<std::endl;
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;
  std::cout << "construct vector" << std::endl;
  Kokkos::vector<Tensor> m1;
  m1.reserve(20);

  int i = 0;
  std::cout << "Construct 1:" << std::endl;
  for(int d = 0; d < 4 ; ++d){
    for(int r = 0; r < 5 ; ++r){
      std::cout << " push back" << std::endl;
      Tensor mm1(d,r);
      m1.push_back(mm1);
      i++;
    }
  }

  std::cout << "Init on device 1:" << std::endl;
  Kokkos::parallel_for(
      m1.size(),
      KOKKOS_LAMBDA(const int& i){
        Tensor& mm1 = m1[i];
        mm1.putScalar(i+1); 
      });

  std::cout << "Print on host:" << std::endl;
  for(int i = 0 ; i < m1.size() ; ++i){
    std::cout << m1[i] << std::endl;
  }
}

