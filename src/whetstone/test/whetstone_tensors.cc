/*
The transport component of the Amanzi code, serial unit tests.
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "UnitTest++.h"

#include "tensor.hh"


/* ****************************************************************
* TBW
**************************************************************** */
TEST(TENSOR_INVERSE) {
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Inversion of Tensors rank 2" << std::endl;
 
  Tensor T(2, 2);
  T(0, 0) = 1.1;
  T(0, 1) =-2.0;
  T(1, 0) = 0.5;
  T(1, 1) = 3.0;

  T.inverse();

  CHECK_CLOSE(T(0, 0), 0.697674418604651, 1e-12);
  CHECK_CLOSE(T(0, 1), 0.465116279069767, 1e-12);
  CHECK_CLOSE(T(1, 0),-0.116279069767442, 1e-12);
  CHECK_CLOSE(T(1, 1), 0.255813953488372, 1e-12);
}


