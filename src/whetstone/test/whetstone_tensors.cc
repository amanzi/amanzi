/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "UnitTest++.h"

#include "Tensor.hh"


/* ****************************************************************
* Test of full tensor inversion
**************************************************************** */
TEST(TENSOR_INVERSE)
{
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Inversion of Tensors rank 2" << std::endl;

  Tensor T(2, 2);
  T(0, 0) = 1.1;
  T(0, 1) = -2.0;
  T(1, 0) = 0.5;
  T(1, 1) = 3.0;

  T.Inverse();

  CHECK_CLOSE(T(0, 0), 0.697674418604651, 1e-12);
  CHECK_CLOSE(T(0, 1), 0.465116279069767, 1e-12);
  CHECK_CLOSE(T(1, 0), -0.116279069767442, 1e-12);
  CHECK_CLOSE(T(1, 1), 0.255813953488372, 1e-12);
}


/* ****************************************************************
* Test of tensor initialization
**************************************************************** */
TEST(TENSOR_CONSTRUCTOR)
{
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Construction of Tensors rank 2" << std::endl;

  double data[4];

  Tensor T1;
  Tensor T2a(2, 1), T2b(2, 2), T2c(3, 1), T2d(3, 2), T2e(2, 4);

  Tensor T3(2, 2, data);

  Tensor T4(T2d);
}


/* ****************************************************************
* Rotation of tensors
**************************************************************** */
TEST(TENSOR_ROTATION_RANK2)
{
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Rotation of tensors of rank 2" << std::endl;

  double data1[4] = { 1.0, -1.0, -1.0, 2.0 };
  Tensor T1(2, 2, data1);

  std::vector<double> angles(3);
  angles[0] = M_PI / 2;

  T1.Rotate(angles.data());
  CHECK_CLOSE(T1(0, 0), 2.0, 1e-12);
  CHECK_CLOSE(T1(0, 1), 1.0, 1e-12);

  double data2[16] = {
   100, 20, 30, 30,
    20, 60, 10, 15,
    30, 10, 80, 25,
    30, 15, 25, 80
  };
  Tensor T2(2, 4, data2);
  angles[0] = M_PI / 3;
  T2.Rotate(angles.data());
  CHECK_CLOSE(T2(0, 0), 57.6466, 1e-3);
}

