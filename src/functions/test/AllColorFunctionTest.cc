/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <vector>

#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include "VerboseObject_objs.hh"

#include "FunctionGridColor.hh"

using namespace Amanzi;

int
main(int argc, char* argv[])
{
  return UnitTest::RunAllTests();
}

SUITE(Grid)
{
  TEST(Grid1Da)
  {
    std::vector<int> count(1, 1);
    std::vector<double> x0(1, 0.0);
    std::vector<double> dx(1, 1.0);
    std::vector<int> array(1, 2);
    FunctionGridColor f(1, count, x0, dx, array);
    double x;
    x = -1.0;
    CHECK_EQUAL(f(&x), 2);
    x = 0.5;
    CHECK_EQUAL(f(&x), 2);
    x = 2.0;
    CHECK_EQUAL(f(&x), 2);
  }
  TEST(Grid1Db)
  {
    std::vector<int> count(1, 3);
    std::vector<double> x0(1, 1.0);
    std::vector<double> dx(1, 1.0);
    std::vector<int> array(3);
    array[0] = 1;
    array[1] = 2;
    array[2] = 3;
    FunctionGridColor f(1, count, x0, dx, array);
    double x;
    x = 0.5;
    CHECK_EQUAL(f(&x), 1);
    x = 1.5;
    CHECK_EQUAL(f(&x), 1);
    x = 2.5;
    CHECK_EQUAL(f(&x), 2);
    x = 3.5;
    CHECK_EQUAL(f(&x), 3);
    x = 4.5;
    CHECK_EQUAL(f(&x), 3);
  }
  TEST(Grid2Da)
  {
    // 1x1, [0,1]x[2,3]
    std::vector<int> count(2, 1);
    std::vector<double> x0(2);
    x0[0] = 1.0;
    x0[1] = 2.0;
    std::vector<double> dx(2);
    dx[0] = -1.0;
    dx[1] = 1.0;
    std::vector<int> array(1, 2);
    FunctionGridColor f(2, count, x0, dx, array);
    double x1[2] = { -1.0, 0.0 };
    CHECK_EQUAL(f(x1), 2);
    double x2[2] = { 2.0, 4.0 };
    CHECK_EQUAL(f(x2), 2);
    double x3[2] = { 0.5, 2.5 };
    CHECK_EQUAL(f(x3), 2);
    double x4[2] = { -1.0, 2.5 };
    CHECK_EQUAL(f(x4), 2);
  }
  TEST(Grid2Db)
  {
    // 2x2, [0,2]x[0,2]
    std::vector<int> count(2);
    count[0] = 2;
    count[1] = 2;
    std::vector<double> x0(2);
    x0[0] = 2.0;
    x0[1] = 0.0;
    std::vector<double> dx(2);
    dx[0] = -1.0;
    dx[1] = 1.0;
    std::vector<int> array(4);
    for (int i = 0; i < 4; ++i) array[i] = 1 + i;
    FunctionGridColor f(2, count, x0, dx, array);
    double x1[2] = { 1.5, 0.5 };
    CHECK_EQUAL(1, f(x1));
    double x2[2] = { 0.5, 0.5 };
    CHECK_EQUAL(2, f(x2));
    double x3[2] = { 1.5, 1.5 };
    CHECK_EQUAL(3, f(x3));
    double x4[2] = { 0.5, 1.5 };
    CHECK_EQUAL(4, f(x4));
    double x5[2] = { -0.5, 1.5 };
    CHECK_EQUAL(4, f(x5));
    double x6[2] = { 2.5, -0.5 };
    CHECK_EQUAL(1, f(x6));
  }
  TEST(Grid3Db)
  {
    // 2x3x4, [-0.5,1.5]x[-0.5,2.5]x[-0.5,3.5]
    std::vector<int> count(3);
    count[0] = 2;
    count[1] = 3;
    count[2] = 4;
    std::vector<double> x0(3, -0.5);
    std::vector<double> dx(3, 1.0);
    std::vector<int> array(24);
    int n = 0;
    for (int iz = 0; iz < count[2]; ++iz) {
      for (int iy = 0; iy < count[1]; ++iy) {
        for (int ix = 0; ix < count[0]; ++ix) {
          array[n] = 2 * ix - iy + iz;
          ++n;
        }
      }
    }
    FunctionGridColor f(3, count, x0, dx, array);
    int fx;
    double x[3];
    for (int iz = 0; iz < count[2]; ++iz) {
      for (int iy = 0; iy < count[1]; ++iy) {
        for (int ix = 0; ix < count[0]; ++ix) {
          fx = 2 * ix - iy + iz;
          x[0] = ix;
          x[1] = iy;
          x[2] = iz;
          CHECK_EQUAL(fx, f(x));
        }
      }
    }
    double x1[3] = { 2.0, 3.0, 4.0 };
    CHECK_EQUAL(3, f(x1));
  }
}
