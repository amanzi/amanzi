/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <cmath>
#include <iostream>

#include "UnitTest++.h"

#include "Brent.hh"

using namespace Amanzi;
using namespace Amanzi::Utils;

struct Tol {
  Tol(double eps) : eps_(eps){};
  bool operator()(double a, double b) const { return std::abs(a - b) <= eps_; }
  double eps_;
};

TEST(BRENT_AND_TOMS748)
{
  int itr;
  double x0, tol(1e-15);

  // original Boost iteration numbers
  Tol tol2(1e-15);
  std::vector<int> itr2({ 3, 9, 11, 9, 20, 22 });

  std::cout << " #        Root     Brent   TOMS748\n\n";

  // Test 1
  itr = 100;
  x0 = findRootBrent([](double x) { return x - 0.5; }, 0.0, 1.0, tol, &itr);
  CHECK_CLOSE(0.5, x0, tol);
  printf(" 1 %11.5f %9i %9i\n", x0, itr, itr2[0]);

  // Test 2
  itr = 100;
  x0 = findRootBrent([](double x) { return std::sin(x) - x / 2; }, 1.57, 3.14, tol, &itr);
  printf(" 2 %11.5f %9i %9i\n", x0, itr, itr2[1]);

  // Test 3
  itr = 100;
  x0 =
    findRootBrent([](double x) { return 2 * (x * std::exp(-20.0) - std::exp(-20 * x)) + 1.0; },
                     0.0,
                     0.9,
                     tol,
                     &itr);
  printf(" 3 %11.5f %9i %9i\n", x0, itr, itr2[2]);

  // Test 4
  itr = 100;
  x0 = findRootBrent(
    [](double x) { return x ? (x / 1.5 + std::sin(x) - 1) / 2 : -0.5; }, -1.0e+4, 1.57, tol, &itr);
  printf(" 4 %11.5f %9i %9i\n", x0, itr, itr2[3]);

  // Test 5
  itr = 100;
  x0 = findRootBrent([](double x) { return (20 * x - 1) / (19 * x); }, 0.01, 0.8, tol, &itr);
  printf(" 5 %11.5f %9i %9i\n", x0, itr, itr2[4]);

  // Test 6
  itr = 100;
  x0 = findRootBrent([](double x) { return std::pow(x, 12) - 0.2; }, 0.0, 5.0, tol, &itr);
  printf(" 6 %11.5f %9i %9i\n", x0, itr, itr2[5]);
}
