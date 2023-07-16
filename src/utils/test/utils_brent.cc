/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>

#include <boost/math/tools/roots.hpp>
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

  boost::uintmax_t itr2;
  Tol tol2(1e-15);
  std::pair<double, double> xx;

  std::cout << " #        Root     Brent   TOMS748\n\n";

  // Test 1
  itr = 100;
  x0 = brent([](double x) { return x - 0.5; }, 0.0, 1.0, tol, &itr);
  CHECK_CLOSE(0.5, x0, tol);

  itr2 = 100;
  xx = boost::math::tools::toms748_solve([](double x) { return x - 0.5; }, 0.0, 1.0, tol2, itr2);
  printf(" 1 %11.5f %9i %9i\n", x0, itr, (int)itr2);

  // Test 2
  itr = 100;
  x0 = brent([](double x) { return std::sin(x) - x / 2; }, 1.57, 3.14, tol, &itr);

  itr2 = 100;
  xx = boost::math::tools::toms748_solve(
    [](double x) { return std::sin(x) - x / 2; }, 1.57, 3.14, tol2, itr2);
  printf(" 2 %11.5f %9i %9i\n", x0, itr, (int)itr2);

  // Test 3
  itr = 100;
  x0 = brent([](double x) { return 2 * (x * std::exp(-20.0) - std::exp(-20 * x)) + 1.0; },
             0.0,
             0.9,
             tol,
             &itr);

  itr2 = 100;
  xx = boost::math::tools::toms748_solve(
    [](double x) { return 2 * (x * std::exp(-20.0) - std::exp(-20 * x)) + 1.0; },
    0.0,
    0.9,
    tol2,
    itr2);
  printf(" 3 %11.5f %9i %9i\n", x0, itr, (int)itr2);

  // Test 4
  itr = 100;
  x0 = brent(
    [](double x) { return x ? (x / 1.5 + std::sin(x) - 1) / 2 : -0.5; }, -1.0e+4, 1.57, tol, &itr);

  itr2 = 100;
  xx = boost::math::tools::toms748_solve(
    [](double x) { return x ? (x / 1.5 + std::sin(x) - 1) / 2 : -0.5; }, -1.0e+4, 1.57, tol2, itr2);
  printf(" 4 %11.5f %9i %9i\n", x0, itr, (int)itr2);

  // Test 5
  itr = 100;
  x0 = brent([](double x) { return (20 * x - 1) / (19 * x); }, 0.01, 0.8, tol, &itr);

  itr2 = 100;
  xx = boost::math::tools::toms748_solve(
    [](double x) { return (20 * x - 1) / (19 * x); }, 0.01, 0.8, tol2, itr2);
  printf(" 5 %11.5f %9i %9i\n", x0, itr, (int)itr2);

  // Test 6
  itr = 100;
  x0 = brent([](double x) { return std::pow(x, 12) - 0.2; }, 0.0, 5.0, tol, &itr);

  itr2 = 100;
  xx = boost::math::tools::toms748_solve(
    [](double x) { return std::pow(x, 12) - 0.2; }, 0.0, 5.0, tol2, itr2);
  printf(" 6 %11.5f %9i %9i\n", x0, itr, (int)itr2);
}
