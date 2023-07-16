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

#include "WRM_vanGenuchten.hh"

using namespace Amanzi;
using namespace Flow;

struct Tol {
  Tol(double eps) : eps_(eps){};
  bool operator()(double a, double b) const { return std::abs(a - b) <= eps_; }
  double eps_;
};

struct F {
  F(double s) : s_(s)
  {
    double m = 0.22;
    double l = 0.5;
    double alpha = 2.0e-3;
    double sr = 0.0;
    std::string krel_function("Mualem");
    double pc0 = 500.0;

    wrm_ = new WRM_vanGenuchten(m, l, alpha, sr, krel_function, pc0);
  }

  double operator()(double pc) const { return wrm_->saturation(pc) - s_; }

  double s_;
  const WRM_vanGenuchten* wrm_;
};

TEST(VAN_GENUCHTEN_BRENT_AND_TOMS748)
{
  int itr, nitr(0), nitr2(0);
  double x0, tol(1e-15);

  boost::uintmax_t itr2;
  Tol tol2(1e-15);
  std::pair<double, double> xx;

  std::cout << "saturation   cap_pressure     Brent   TOMS748\n\n";

  for (double s = 0.15; s < 1.0; s += 0.01) {
    F f(s);

    // Test 1
    itr = 100;
    x0 = Utils::brent(f, 0.0, 1.0e+6, tol, &itr);
    nitr += itr;

    itr2 = 100;
    xx = boost::math::tools::toms748_solve(f, 0.0, 1.0e+6, tol2, itr2);
    nitr2 += itr2;
    printf("%10.3f %14.5f %9i %9i\n", s, x0, itr, (int)itr2);
  }

  printf("\n Total:  %26i %9i\n", nitr, nitr2);
}
