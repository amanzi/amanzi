/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone
*/

#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "UnitTest++.h"

// Amanzi::WhetStone
#include "SplineCubicNotAKnot1D.hh"
#include "SplineCubicNotAKnot2D.hh"
#include "TensorCubicBSpline2D.hh"

double f_exact(double x) {
  return 1.0 - 2.0 * x + 0.5 * x * x + 3.0 * x * x * x;
}

double df_exact(double x) {
  return -2.0 + x + 9.0 * x * x;
}

double d2f_exact(double x) {
  return 1.0 + 18.0 * x;
}


TEST(CUBIC_EXACTNESS_1D)
{
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  std::vector<double> x = { -1.0, -0.4, 0.1, 0.9, 1.7, 2.5 };

  std::vector<double> y(x.size());
  for (int i = 0; i < x.size(); ++i) y[i] = f_exact(x[i]);

  SplineCubicNotAKnot1D spline(x, y);

  std::vector<double> test_xp = { -1.0, -0.8, -0.4, -0.15, 0.1,
                                   0.45, 0.9, 1.3, 1.7, 2.1, 2.5 };

  double tol = 1.0e-12;
  for (double xp : test_xp) {
    auto out = spline.evaluate(xp);

    double f_err   = std::fabs(out[0] - f_exact(xp));
    double df_err  = std::fabs(out[1] - df_exact(xp));
    double d2f_err = std::fabs(out[2] - d2f_exact(xp));

    CHECK_CLOSE(f_err, 0.0, tol);
    CHECK_CLOSE(df_err, 0.0, tol);
    CHECK_CLOSE(d2f_err, 0.0, tol);
  }
}


double g_exact(double d, double t) {
  return 1.0
       + 2.0 * d - 3.0 * t
       + 0.5 * d * d - 0.25 * t * t + 0.75 * d * t
       + 0.1 * d * d * d - 0.2 * t * t * t + 0.05 * d * d * t - 0.07 * d * t * t
       + 0.03 * d * d * d * t * t;
}

double gd_exact(double d, double t) {
  return 2.0
       + d + 0.75 * t
       + 0.3 * d * d + 0.1 * d * t - 0.07 * t * t
       + 0.09 * d * d * t * t;
}

double gt_exact(double d, double t) {
  return -3.0
        - 0.5 * t + 0.75 * d
        - 0.6 * t * t + 0.05 * d * d - 0.14 * d * t
        + 0.06 * d * d * d * t;
}

double gdd_exact(double d, double t) {
  return 1.0 + 0.6 * d + 0.1 * t + 0.18 * d * t * t;
}

double gdt_exact(double d, double t) {
  return 0.75 + 0.1 * d - 0.14 * t + 0.18 * d * d * t;
}

double gtt_exact(double d, double t) {
  return -0.5 - 1.2 * t - 0.14 * d + 0.06 * d * d * d;
}


TEST(CUBIC_EXACTNESS_2D)
{
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  std::vector<double> x = { -1.0, -0.4, 0.2, 0.9, 1.7, 2.5 };
  std::vector<double> y = { -0.8, -0.2, 0.4, 1.1, 1.9 };

  std::vector<std::vector<double>> values;
  values.resize(y.size(), std::vector<double>(x.size()));

  for (int j = 0; j < y.size(); ++j) {
    for (int i = 0; i < x.size(); ++i) {
      values[j][i] = g_exact(x[i], y[j]);
    }
  }

  SplineCubicNotAKnot2D spline(x, y, values);

  std::vector<double> test_xp = { -1.0, -0.75, -0.4, -0.1, 0.2, 0.55, 0.9, 1.3, 1.7, 2.1, 2.5 };
  std::vector<double> test_yp = { -0.8, -0.5, -0.2, 0.1, 0.4, 0.75, 1.1, 1.5, 1.9 };

  double tol = 1.0e-10;
  for (double xp : test_xp) {
    for (double yp : test_yp) {
      auto out = spline.evaluate(xp, yp);

      double exact[6] = { g_exact(xp, yp),
                          gd_exact(xp, yp),
                          gt_exact(xp, yp),
                          gdd_exact(xp, yp),
                          gdt_exact(xp, yp),
                          gtt_exact(xp, yp) };

      for (int k = 0; k < 6; ++k) {
        double err = std::fabs(out[k] - exact[k]);
        CHECK_CLOSE(err, 0.0, tol);
      }
    }
  }
}


TEST(CUBIC_BSPLINE_2D)
{
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  std::vector<double> x = { -1.0, -0.4, 0.2, 0.9, 1.7, 2.5 };
  std::vector<double> y = { -0.8, -0.2, 0.4, 1.1, 1.9 };

  int nx = x.size();
  int ny = y.size();

  std::vector<std::vector<double>> values(ny, std::vector<double>(nx));

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
       // values[j][i] = std::sin(x[i]) * std::cos(y[j]);
       values[j][i] = g_exact(x[i], y[j]);
    }
  }

  TensorCubicBSpline2D spline(x, y);
  spline.build(values);

  std::vector<double> test_xp = { -1.0, -0.75, -0.4, -0.1, 0.2, 0.55, 0.9, 1.3, 1.7, 2.1, 2.5 };
  std::vector<double> test_yp = { -0.8, -0.5, -0.2, 0.1, 0.4, 0.75, 1.1, 1.5, 1.9 };

  double tol = 1.0e-10;
  for (double xp : test_xp) {
    for (double yp : test_yp) {
      auto out = spline.evaluate(xp, yp);

      double exact[6] = { g_exact(xp, yp),
                          gd_exact(xp, yp),
                          gt_exact(xp, yp),
                          gdd_exact(xp, yp),
                          gdt_exact(xp, yp),
                          gtt_exact(xp, yp) };

      for (int k = 0; k < 6; ++k) {
        double err = std::fabs(out[k] - exact[k]);
        CHECK_CLOSE(0.0, err, tol);
      }
    }
  }
}
