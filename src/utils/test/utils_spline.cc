#include <iostream>
#include "UnitTest++.h"
 
#include "Spline.hh"


TEST(SPLINE_LINEAR_EXACT) {
  double m = (-1.3 - 2.2) / (1.3 - 0.1);
  Amanzi::Utils::Spline s(0.1, 2.2, m, 1.3, -1.3, m);

  // check endpoints
  CHECK_CLOSE(2.2, s(0.1), 1.e-10);
  CHECK_CLOSE(-1.3, s(1.3), 1.e-10);
  CHECK_CLOSE(m, s.Derivative(0.1), 1.e-10);
  CHECK_CLOSE(m, s.Derivative(1.3), 1.e-10);

  // check others
  CHECK_CLOSE((2.2 -1.3)/2, s((1.3+0.1)/2), 1.e-10);
  double x = 0.1 + 0.2;
  double y = 2.2 + m*0.2;
  CHECK_CLOSE(y, s(x), 1.e-10);

  CHECK_CLOSE(m, s.Derivative(0.5), 1.e-10);

}

TEST(SPLINE_NONLINEAR) {
  double m = (-1.3 - 2.2) / (1.3 - 0.1);
  Amanzi::Utils::Spline s(0.1, 2.2, m, 1.3, -1.3, 0.5*m);

  // check endpoints
  CHECK_CLOSE(2.2, s(0.1), 1.e-10);
  CHECK_CLOSE(-1.3, s(1.3), 1.e-10);
  CHECK_CLOSE(m, s.Derivative(0.1), 1.e-10);
  CHECK_CLOSE(0.5*m, s.Derivative(1.3), 1.e-10);

  // monotonic?
  for (int i=0; i!=1000; ++i) {
    CHECK(s.Derivative(0.1 + (1.3-0.1)*i/1000) < 0.);
  }

  // random points from python implementation
  CHECK_CLOSE(0.9037037037037037, s.Value(0.5), 1.e-10);
  CHECK_CLOSE(-0.7574675325608715, s.Value(1.03333), 1.e-10);
  CHECK_CLOSE(-3.4027777777777786, s.Derivative(0.5), 1.e-10);
  CHECK_CLOSE(-2.5385910493489576, s.Derivative(1.03333), 1.e-10);
  
}
