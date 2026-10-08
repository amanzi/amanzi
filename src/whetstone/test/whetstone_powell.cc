/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <vector>
#include <iostream>

#include "UnitTest++.h"

#include "DenseVector.hh"
#include "PowellHybrid.hh"

using namespace Amanzi::Utils;

class Vector : public Amanzi::WhetStone::DenseVector {
 public:
  Vector(int n) { Reshape(n, 0.0); }

  double& operator[](int i) { return this->operator()(i); }
  const double& operator[](int i) const { return this->operator()(i); }

  friend Vector operator*(double c, const Vector& v)
  {
    Vector r(v);
    for (int i = 0; i < v.size(); ++i) r[i] *= c;
    return r;
  }

  friend Vector operator+(const Vector& u, const Vector& v)
  {
    Vector r(u);
    for (int i = 0; i < v.size(); ++i) r[i] += v[i];
    return r;
  }

  friend Vector operator-(const Vector& u, const Vector& v)
  {
    Vector r(u);
    for (int i = 0; i < v.size(); ++i) r[i] -= v[i];
    return r;
  }

  int size() const { return NumRows(); }
  friend double norm(const Vector& v) { double a; v.Norm2(&a); return a; }
};


TEST(POWELL_HYBRID_1)
{
  auto F1 = [](const Vector& x) {
    Vector r(x.size());
    r[0] = x[0] * x[0] + x[1] - 1.0;
    r[1] = x[0] + x[1] * x[1] - 1.0;
    return r;
  };

  Vector x1(2);
  x1[0] = 2.0;
  x1[1] = 1.5;
  int itrs = 10;
  Vector sol = PowellHybrid(x1, F1, &itrs);

  double exact = (-1.0 + std::sqrt(5.0)) / 2;
  CHECK_CLOSE(exact, sol[0], 1e-8);
  CHECK_CLOSE(exact, sol[1], 1e-8);
  std::cout << "sol: " << sol[0] << " " << sol[1] << ", itrs=" << itrs << "\n";
}


TEST(POWELL_HYBRID_2)
{
  auto F2 = [](const Vector& x) {
    Vector r(x.size());
    r[0] = x[0] * x[0] + x[1] + x[2] - 3.0;
    r[1] = x[0] + x[1] * x[1] + x[2] - 3.0;
    r[2] = x[0] + x[1] + x[2] * x[2] - 3.0;
    return r;
  };

  Vector x0(3);
  x0[0] = 4.0;
  x0[1] =-2.0;
  x0[2] = 0.0;
  int itrs = 10;
  Vector sol = PowellHybrid(x0, F2, &itrs);

  double tmp = std::sqrt(2.0);
  CHECK_CLOSE(1 + tmp, sol[0], 1e-8);
  CHECK_CLOSE(-tmp, sol[1], 1e-8);
  CHECK_CLOSE(-tmp, sol[2], 1e-8);
  std::cout << "sol: " << sol[0] << " " << sol[1] << " " << sol[2] << ", itrs=" << itrs << "\n";
}


TEST(POWELL_HYBRID_3)
{
  auto F3 = [](const Vector& x) {
    Vector r(x.size());
    r[0] = x[0] * x[0] + x[1] + x[2] + x[3] - 4.0;
    r[1] = x[0] + x[1] * x[1] + x[2] + x[3] - 4.0;
    r[2] = x[0] + x[1] + x[2] * x[2] + x[3] - 4.0;
    r[3] = x[0] + x[1] + x[2] + x[3] * x[3] - 4.0;
    return r;
  };

  Vector x0(4);
  x0[0] = 4.0;
  x0[1] =-2.0;
  x0[2] = 0.0;
  x0[3] = 0.0;
  int itrs = 10;
  Vector sol = PowellHybrid(x0, F3, &itrs);

  Vector r = F3(sol);
  CHECK_CLOSE(0.0, r[0], 1e-8);
  CHECK_CLOSE(0.0, r[1], 1e-8);
  CHECK_CLOSE(0.0, r[2], 1e-8);
  CHECK_CLOSE(0.0, r[3], 1e-8);
  std::cout << "sol: " << sol[0] << " " << sol[1] << " " << sol[2] << " " << sol[3] << ", itrs=" << itrs << "\n";

  // repeat
  x0 = sol;
  x0[2] += 1e-5;
  sol = PowellHybrid(x0, F3, &itrs);
  std::cout << "sol: " << sol[0] << " " << sol[1] << " " << sol[2] << " " << sol[3] << ", itrs=" << itrs << "\n";
}
