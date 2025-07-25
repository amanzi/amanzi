/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Non-polynomial solution and a full non-constant tensor:
  Solution: p = x^3y^2 + x sin(2 PI xy) sin(2 PI y) - gy * y
  Diffusion: K = 1, k = (x+1)^2 + y^2
  Velocity: v = [0, 0]
  Source: f = -div(K grad(p))
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_01C_HH_
#define AMANZI_OPERATOR_ANALYTIC_01C_HH_

#include "AnalyticBase.hh"

class Analytic01c : public AnalyticBase {
 public:
  Analytic01c(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : AnalyticBase(mesh), g_(0.0) {};
  Analytic01c(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, double g)
    : AnalyticBase(mesh), g_(g) {};
  ~Analytic01c() {};

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

  double ScalarDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    return (x + 1) * (x + 1) + y * y;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const
  {
    double x = p[0];
    double y = p[1];
    double xy = x * y;
    return x * xy * xy + x * sin(2 * M_PI * xy) * sin(2 * M_PI * y) - g_ * y;
  }

  // Gradient of potential, since the base class does not handle gravity.
  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];

    double t01, t02, t03, t12, t13;
    double px, py;

    t01 = x * x * y;
    t02 = sin(2 * M_PI * x * y);
    t03 = sin(2 * M_PI * y);

    t12 = cos(2 * M_PI * x * y);
    t13 = cos(2 * M_PI * y);

    px = 3 * y * t01 + t03 * (t02 + 2 * M_PI * y * x * t12);
    py = 2 * x * t01 + x * 2 * M_PI * (x * t12 * t03 + t02 * t13) - g_;

    return Amanzi::AmanziGeometry::Point(px, py);
  }

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(2);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];

    double t01, t02, t03, t12, t13;
    double kr, px, py, pxx, pyy;

    kr = (x + 1) * (x + 1) + y * y;

    t01 = x * x * y;
    t02 = sin(2 * M_PI * x * y);
    t03 = sin(2 * M_PI * y);

    t12 = cos(2 * M_PI * x * y);
    t13 = cos(2 * M_PI * y);

    px = 3 * y * t01 + t03 * (t02 + 2 * M_PI * y * x * t12);
    py = 2 * x * t01 + x * 2 * M_PI * (x * t12 * t03 + t02 * t13) - g_;

    pxx = 6 * x * y * y + 4 * M_PI * t03 * (y * t12 - M_PI * y * y * x * t02);
    pyy =
      2 * x * x * x + x * 4 * M_PI * M_PI * (-x * x * t02 * t03 + 2 * x * t12 * t13 - t02 * t03);

    return -2 * (x + 1) * px - 2 * y * py - kr * (pxx + pyy);
  }

 private:
  double g_;
};

#endif
