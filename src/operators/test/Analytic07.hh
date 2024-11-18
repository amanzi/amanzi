/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Trigonometric solution
  Solution: p = cos(a * x * y)
  Diffusion: K = 1, k = 1
  Velocity: v = [0, 0]
  Source: f = -div(K grad(p))
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_07_HH_
#define AMANZI_OPERATOR_ANALYTIC_07_HH_

#include "AnalyticBase.hh"

const double a = 6.1;

class Analytic07 : public AnalyticBase {
 public:
  Analytic07(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : AnalyticBase(mesh) {};
  ~Analytic07(){};

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

  double ScalarDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return 1.0;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const
  {
    double x = p[0];
    double y = p[1];
    return std::cos(a * x * y);
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];

    return -a * Amanzi::AmanziGeometry::Point(std::sin(a * x * y), std::sin(a * x * y));
  }

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(2);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];

    return a * a * (x * x + y * y) * std::cos(a * x * y);
  }
};

#endif
