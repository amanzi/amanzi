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
  Solution: p = sin(PI x) sin(2 PI y)
  Diffusion: K = R [1 0.1; 0.1 3] R^T
  Velocity: v = [0, 0]
  Source: f = -div(K grad(p))
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_09_HH_
#define AMANZI_OPERATOR_ANALYTIC_09_HH_

#include "AnalyticBase.hh"

class Analytic09 : public AnalyticBase {
 public:
  Analytic09(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, double theta = 0.0)
    : AnalyticBase(mesh), theta_(theta) {};
  ~Analytic09() {};

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 2);
    K(0, 0) = 1.0;
    K(0, 1) = K(1, 0) = 0.1;
    K(1, 1) = 100.0;

    K.Rotate(&theta_);
    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const
  {
    double x = p[0];
    double y = p[1];
    return sin(2 * x) * sin(3 * y);
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];

    double px = 2 * cos(2 * x) * sin(3 * y);
    double py = 3 * sin(2 * x) * cos(2 * y);
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
    auto K =  TensorDiffusivity(p, t);
    return (4 * K(0, 0) + 9 * K(1, 1)) * std::sin(2 * x) * std::sin(3 * y) - 12 * K(0, 1) * std::cos(2 * x) * std::cos(3 * y);
  }

 private:
  double theta_;
};

#endif
