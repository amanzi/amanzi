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
  Solution: p = sin(xyz) - gz * z
  Diffusion: K = [2+x^2   -y     -x ]
                 [ -y    2+y^2   -z ]
                 [ -x     -z   2+z^2]
  Velocity: v = [0, 0, 0]
  Source: f = -div(K grad(p))
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_01B_HH_
#define AMANZI_OPERATOR_ANALYTIC_01B_HH_

#include "AnalyticBase.hh"

class Analytic01b : public AnalyticBase {
 public:
  Analytic01b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : AnalyticBase(mesh), g_(0.0){};
  Analytic01b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, double g)
    : AnalyticBase(mesh), g_(g){};
  ~Analytic01b(){};

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    double z = p[2];
    Amanzi::WhetStone::Tensor K(3, 2);
    K(0, 0) = 2 + x * x;
    K(1, 1) = 2 + y * y;
    K(2, 2) = 2 + z * z;
    K(0, 1) = K(1, 0) = -y;
    K(0, 2) = K(2, 0) = -x;
    K(1, 2) = K(2, 1) = -z;

    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const
  {
    double x = p[0];
    double y = p[1];
    double z = p[2];
    return sin(x * y * z) - g_ * z;
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    double z = p[2];

    double px, py, pz;
    double tmp = cos(x * y * z);
    px = tmp * y * z;
    py = tmp * x * z;
    pz = tmp * x * y - g_;

    return Amanzi::AmanziGeometry::Point(px, py, pz);
  }

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(3);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    double z = p[2];

    double s = sin(x * y * z);
    double c = cos(x * y * z);

    const auto& K = TensorDiffusivity(p, t);
    const auto& v = gradient_exact(p, t);

    double f1 = 2 * (p * v) - v[0] - v[1] - v[2];
    double f2 = -s * (K(0, 0) * y * y * z * z + K(1, 1) * x * x * z * z + K(2, 2) * x * x * y * y);
    f2 += 2 * K(0, 1) * (c * z - s * x * y * z * z);
    f2 += 2 * K(0, 2) * (c * y - s * x * y * y * z);
    f2 += 2 * K(1, 2) * (c * x - s * x * x * y * z);

    return -(f1 + f2);
  }

 private:
  double g_;
};

#endif
