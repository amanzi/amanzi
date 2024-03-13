/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Linear solution for problem with constant tensorial coefficient
  working in 2D and 3D
  Solution: p = x + 2y - gy y
  Diffusion: K = [1   0.1]
                 [0.1   3]
  Velocity: v = [0, 0]
  Source: f = -div(K grad(p))
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_02_HH_
#define AMANZI_OPERATOR_ANALYTIC_02_HH_

#include "AnalyticBase.hh"

class Analytic02 : public AnalyticBase {
 public:
  Analytic02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : AnalyticBase(mesh), g_(0.0), v_(d_)
  {
    v_[0] = 1.0, v_[1] = 2.0;
  }
  Analytic02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, double g)
    : AnalyticBase(mesh), g_(g), v_(d_)
  {
    v_[0] = 1.0, v_[1] = 2.0;
  }
  Analytic02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
             const Amanzi::AmanziGeometry::Point& v,
             double g,
             const Amanzi::WhetStone::Tensor& K)
    : AnalyticBase(mesh), g_(g), K_(K), v_(v){};
  ~Analytic02(){};

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(d_, 2);
    if (K_.size() == 0) {
      K(0, 0) = 1.0;
      K(1, 1) = 3.0;
      K(0, 1) = 0.1;
      K(1, 0) = 0.1;
      if (d_ == 3) K(2, 2) = 1.0;
      return K;
    }
    return K_;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const
  {
    return p * v_ - g_ * p[d_ - 1];
  }

  // Gradient of potential, since the base class does not handle gravity.
  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return v_;
  }

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(d_);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { return 0.0; }

 private:
  double g_;
  Amanzi::WhetStone::Tensor K_;
  Amanzi::AmanziGeometry::Point v_;
};

#endif
