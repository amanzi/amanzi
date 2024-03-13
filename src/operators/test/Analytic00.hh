/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Polynomial solution and constant coefficient is defined by
  the user-provided gradient and polynomial order:
  Solution: p = 1  order=0
            p = 1 + x + 2y  order=1
            p = 1 + x + 2y + 3x^2 + 4xy - 3y^2  order=2
            p = 1 + x + 2y + 3x^2 + 4xy - 3y^2 + x^3 + 6x^2y - 3xy^2 - 3y^3  order=3
  Diffusion: K = 1
  Velocity: v = [vx, vy]
  Source: f = -Laplacian(p)
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_00_HH_
#define AMANZI_OPERATOR_ANALYTIC_00_HH_

#include "Polynomial.hh"
#include "VectorObjects.hh"

#include "AnalyticBase.hh"

class Analytic00 : public AnalyticBase {
 public:
  Analytic00(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
             int order,
             double g = 0.0,
             const Amanzi::AmanziGeometry::Point v = Amanzi::AmanziGeometry::Point(2),
             double K = 1.0,
             double kr = 1.0)
    : AnalyticBase(mesh), g_(g), K_(K), kr_(kr), v_(v), poly_(2, order)
  {
    poly_(0, 0) = 1.0;

    if (order > 0) {
      poly_(1, 0) = 1.0;
      poly_(1, 1) = 2.0;
    }

    if (order > 1) {
      poly_(2, 0) = 3.0;
      poly_(2, 1) = 4.0;
      poly_(2, 2) = -3.0;
    }

    if (order > 2) {
      poly_(3, 0) = 1.0;
      poly_(3, 1) = 6.0;
      poly_(3, 2) = -3.0;
      poly_(3, 3) = -2.0;
    }

    grad_ = Gradient(poly_);
    // grad_[0] += v_[0] * poly_;
    // grad_[1] += v_[1] * poly_;
    grad_[1](0) += g_;

    Amanzi::WhetStone::VectorPolynomial tmp(2, 2);
    for (int i = 0; i < 2; ++i) { tmp[i] = v_[i] * poly_; }
    rhs_ = kr_ * (Amanzi::WhetStone::Divergence(tmp) - poly_.Laplacian());
  }
  ~Analytic00(){};

  // diffusivity
  virtual Amanzi::WhetStone::Tensor
  TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = K_;
    return K;
  }

  virtual double ScalarDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    return kr_;
  }

  // exact solution
  virtual double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override
  {
    return poly_.Value(p);
  }

  virtual Amanzi::AmanziGeometry::Point
  velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    Amanzi::AmanziGeometry::Point v(2);
    v[0] = -grad_[0].Value(p);
    v[1] = -grad_[1].Value(p);
    return kr_ * (K_ * v);
  }

  virtual Amanzi::AmanziGeometry::Point
  gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    return -velocity_exact(p, t);
  }

  virtual Amanzi::AmanziGeometry::Point
  advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    return v_;
  }

  virtual double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    return rhs_.Value(p);
  }

 private:
  double g_, K_, kr_;
  Amanzi::AmanziGeometry::Point v_;
  Amanzi::WhetStone::Polynomial poly_, rhs_;
  Amanzi::WhetStone::VectorPolynomial grad_;
};

#endif
