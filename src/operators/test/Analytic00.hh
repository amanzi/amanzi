/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_ANALYTIC_00_HH_
#define AMANZI_OPERATOR_ANALYTIC_00_HH_

#include "Polynomial.hh"
#include "VectorPolynomial.hh"

#include "AnalyticBase.hh"

class Analytic00 : public AnalyticBase {
 public:
  Analytic00(int order, double kr, double K, double gravity,
             const Amanzi::AmanziGeometry::Point v=Amanzi::AmanziGeometry::Point(2))
      : AnalyticBase(2),
        order_(order),
        poly_(2, order),
        kr_(kr),
        K_(2, 1),
        gravity_(gravity),
        v_(v)
  {
    K_(0,0) = K;
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

    grav_ = Amanzi::WhetStone::VectorPolynomial(2, 2);
    grav_[0](0,0) = 0.;
    grav_[1](0,0) = -gravity;
 
    grad_ = Gradient(poly_) - grav_;

    if (order > 2) {
      poly_(3, 0) = 1.0;
      poly_(3, 1) = 6.0;
      poly_(3, 2) = -3.0;
      poly_(3, 3) = -2.0;
    }

    rhs_ = Amanzi::WhetStone::Divergence((K_ * grad_) * -kr_);
  }
  ~Analytic00(){};

  std::string name() const override {
    std::stringstream str;
    str << "Analytic00";
    if (gravity_ != 0.0) str << "_Grav";
    if (K_(0,0) != 1.0) str << "_K";
    if (kr_ != 1.0) str << "_kr";
    str << "_PolyOrder" << order_;
    return str.str();
  }

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    return K_;
  }

  double ScalarDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    return kr_;
  }
  
  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override { 
    return poly_.Value(p);
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override { 
    return (K_ * gradient_exact(p, t)) * -kr_;
  }
 
  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override { 
    Amanzi::AmanziGeometry::Point grad(2);
    grad[0] = grad_[0].Value(p);
    grad[1] = grad_[1].Value(p);
    return grad;
  }

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    return v_;
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    return rhs_.Value(p);
  }

 private:
  Amanzi::AmanziGeometry::Point v_;
  Amanzi::WhetStone::Polynomial poly_, rhs_;
  Amanzi::WhetStone::VectorPolynomial grad_, grav_;
  double gravity_;
  int order_;
  double kr_;
  WhetStone::Tensor K_;
};

  
#endif
