/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

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
  Analytic02(int dim) :
      AnalyticBase(dim),
      g_(0.0), v_(d_) { v_[0] = 1.0, v_[1] = 2.0; }
  Analytic02(int dim, double g) :
      AnalyticBase(dim),
      g_(g), v_(d_) { v_[0] = 1.0, v_[1] = 2.0; }
  Analytic02(const Amanzi::AmanziGeometry::Point& v, double g, const Amanzi::WhetStone::Tensor<Kokkos::HostSpace>& K) :
      AnalyticBase(v.dim()),
      g_(g),
      v_(v)
      {
        K_ = K; 
      };
  ~Analytic02() {};

  virtual std::string name() const override { return "Analytic02"; }
  
  const Amanzi::WhetStone::Tensor<Kokkos::HostSpace>&
  TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) const override{
    Amanzi::WhetStone::Tensor<Kokkos::HostSpace> K(d_, 2);
    if (K_.size() == 0) {
      K(0, 0) = 1.0;
      K(1, 1) = 3.0;
      K(0, 1) = 0.1;
      K(1, 0) = 0.1;
      if (d_ == 3) K(2, 2) = 1.0;
      return std::move(K);
    }
    return K_;
  }

  double
  pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override { 
    return p * v_ - g_ * p[d_ - 1];
  }

  // Gradient of potential, since the base class does not handle gravity. 
  Amanzi::AmanziGeometry::Point
  gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override { 
    return v_;
  }

  Amanzi::AmanziGeometry::Point
  advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    return Amanzi::AmanziGeometry::Point(d_);
  }

  double
  source_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    return 0.0;
  }

 private:
  double g_;
  Amanzi::AmanziGeometry::Point v_;
};

#endif

