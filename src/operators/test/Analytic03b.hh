/*
  This is the operators component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)

  Non-polynomial solution plus discontinous (scalar) coefficient. This solution
  has discontinuous tangential flux.

  Same as Analytic03, but moves the coefficient from the tensor to the scalar.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_03B_HH_
#define AMANZI_OPERATOR_ANALYTIC_03B_HH_

#include "AnalyticBase.hh"

class Analytic03b : public AnalyticBase {
 public:
  Analytic03b() :
      AnalyticBase(2)
  {
    K_(0,0) = 1.0;
    k1 = 1.0;
    k2 = 20.0;
    a1 = 1.0 / k1;
    a2 = 1.0 / k2;
    b2 = (a1 - a2) / 4;
  }
  ~Analytic03b() {};

  std::string name() const override { return "Analytic03b"; }
  
  const KOKKOS_INLINE_FUNCTION Amanzi::WhetStone::Tensor<DefaultExecutionSpace>&
  TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    return std::move(K_device_);
  }

  const Amanzi::WhetStone::Tensor<Kokkos::HostSpace>& 
  TensorDiffusivity_host(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    return std::move(K_);
  }

  double ScalarDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    double x = p[0];
    double y = p[1];
    double kr;
    if (x < 0.5) { 
      kr = k1 * (1.0 + x * sin(y));
    } else {
      kr = k2 * (1.0 + 2 * x * x * sin(y));
    }
    return kr;
  }

  // gradient of scalar factor of the tensor
  Amanzi::AmanziGeometry::Point
  ScalarTensorGradient(const Amanzi::AmanziGeometry::Point& p, double t) const {
    double x = p[0];
    double y = p[1];
    Amanzi::AmanziGeometry::Point v(dim);
    if (x < 0.5) { 
      v[0] = k1 * sin(y);
      v[1] = k1 * x * cos(y);
    } else {
      v[0] = k2 * 4 * x * sin(y);
      v[1] = k2 * 2 * x * x * cos(y);
    }
    return v;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override { 
    double x = p[0];
    double y = p[1];
    if (x < 0.5) { 
      return a1 * x * x + y * y;
    } else {
      return a2 * x * x + y * y + b2;
    }
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override { 
    double x = p[0];
    double y = p[1];
    Amanzi::AmanziGeometry::Point v(dim);
    if (x < 0.5) { 
      v[0] = 2 * a1 * x;
      v[1] = 2 * y;
    } else {
      v[0] = 2 * a2 * x;
      v[1] = 2 * y;
    }
    return v;
  }

  Amanzi::AmanziGeometry::Point
  advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override {
    return Amanzi::AmanziGeometry::Point(2);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override { 
    double x = p[0];
    double y = p[1];

    double plaplace, kmean;
    Amanzi::AmanziGeometry::Point pgrad(dim), kgrad(dim);

    kmean = ScalarDiffusivity(p,t);
    kgrad = ScalarTensorGradient(p, t);

    pressure_exact(p, t);
    pgrad = gradient_exact(p, t);

    if (x < 0.5) { 
      plaplace = 2 * (1.0 + a1);
    } else {
      plaplace = 2 * (1.0 + a2);
    }

    return -kgrad * pgrad - kmean * plaplace;
  }

 private:
  int dim;
  double k1, k2;
  double a1, a2, b2;  
};

#endif

