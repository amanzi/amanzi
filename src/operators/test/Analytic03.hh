/*
  This is the operators component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Non-polynomial solution plus discontinous (scalar) coefficient. This solution
  has discontinuous tangential flux:

  Solution:  p = x^2 / k1 + y^2       if x < 0.5,
             p = x^x / k2 + y^2 + b2  otherwise
  Diffusion: K = k1 (1 + x sin(y))      if x < 0.5,
             K = k2 (1 + 2 x^2 sin(y))  otherwise
  Velocity: v = [0, 0]
  Source: f = -div(K grad(p))
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_03_HH_
#define AMANZI_OPERATOR_ANALYTIC_03_HH_

#include "AnalyticBase.hh"

class Analytic03 : public AnalyticBase {
 public:
  Analytic03(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : AnalyticBase(mesh) {
    k1 = 1.0;
    k2 = 20.0;
    a1 = 1.0 / k1;
    a2 = 1.0 / k2;
    b2 = (a1 - a2) / 4;

    dim = mesh_->getSpaceDimension();
  }
  ~Analytic03() {};

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) {
    double x = p[0];
    double y = p[1];
    Amanzi::WhetStone::Tensor K(dim, 1);
    if (x < 0.5) { 
      K(0, 0) = k1 * (1.0 + x * sin(y));
    } else {
      K(0, 0) = k2 * (1.0 + 2 * x * x * sin(y));
    }
    return K;
  }

  // gradient of scalar factor of the tensor
  Amanzi::AmanziGeometry::Point ScalarTensorGradient(const Amanzi::AmanziGeometry::Point& p, double t) {
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

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    double x = p[0];
    double y = p[1];
    if (x < 0.5) { 
      return a1 * x * x + y * y;
    } else {
      return a2 * x * x + y * y + b2;
    }
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
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

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) {
    return Amanzi::AmanziGeometry::Point(2);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];

    double plaplace, kmean;
    Amanzi::AmanziGeometry::Point pgrad(dim), kgrad(dim);

    kmean = (TensorDiffusivity(p, t))(0, 0);
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

