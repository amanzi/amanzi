/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  This is the operators component of the Amanzi code.

  Non-polynomial solution plus discontinous (scalar) coefficient. This solution
  has discontinuous tangential flux.

  Same as Analytic03, but moves the coefficient from the tensor to the scalar.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_03_HH_
#define AMANZI_OPERATOR_ANALYTIC_03_HH_

#include "AnalyticBase.hh"

class Analytic03 : public AnalyticBase {
 public:
  Analytic03(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : AnalyticBase(mesh)
  {
    k1 = 1.0;
    k2 = 20.0;
    a1 = 1.0 / k1;
    a2 = 1.0 / k2;
    b2 = (a1 - a2) / 4;

    dim = mesh_->getSpaceDimension();
  }
  ~Analytic03(){};

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
  ScalarTensorGradient(const Amanzi::AmanziGeometry::Point& p, double t)
  {
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

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const
  {
    double x = p[0];
    double y = p[1];
    if (x < 0.5) {
      return a1 * x * x + y * y;
    } else {
      return a2 * x * x + y * y + b2;
    }
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
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

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(2);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];

    double plaplace, kmean;
    Amanzi::AmanziGeometry::Point pgrad(dim), kgrad(dim);

    kmean = ScalarDiffusivity(p, t);
    kgrad = ScalarTensorGradient(p, t);
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
