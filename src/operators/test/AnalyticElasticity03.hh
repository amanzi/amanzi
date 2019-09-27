/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Quartic displacement which is zero on the boundary.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELASTICITY_03_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELASTICITY_03_HH_

#include "AnalyticElasticityBase.hh"

class AnalyticElasticity03 : public AnalyticElasticityBase {
 public:
  const double lambda = 0.0;
  const double mu = 0.5;
  const double as = 1.0;
 
  AnalyticElasticity03(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticElasticityBase(mesh) {};
  ~AnalyticElasticity03() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(2, 4);
    K(0, 0) = K(1, 1) = lambda + 2 * mu;
    K(0, 1) = K(1, 0) = lambda;
    K(2, 2) = K(3, 3) = 2 * mu;
    return K;
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    double tmp = x * (x - 1.0) * y * (y - 1.0);
    return Amanzi::AmanziGeometry::Point(as * tmp, tmp);
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return 0.0;
  }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];

    double a1 = 2 * as * (lambda + 2 * mu), a2 = 2 * (lambda + 2 * mu);
    double b1 = 2 * as * mu, b2 = 2 * mu;
    double c1 = lambda + mu, c2 = as * (lambda + mu);

    double fx = a1 * y * (y - 1.0) + b1 * x * (x - 1.0) + c1 * (2 * x - 1) * (2 * y - 1);
    double fy = a2 * x * (x - 1.0) + b2 * y * (y - 1.0) + c2 * (2 * x - 1) * (2 * y - 1);
    return Amanzi::AmanziGeometry::Point(-fx, -fy);
  }
};

#endif
