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
  AnalyticElasticity03(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticElasticityBase(mesh) {};
  ~AnalyticElasticity03() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    double tmp = x * (x - 1.0) * y * (y - 1.0);
    return Amanzi::AmanziGeometry::Point(tmp, tmp);
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return 0.0;
  }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    double fx = 2 * y * (y - 1.0) + x * (x - 1.0) + (2 * x - 1) * (2 * y - 1) / 2;
    double fy = 2 * x * (x - 1.0) + y * (y - 1.0) + (2 * x - 1) * (2 * y - 1) / 2;
    return Amanzi::AmanziGeometry::Point(-fx, -fy);
  }
};

#endif
