/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Nonlinear fields.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELASTICITY_02_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELASTICITY_02_HH_

#include "AnalyticElasticityBase.hh"

class AnalyticElasticity02 : public AnalyticElasticityBase {
 public:
  AnalyticElasticity02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticElasticityBase(mesh) {};
  ~AnalyticElasticity02() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    return Amanzi::AmanziGeometry::Point(x + y*y*y, x*x*x - y);
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    return -3*x*y + 0.75;
  }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return Amanzi::AmanziGeometry::Point(0.0, 0.0);
  }
};

#endif
