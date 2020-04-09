/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_ANALYTIC_ELASTICITY_01_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELASTICITY_01_HH_

#include "AnalyticElasticityBase.hh"

class AnalyticElasticity01 : public AnalyticElasticityBase {
 public:
  AnalyticElasticity01(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : AnalyticElasticityBase(mesh){};
  ~AnalyticElasticity01(){};

  Amanzi::WhetStone::Tensor
  Tensor(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 2.0;
    return K;
  }

  Amanzi::AmanziGeometry::Point
  velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    return Amanzi::AmanziGeometry::Point(x + y, x - y);
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return 0.0;
  }

  Amanzi::AmanziGeometry::Point
  source_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(0.0, 0.0);
  }
};

#endif
