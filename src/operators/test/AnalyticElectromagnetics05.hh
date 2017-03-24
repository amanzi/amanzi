/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Convergence analysis in space and time
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_05_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_05_HH_

#include "AnalyticElectromagneticsBase.hh"

class AnalyticElectromagnetics05 : public AnalyticElectromagneticsBase {
 public:
  AnalyticElectromagnetics05(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticElectromagneticsBase(mesh) {};
  ~AnalyticElectromagnetics05() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(p.dim(), 1);
    K(0, 0) = 1.0;
    return K;
  }

  Amanzi::AmanziGeometry::Point electric_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];

    double tmp = std::exp(t + 0.6 * x + 0.8 * y);
    return Amanzi::AmanziGeometry::Point(0.0, 0.0, tmp);
  }

  Amanzi::AmanziGeometry::Point magnetic_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];

    double tmp = std::exp(t + 0.6 * x + 0.8 * y);
    return Amanzi::AmanziGeometry::Point(-0.8 * tmp, 0.6 * tmp, 0.0);
  }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return Amanzi::AmanziGeometry::Point(0.0, 0.0, 0.0);
  }
};

#endif
