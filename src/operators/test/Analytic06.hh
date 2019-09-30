/*
  This is the operators component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Description is missing.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_06_HH_
#define AMANZI_OPERATOR_ANALYTIC_06_HH_

#include "AnalyticBase.hh"
#include "AnalyticNonlinearCoupled00.hh"

class Analytic06 : public AnalyticBase {
 public:
  Analytic06(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
      : AnalyticBase(mesh),
        ana_(mesh) {}

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) {
    double u = ana_.exact0(p,t);
    double v = ana_.exact1(p,t);
    Amanzi::WhetStone::Tensor K = ana_.Tensor11(p,t);
    double kr = ana_.ScalarCoefficient11(u,v);
    K *= kr;
    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    return ana_.exact1(p,t);
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return ana_.gradient_exact1(p,t);
  }

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) {
    return Amanzi::AmanziGeometry::Point(2);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return ana_.source_exact1(p,t);
  }

 private:
  AnalyticNonlinearCoupled00 ana_;
};

#endif

