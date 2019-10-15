/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_ANALYTIC_05_HH_
#define AMANZI_OPERATOR_ANALYTIC_05_HH_

#include "AnalyticBase.hh"

class Analytic05 : public AnalyticBase {
 public:
  Analytic05(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : AnalyticBase(mesh){};
  ~Analytic05(){};

  Amanzi::WhetStone::Tensor
  TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    Amanzi::WhetStone::Tensor K(2, 2);
    K(0, 0) = 1.0;
    K(1, 1) = 0.5;
    K(0, 1) = -x * x;
    K(1, 0) = x * x;
    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    return 2 * x + y;
  }

  Amanzi::AmanziGeometry::Point
  gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::AmanziGeometry::Point v(2);
    v[0] = 2.0;
    v[1] = 1.0;
    return v;
  }

  Amanzi::AmanziGeometry::Point
  advection_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(2);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    return 2 * x;
  }
};

#endif
