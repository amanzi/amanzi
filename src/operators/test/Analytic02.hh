/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Linear solution for problem with constant tensorial coefficient.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_02_HH_
#define AMANZI_OPERATOR_ANALYTIC_02_HH_

#include "AnalyticBase.hh"

class Analytic02 : public AnalyticBase {
 public:
  Analytic02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : AnalyticBase(mesh) {};
  ~Analytic02() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(2, 2);
    K(0, 0) = 3.0;
    K(1, 1) = 1.0;
    K(0, 1) = 1.0;
    K(1, 0) = 1.0;
    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    return x + 2 * y;
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    Amanzi::AmanziGeometry::Point v(2);
    v[0] = 1.0;
    v[1] = 2.0;
    return v;
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return 0.0;
  }
};

#endif

