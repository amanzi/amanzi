/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Linear vector field and constant tensor, time-dependent
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_MHD_03_HH_
#define AMANZI_OPERATOR_ANALYTIC_MHD_03_HH_

#include "AnalyticMHD_Base.hh"

class AnalyticMHD_03 : public AnalyticMHD_Base {
 public:
  AnalyticMHD_03(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticMHD_Base(mesh) {};
  ~AnalyticMHD_03() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(3, 2);
    K(0, 0) = 1.0;
    K(1, 1) = 2.0;
    K(2, 2) = 3.0;
    K(0, 1) = K(1, 0) = -0.5;
    K(0, 2) = K(2, 0) = -0.1;
    K(1, 2) = K(2, 1) = -0.2;
    return K;
  }

  Amanzi::AmanziGeometry::Point electric_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    double z = p[2];
    return Amanzi::AmanziGeometry::Point(z + t, x + t, y + t);
  }

  Amanzi::AmanziGeometry::Point magnetic_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return Amanzi::AmanziGeometry::Point(0.0, 0.0, 0.0);
  }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return Amanzi::AmanziGeometry::Point(1.0, 1.0, 1.0);
  }
};

#endif
