/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Relaxation of magnetic lines, trivial steady-state solution 
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_04_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_04_HH_

#include "AnalyticElectromagneticsBase.hh"

class AnalyticElectromagnetics04 : public AnalyticElectromagneticsBase {
 public:
  AnalyticElectromagnetics04(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticElectromagneticsBase(mesh) {};
  ~AnalyticElectromagnetics04() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(p.dim(), 1);
    K(0, 0) = 1.0;
    return K;
  }

  Amanzi::AmanziGeometry::Point electric_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    return Amanzi::AmanziGeometry::Point(0.0, 0.0, 0.0);
  }

  Amanzi::AmanziGeometry::Point magnetic_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    if (t > 0.0) 
      return Amanzi::AmanziGeometry::Point(0.0, 0.0, 1.0);

    double x = p[0];
    double y = p[1];
    double z = (p.dim() == 2) ? 1.0 : p[2];

    double phi = 3.1415926 / 4;
    double tmp = phi * std::exp(-(x * x + y * y) / 2 - z * z / 4);
    return Amanzi::AmanziGeometry::Point(-tmp * y * z, tmp * x * z, 1.0);
  }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    return Amanzi::AmanziGeometry::Point(0.0, 0.0, 0.0);
  }
};

#endif
