/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Identity tensor and non-zero source term.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_02_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_02_HH_

#include "AnalyticElectromagneticsBase.hh"

class AnalyticElectromagnetics02 : public AnalyticElectromagneticsBase {
 public:
  AnalyticElectromagnetics02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : AnalyticElectromagneticsBase(mesh){};
  ~AnalyticElectromagnetics02(){};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(3, 1);
    K(0, 0) = 1.0;
    return K;
  }

  Amanzi::AmanziGeometry::Point electric_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    double z = p[2];
    double Ex = y * y * sin(z);
    double Ey = z * z * sin(x);
    double Ez = x * x * sin(y);
    return Amanzi::AmanziGeometry::Point(Ex, Ey, Ez);
  }

  Amanzi::AmanziGeometry::Point magnetic_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(0.0, 0.0, 0.0);
  }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    double z = p[2];
    double fx = -2 * sin(z) + y * y * sin(z);
    double fy = -2 * sin(x) + z * z * sin(x);
    double fz = -2 * sin(y) + x * x * sin(y);
    return Amanzi::AmanziGeometry::Point(fx, fy, fz);
  }
};

#endif
