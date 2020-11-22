/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The problem is c E + curl K curl E = Q, where K is the identity tensor.
  We define B = curl E for other tests.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_02_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_02_HH_

#include "AnalyticElectromagneticsBase.hh"

class AnalyticElectromagnetics02 : public AnalyticElectromagneticsBase {
 public:
  AnalyticElectromagnetics02(double c,
                             Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      c_(c),
      AnalyticElectromagneticsBase(mesh) {};
  ~AnalyticElectromagnetics02() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(3, 1);
    K(0, 0) = 1.0;
    return K;
  }

  Amanzi::AmanziGeometry::Point electric_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    double x = p[0];
    double y = p[1];
    double z = p[2];

    double Ex = y * y * sin(z);
    double Ey = z * z * sin(x);
    double Ez = x * x * sin(y);
    return Amanzi::AmanziGeometry::Point(Ex, Ey, Ez);
  }

  Amanzi::AmanziGeometry::Point magnetic_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    double x = p[0];
    double y = p[1];
    double z = p[2];

    double Bx = x * x * cos(y) - 2 * z * sin(x);
    double By = y * y * cos(z) - 2 * x * sin(y);
    double Bz = z * z * cos(x) - 2 * y * sin(z);
    return Amanzi::AmanziGeometry::Point(Bx, By, Bz);
  }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    double x = p[0];
    double y = p[1];
    double z = p[2];

    double fx = -2 * sin(z) + y * y * sin(z);
    double fy = -2 * sin(x) + z * z * sin(x);
    double fz = -2 * sin(y) + x * x * sin(y);
    return Amanzi::AmanziGeometry::Point(fx, fy, fz) + c_ * electric_exact(p, t);
  }

 private:
  double c_; 

};

#endif
