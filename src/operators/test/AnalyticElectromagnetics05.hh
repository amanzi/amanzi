/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Convergence analysis in space and time.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_05_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELECTROMAGNETICS_05_HH_

#include "AnalyticElectromagneticsBase.hh"

const double ct = 0.5;

class AnalyticElectromagnetics05 : public AnalyticElectromagneticsBase {
 public:
  AnalyticElectromagnetics05(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticElectromagneticsBase(mesh) {};
  ~AnalyticElectromagnetics05() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(p.dim(), 1);
    K(0, 0) = 1.0 / ct;
    return K;
  }

  Amanzi::AmanziGeometry::Point electric_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    Amanzi::AmanziGeometry::Point a(p.dim()), b(3), e(3);
    OrthonormalSystem(a, b, e);

    double tmp = std::exp(ct * (t + p * a));
    return e * tmp;
  }

  Amanzi::AmanziGeometry::Point magnetic_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    Amanzi::AmanziGeometry::Point a(p.dim()), b(3), e(3);
    OrthonormalSystem(a, b, e);

    double tmp = std::exp(ct * (t + p * a));
    return b * tmp;
  }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    return Amanzi::AmanziGeometry::Point(0.0, 0.0, 0.0);
  }

 private: 
  void OrthonormalSystem(Amanzi::AmanziGeometry::Point& a,
                         Amanzi::AmanziGeometry::Point& b, 
                         Amanzi::AmanziGeometry::Point& e) const {
    if (a.dim() == 2) {
      a = Amanzi::AmanziGeometry::Point(0.6, 0.8);
      b = Amanzi::AmanziGeometry::Point(-0.8, 0.6, 0.0);
      e = Amanzi::AmanziGeometry::Point(0.0, 0.0, 1.0);
    } else {
      double tmp = std::pow(0.71, 0.5);
      a = Amanzi::AmanziGeometry::Point(0.2, 0.5, tmp);
      b = Amanzi::AmanziGeometry::Point(-4 * tmp, tmp, 0.3);

      b /= norm(b); 
      e = a^b;
    }
  }
};

#endif
