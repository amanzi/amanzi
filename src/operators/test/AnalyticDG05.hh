/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Nonlinear function: f = sin(3x) sin(6y) * step(x),
  where
     step(x) = 1 if x < 0.5
     step(x) =-1 otherwise.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_05_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_05_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG05 : public AnalyticDGBase {
 public:
  AnalyticDG05(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, advection)
    : AnalyticDGBase(mesh, order, advection) {};
  ~AnalyticDG05() {};

  // analytic solution in conventional Taylor basis
  virtual void TaylorCoefficients(const Amanzi::AmanziGeometry::Point& p, double t,
                                  Amanzi::WhetStone::Polynomial& coefs) override {
    coefs.Reshape(d_, order_, true); 
    coefs(0, 0) = std::sin(3 * p[0]) * std::sin(6 * p[1]);

    if (order_ > 0) {
      coefs(1, 0) = 3 * std::cos(3 * p[0]) * std::sin(6 * p[1]);
      coefs(1, 1) = 6 * std::sin(3 * p[0]) * std::cos(6 * p[1]);
    }

    if (order_ > 1) {
      int k = (d_ == 2) ? 2 : 3;
      coefs(2, 0) =  -4.5 * std::sin(3 * p[0]) * std::sin(6 * p[1]);
      coefs(2, 1) =  18.0 * std::cos(3 * p[0]) * std::cos(6 * p[1]);
      coefs(2, k) = -18.0 * std::sin(3 * p[0]) * std::sin(6 * p[1]);
    }

    // multiply by a step function (1, -1)
    if (p[0] > 0.5) {
      coefs *= -1.0;
    }
  }
};

#endif

