/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Nonlinear function: f = sin(3x) sin(6y).
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_04_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_04_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG04 : public AnalyticDGBase {
 public:
  AnalyticDG04(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order)
    : AnalyticDGBase(mesh, order) {};
  ~AnalyticDG04() {};

  // diffusion tensor
  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) override {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

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

    if (order_ > 2) {
      coefs(3, 0) =  -4.5 * std::cos(3 * p[0]) * std::sin(6 * p[1]);
      coefs(3, 1) = -27.0 * std::sin(3 * p[0]) * std::cos(6 * p[1]);
      coefs(3, 2) = -54.0 * std::cos(3 * p[0]) * std::sin(6 * p[1]);
      coefs(3, 3) = -36.0 * std::sin(3 * p[0]) * std::cos(6 * p[1]);
    }
  }

  // source term
  virtual void SourceTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                            Amanzi::WhetStone::Polynomial& src) override {
    src.Reshape(d_, 0, true);
  }
};

#endif

