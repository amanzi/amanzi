/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Constant function: f = 1 + x^2 - y^2.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_02_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_02_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG02 : public AnalyticDGBase {
 public:
  AnalyticDG02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order)
    : AnalyticDGBase(mesh, order) {};
  ~AnalyticDG02() {};

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
    coefs(0, 0) = 1.0 + p[0] * p[0] - p[1] * p[1];

    coefs(1, 0) = 2.0 * p[0];
    coefs(1, 1) =-2.0 * p[1];

    coefs(2, 0) = 1.0;
    coefs(2, 2) =-1.0;
  }

  // source term
  virtual void SourceTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                            Amanzi::WhetStone::Polynomial& src) override {
    src.Reshape(d_, 0, true);
  }
};

#endif

