/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution: u = 1 + x^2 + y^2
  Diffusion: K = [1 0.5 0; 0.5 2 0; 0 0 1]
  Accumulation: a = 0
  Reaction: r = 0
  Velocity: v = [0.1 + x - x^2, y - y^2, 0]
  Source: f = -6
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
    Amanzi::WhetStone::Tensor K(d_, 2);
    if (d_ == 3) {
      K.PutScalar(0.0);
      K(2, 2) = 1.0;
    }
    K(0, 0) = 1.0;
    K(1, 1) = 2.0;
    K(1, 0) = K(0, 1) = 0.5;
    return K;
  }

  // analytic data in conventional Taylor basis
  // -- solution
  virtual void SolutionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& sol) override {
    sol.Reshape(d_, order_, true); 
    sol(0, 0) = 1.0 + p[0] * p[0] + p[1] * p[1];

    sol(1, 0) = 2.0 * p[0];
    sol(1, 1) = 2.0 * p[1];

    sol(2, 0) = 1.0;
    sol(2, 2) = 1.0;
    sol.set_origin(p);
  }

  // -- accumulation
  virtual void AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                                  Amanzi::WhetStone::Polynomial& a) override {
    a.Reshape(d_, 0, true); 
  }

  // -- velocity
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override {
    v.resize(d_);
    for (int i = 0; i < 2; ++i) {
      v[i].Reshape(d_, 2, true); 
      v[i](1, i) = 1.0;
      v[i](2, 2*i) = -1.0;
    }
    v[0](0, 0) = 0.1;

    if (d_ == 3) { 
      v[2].Reshape(d_, 0, true); 
    }
  }

  // -- reaction
  virtual void ReactionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& r) override {
    r.Reshape(d_, 0, true); 
  }

  // -- source term
  virtual void SourceTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                            Amanzi::WhetStone::Polynomial& src) override {
    src.Reshape(d_, 0, true);
    src(0, 0) = -6.0;
    src.set_origin(p);
  }
};

#endif

