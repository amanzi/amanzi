/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution: u = 1 + x + 2 y + 3 z
  Diffusion: K = 1
  Accumulation: a = 0
  Reaction: r = 0
  Velocity: v = [1.0, 0, 0]
  Source: f = 0 or 1
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_01_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_01_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG01 : public AnalyticDGBase {
 public:
  AnalyticDG01(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection) {};
  ~AnalyticDG01() {};

  // analytic data in conventional Taylor basis
  // -- diffusion tensor
  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) override {
    Amanzi::WhetStone::Tensor K(d_, 1);
    K(0, 0) = 1.0;
    return K;
  }

  // -- solution
  virtual void SolutionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& sol) override {
    sol.Reshape(d_, order_, true); 
    sol.set_origin(p);

    sol(0) = 1.0 + p[0] + 2 * p[1];
    sol(1) = 1.0;
    sol(2) = 2.0;

    if (d_ == 3) {
      sol(0) += 3 * p[2];
      sol(3) += 3.0;
    }
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
    v.set_origin(p);

    for (int i = 0; i < d_; ++i) {
      v[i].Reshape(d_, 0, true); 
    }
    if (advection_) v[0](0, 0) = 1.0;
  }

  // -- reaction
  virtual void ReactionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& r) override {
    r.Reshape(d_, 0, true); 
    r.set_origin(p);
  }

  // -- source term
  virtual void SourceTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                            Amanzi::WhetStone::Polynomial& src) override {
    src.Reshape(d_, 0, true);
    if (advection_) src(0, 0) = 1.0;
    src.set_origin(p);
  }
};

#endif

