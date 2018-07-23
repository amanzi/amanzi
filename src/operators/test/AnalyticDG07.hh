/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution: u = 0.3 + t - [(x - 0.5)^2 + (y - 0.5)^2])^0.5
  Diffusion: K = 1
  Accumulation: a = 1
  Reaction: r = 0
  Velocity: v = [0, 0]
  Source: f = 0
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_07_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_07_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG07 : public AnalyticDGBase {
 public:
  AnalyticDG07(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection) {};
  ~AnalyticDG07() {};

  // analytic data in conventional Taylor basis
  // -- diffusion tensor
  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) override {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

  // -- solution
  virtual void SolutionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& sol) override {
    sol.Reshape(d_, order_, true); 
    sol.set_origin(p);

    double dx = p[0] - 0.5;
    double dy = p[1] - 0.5;
    double dist2 = dx * dx + dy * dy;
    double dist = std::pow(dist2, 0.5);

    sol(0, 0) = 0.3 + t - dist;

    if (order_ > 0) {
      sol(1, 0) = -dx / dist;
      sol(1, 1) = -dy / dist;
    }

    if (order_ > 1) AMANZI_ASSERT(false);
  }

  // -- accumulation
  virtual void AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                                  Amanzi::WhetStone::Polynomial& a) override {
    a.Reshape(d_, 0, true); 
    a(0, 0) = 1.0;
    a.set_origin(p);
  }

  // -- velocity
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override {
    v.resize(d_);
    for (int i = 0; i < d_; ++i) v[i].Reshape(d_, 0, true);
    v.set_origin(p);
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
    src.set_origin(p);
  }

 private:
  double x0_, y0_;
};

#endif

