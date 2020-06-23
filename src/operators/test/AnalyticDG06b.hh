/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution: u = sqrt((x - x0)^2 + (y - y0)^2 + 0.05) - 0.3
            x0 = 0.5 + 0.25 cos(t)
            y0 = 0.5 + 0.25 sin(t)
  Diffusion: K = 1
  Accumulation: a = 0
  Reaction: r = 0
  Velocity: v = [0.5 - y, x - 0.5]
  Source: f = 0
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_06b_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_06b_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG06b : public AnalyticDGBase {
 public:
  const double a = 20.0;

 public:
  AnalyticDG06b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection) {};
  ~AnalyticDG06b() {};

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

    double x0 = 0.5 + 0.25 * std::cos(t);
    double y0 = 0.5 + 0.25 * std::sin(t);

    double dx = p[0] - x0;
    double dy = p[1] - y0;
    double dist2 = dx * dx + dy * dy;
    double u = std::exp(-a * dist2);

    dx *= -2 * a;
    dy *= -2 * a;

    sol(0, 0) = u;

    if (order_ > 0) {
      sol(1, 0) = u * dx;
      sol(1, 1) = u * dy;
    }

    double dx2, dy2;
    if (order_ > 1) {
      dx2 = dx * dx;
      dy2 = dy * dy;
      sol(2, 0) =  u * (dx2 - 2 * a) / 2;
      sol(2, 1) =  u * dx * dy;
      sol(2, 2) =  u * (dy2 - 2 * a) / 2;
    }

    if (order_ > 2) {
      sol(3, 0) = u * (dx2 * dx + 3 * dx2) / 6;
      sol(3, 1) = u * (dx2 * dy - 2 * a * dy) / 2;
      sol(3, 2) = u * (dy2 * dx - 2 * a * dx) / 2;
      sol(3, 3) = u * (dy2 * dy + 3 * dy2) / 6;
    }

    if (order_ > 3) AMANZI_ASSERT(false);
  }

  // -- accumulation
  virtual void AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                                  Amanzi::WhetStone::Polynomial& acc) override {
    acc.Reshape(d_, 0, true); 
    acc.set_origin(p);
  }

  // -- velocity
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override {
    v.resize(d_);
    for (int i = 0; i < d_; ++i) {
      v[i].Reshape(d_, 1, true); 
      v[i].set_origin(p);
    }
    v[0](0, 0) = 0.5 - p[1];
    v[1](0, 0) = p[0] - 0.5;

    v[0](1, 1) =-1.0;
    v[1](1, 0) = 1.0;
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
};

#endif

