/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution: u = sqrt((x - x0)^2 + (y - y0)^2) - 0.2
            x0 = 0.5 + 0.25 cos(t)
            y0 = 0.5 + sin(t)
  Diffusion: K = 1
  Accumulation: a = 0
  Reaction: r = 0
  Velocity: v = [0.5 - y, x - 0.5]
  Source: f = 0
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_06_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_06_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG06 : public AnalyticDGBase {
 public:
  AnalyticDG06(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order)
    : AnalyticDGBase(mesh, order) {};
  ~AnalyticDG06() {};

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
    double x0 = 0.5 + 0.25 * std::cos(t);
    double y0 = 0.5 + std::sin(t);

    double dx = p[0] - x0;
    double dy = p[1] - y0;

    double dist, dist2, dist3, dist5;
    dist2 = dx * dx + dy * dy;
    dist = std::pow(dist2, 0.5);
    sol(0, 0) = dist - 0.2;

    if (order_ > 0) {
      sol(1, 0) = dx / dist;
      sol(1, 1) = dy / dist;
    }

    if (order_ > 1) {
      dist3 = dist2 * dist;
      sol(2, 0) =  0.5 * dy * dy / dist3;
      sol(2, 1) =        dx * dy / dist3;
      sol(2, 2) =  0.5 * dx * dx / dist3;
    }

    if (order_ > 2) {
      dist5 = dist2 * dist3;
      sol(3, 0) = -0.5 * dx * dy * dy / dist5;
      sol(3, 1) = dy * (dx * dx - dy * dy / 2) / dist5;
      sol(3, 2) = dx * (dy * dy - dx * dx / 2) / dist5;
      sol(3, 3) = -0.5 * dy * dx * dx / dist5;
    }

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
    for (int i = 0; i < d_; ++i) {
      v[i].Reshape(d_, 1, true); 
      v[i].set_origin(p);
    }
    v[0](0, 0) = 0.5;
    v[0](1, 1) =-1.0;

    v[1](0, 0) =-0.5;
    v[1](1, 0) = 1.0;
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
  }
};

#endif

