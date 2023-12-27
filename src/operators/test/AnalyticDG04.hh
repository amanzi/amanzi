/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution: u = t sin(3x) sin(6y)
  Diffusion: K = 1
  Accumulation: a = 0
  Reaction: r = 0
  Velocity: v = [xy(1-x)/2, xy(1-y)/2]
  Source: f = u / t + v . \grad u  
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_04_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_04_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG04 : public AnalyticDGBase {
 public:
  AnalyticDG04(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection){};
  ~AnalyticDG04(){};

  // analytic data in conventional Taylor basis
  // -- diffusion tensor
  virtual Amanzi::WhetStone::Tensor
  Tensor(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

  // -- solution
  virtual void SolutionTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::Polynomial& sol) override
  {
    sol.Reshape(d_, order_, true);
    sol(0, 0) = std::sin(3 * p[0]) * std::sin(6 * p[1]);

    if (order_ > 0) {
      sol(1, 0) = 3 * std::cos(3 * p[0]) * std::sin(6 * p[1]);
      sol(1, 1) = 6 * std::sin(3 * p[0]) * std::cos(6 * p[1]);
    }

    if (order_ > 1) {
      int k = (d_ == 2) ? 2 : 3;
      sol(2, 0) = -4.5 * std::sin(3 * p[0]) * std::sin(6 * p[1]);
      sol(2, 1) = 18.0 * std::cos(3 * p[0]) * std::cos(6 * p[1]);
      sol(2, k) = -18.0 * std::sin(3 * p[0]) * std::sin(6 * p[1]);
    }

    if (order_ > 2) {
      sol(3, 0) = -4.5 * std::cos(3 * p[0]) * std::sin(6 * p[1]);
      sol(3, 1) = -27.0 * std::sin(3 * p[0]) * std::cos(6 * p[1]);
      sol(3, 2) = -54.0 * std::cos(3 * p[0]) * std::sin(6 * p[1]);
      sol(3, 3) = -36.0 * std::sin(3 * p[0]) * std::cos(6 * p[1]);
    }

    if (order_ > 3) {
      sol(4, 0) = 3.375 * std::sin(3 * p[0]) * std::sin(6 * p[1]);
      sol(4, 1) = -27.0 * std::cos(3 * p[0]) * std::cos(6 * p[1]);
      sol(4, 2) = 81.0 * std::sin(3 * p[0]) * std::sin(6 * p[1]);
      sol(4, 3) = -108.0 * std::cos(3 * p[0]) * std::cos(6 * p[1]);
      sol(4, 4) = 54.0 * std::sin(3 * p[0]) * std::sin(6 * p[1]);
    }
    sol.set_origin(p);
    sol *= t;
  }

  // -- accumulation
  virtual void AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p,
                                  double t,
                                  Amanzi::WhetStone::Polynomial& a) override
  {
    a.Reshape(d_, 0, true);
  }

  // -- velocity
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override
  {
    int order_tmp = std::min(order_, 3);
    v.resize(d_);

    for (int i = 0; i < d_; ++i) {
      v[i].Reshape(d_, order_tmp, true);
      v[i].set_origin(p);
    }

    double x(p[0]), y(p[1]);
    v[0](0) = y * (x - x * x) / 2;
    v[1](0) = x * (y - y * y) / 2;

    if (order_ > 0) {
      v[0](1) = y * (1 - 2 * x) / 2;
      v[0](2) = (x - x * x) / 2;

      v[1](1) = (y - y * y) / 2;
      v[1](2) = x * (1 - 2 * y) / 2;
    }

    if (order_ > 1) {
      v[0](3) = -y / 2;
      v[0](4) = (1 - 2 * x) / 2;

      v[1](4) = (1 - 2 * y) / 2;
      v[1](5) = -x / 2;
    }

    if (order_ > 2) {
      v[0](7) = -1.0;
      v[1](8) = -1.0;
    }
  }

  // -- reaction
  virtual void ReactionTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::Polynomial& r) override
  {
    r.Reshape(d_, 0, true);
  }

  // -- source term
  virtual void SourceTaylor(const Amanzi::AmanziGeometry::Point& p,
                            double t,
                            Amanzi::WhetStone::Polynomial& src) override
  {
    Amanzi::WhetStone::Polynomial sol;
    Amanzi::WhetStone::VectorPolynomial v;

    SolutionTaylor(p, 1.0, sol);
    VelocityTaylor(p, t, v);

    v[0].ChangeOrigin(p);
    v[1].ChangeOrigin(p);

    src = sol + (v * Gradient(sol)) * t;
  }
};

#endif
