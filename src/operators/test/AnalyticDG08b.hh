/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Analytic solution is the step function consisting of three
  geometric shapes: cone centered at x0 of radius r0, hump 
  centered at x1 of radius r1, and notched cylinder centered 
  at x2 of radius r2. The width of the notch is w = r / 4.

  Solution: u = step function(x0, x1, x2, r, w)
  Diffusion: K = 1
  Accumulation: a = 1
  Reaction: r = 0
  Velocity: v = (R(2 \pi t) - I) [x, y]
  Source: f = 0
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_08B_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_08B_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG08b : public AnalyticDGBase {
 public:
  // cone
  const bool shape0 = true;
  const double phi = 2 * M_PI / 3;
  const Amanzi::AmanziGeometry::Point r0 = Amanzi::AmanziGeometry::Point(0.5 * std::cos(phi), 0.5 * std::sin(phi));
  const double R0 = 0.35;

  // hump
  const bool shape1 = true;
  const Amanzi::AmanziGeometry::Point r1 = Amanzi::AmanziGeometry::Point(0.5 * std::cos(phi), -0.5 * std::sin(phi));
  const double R1 = 0.35;

  // notched cylinder
  const bool shape2 = true;
  const Amanzi::AmanziGeometry::Point r2 = Amanzi::AmanziGeometry::Point(0.5, 0.0);
  const double R2 = 0.4;
  const double W2 = 0.05;

 public:
  AnalyticDG08b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection) {};
  ~AnalyticDG08b() {};

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

    bool ok(false);

    if (shape0) {
      auto dp = p - r0;
      double tmp = Amanzi::AmanziGeometry::norm(dp);
      if (tmp <= R0) {
        sol(0) += 1.0 - tmp / R0;

        if (order_ > 0) {
          double factor = tmp * R0;
          sol(1) += -dp[0] / factor;
          sol(2) += -dp[1] / factor;
        }

        // if (order_ > 1) AMANZI_ASSERT(false);
        ok = true;
      }
    }
    if (shape1) {
      auto dp = p - r1;
      double tmp = Amanzi::AmanziGeometry::norm(dp);
      double arg = M_PI * tmp / R1;

      if (tmp <= R1) {
        sol(0) = (1.0 + std::cos(arg)) / 4;

        if (order_ > 0) {
          double sn = std::sin(arg);
          double factor = M_PI / (4 * R1 * tmp);
          sol(1) += -sn * factor * dp[0];
          sol(2) += -sn * factor * dp[1];
        }

        // if (order_ > 1) AMANZI_ASSERT(false);
        ok = true;
      }
    }
    if (shape2) {
      auto dp = p - r2;
      double tmp = Amanzi::AmanziGeometry::norm(dp);
      if (tmp <= R2 && !(dp[1] >= 0.0 && std::fabs(dp[0]) <= W2)) {
        sol(0) += 1.0;
        ok = true;
      }
    }
  }

  // -- accumulation
  virtual void AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                                  Amanzi::WhetStone::Polynomial& a) override {
    a.Reshape(d_, 0, true); 
    a(0, 0) = 1.0;
    a.set_origin(p);
  }

  // -- velocity is not used yet in tests
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override {
    v.resize(d_);
    for (int i = 0; i < d_; ++i) {
      v[i].Reshape(d_, 0, true); 
      v[i].set_origin(p);
    }
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
  double StepNotchedCircle_(double x, double y);
};

#endif

