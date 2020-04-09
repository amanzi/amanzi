/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_08B_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_08B_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG08b : public AnalyticDGBase {
 public:
  // cone
  const double phi = 2 * M_PI / 3;
  const Amanzi::AmanziGeometry::Point r0 =
    Amanzi::AmanziGeometry::Point(0.5 * std::cos(phi), 0.5 * std::sin(phi));
  const double R0 = 0.35;

  // hump
  const Amanzi::AmanziGeometry::Point r1 =
    Amanzi::AmanziGeometry::Point(0.5 * std::cos(phi), -0.5 * std::sin(phi));
  const double R1 = 0.35;

  // notched cylinder
  const Amanzi::AmanziGeometry::Point r2 =
    Amanzi::AmanziGeometry::Point(0.5, 0.0);
  const double R2 = 0.4;
  const double W2 = 0.05;

 public:
  AnalyticDG08b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order,
                bool advection)
    : AnalyticDGBase(mesh, order, advection),
      cone_(true),
      hump_(true),
      cylinder_(true),
      hump_amp_(1.0){};
  ~AnalyticDG08b(){};

  // control of details
  void set_shapes(bool cone, bool hump, bool cylinder)
  {
    cone_ = cone;
    hump_ = hump;
    cylinder_ = cylinder;
  }
  void set_amplitudes(double hump_amp) { hump_amp_ = hump_amp; }

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
  virtual void SolutionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& sol) override
  {
    sol.Reshape(d_, order_, true);
    sol.set_origin(p);

    if (cone_) {
      auto dp = p - r0;
      double tmp = Amanzi::AmanziGeometry::norm(dp);
      if (tmp <= R0) {
        sol(0) += 1.0 - tmp / R0;

        if (order_ > 0) {
          double factor = tmp * R0;
          sol(1) += -dp[0] / factor;
          sol(2) += -dp[1] / factor;
        }

        if (order_ > 1) {
          double factor = tmp * tmp * tmp * R0;
          sol(3) += -dp[1] * dp[1] / factor / 2;
          sol(4) += dp[0] * dp[1] / factor;
          sol(5) += -dp[0] * dp[0] / factor / 2;
        }

        if (order_ > 2) AMANZI_ASSERT(false);
      }
    }
    if (hump_) {
      auto dp = p - r1;
      double tmp = Amanzi::AmanziGeometry::norm(dp);
      double sn, cs, arg = M_PI * tmp / R1;

      if (tmp <= R1) {
        cs = std::cos(arg);
        sol(0) = (1.0 + cs) / 4;

        if (order_ > 0) {
          sn = std::sin(arg);
          double factor = M_PI / (4 * R1 * tmp);
          sol(1) += -sn * factor * dp[0];
          sol(2) += -sn * factor * dp[1];
        }

        if (order_ > 1) {
          double fac1 = M_PI / (4 * R1 * tmp * tmp * tmp);
          double fac2 = M_PI * tmp / R1;
          sol(3) +=
            -fac1 * (dp[1] * dp[1] * sn + dp[0] * dp[0] * cs * fac2) / 2;
          sol(4) += fac1 * dp[0] * dp[1] * (sn - cs * fac2);
          sol(5) +=
            -fac1 * (dp[0] * dp[0] * sn + dp[1] * dp[1] * cs * fac2) / 2;
        }

        if (order_ > 2) AMANZI_ASSERT(false);

        if (hump_amp_ != 1.0) sol *= hump_amp_;
      }
    }
    if (cylinder_) {
      auto dp = p - r2;
      double tmp = Amanzi::AmanziGeometry::norm(dp);
      if (tmp <= R2 && !(dp[1] >= 0.0 && std::fabs(dp[0]) <= W2)) {
        sol(0) += 1.0;
      }
    }
  }

  // -- accumulation
  virtual void
  AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                     Amanzi::WhetStone::Polynomial& a) override
  {
    a.Reshape(d_, 0, true);
    a(0, 0) = 1.0;
    a.set_origin(p);
  }

  // -- velocity is not used yet in tests
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override
  {
    v.resize(d_);
    for (int i = 0; i < d_; ++i) {
      v[i].Reshape(d_, 0, true);
      v[i].set_origin(p);
    }
  }

  // -- reaction
  virtual void ReactionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& r) override
  {
    r.Reshape(d_, 0, true);
    r.set_origin(p);
  }

  // -- source term
  virtual void SourceTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                            Amanzi::WhetStone::Polynomial& src) override
  {
    src.Reshape(d_, 0, true);
    src.set_origin(p);
  }

 private:
  bool cone_, hump_, cylinder_;
  double hump_amp_;

 private:
  double StepNotchedCircle_(double x, double y);
};

#endif
