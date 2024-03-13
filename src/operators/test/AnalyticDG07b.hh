/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Analytic solution for level-set algorithms: two expanding and
  merging spheres.

  Solution: u = 0.1 + t - (x^2 + y^2)^0.5      if x < 0.5
            u = 0.1 + t - ((x-1)^2 + y^2)^0.5  otherwise
  Diffusion: K = 1
  Accumulation: a = 1
  Reaction: r = 0
  Velocity: v = [0, 0]
  Source: f = 0
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_07B_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_07B_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG07b : public AnalyticDGBase {
 public:
  AnalyticDG07b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection){};
  ~AnalyticDG07b(){};

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
    sol.set_origin(p);

    double dx, dy, dist, dist2, dist3, dist5;
    dx = (p[0] < 0.5) ? p[0] : p[0] - 1.0;
    dy = p[1];
    dist2 = dx * dx + dy * dy;
    dist = std::pow(dist2, 0.5);

    sol(0, 0) = 0.1 + t - dist;

    if (order_ > 0) {
      sol(1, 0) = -dx / dist;
      sol(1, 1) = -dy / dist;
    }

    if (order_ > 1) {
      dist3 = dist2 * dist;
      sol(2, 0) = -dy * dy / dist3 / 2;
      sol(2, 1) = dx * dy / dist3;
      sol(2, 2) = -dx * dx / dist3 / 2;
    }

    if (order_ > 2) {
      dist5 = dist3 * dist2;
      sol(3, 0) = dx * dy * dy / dist5 / 2;
      sol(3, 1) = dy * (dy * dy - 2 * dx * dx) / dist5 / 2;
      sol(3, 2) = dx * (dx * dx - 2 * dy * dy) / dist5 / 2;
      sol(3, 3) = dy * dx * dx / dist5 / 2;
    }

    if (order_ > 3) AMANZI_ASSERT(false);
  }

  // -- accumulation
  virtual void AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p,
                                  double t,
                                  Amanzi::WhetStone::Polynomial& a) override
  {
    a.Reshape(d_, 0, true);
    a(0, 0) = 1.0;
    a.set_origin(p);
  }

  // -- velocity is defined as v = -grad u / |grad u|
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override
  {
    v.resize(d_);
    for (int i = 0; i < d_; ++i) v[i].Reshape(d_, order_ - 1, true);
    v.set_origin(p);

    double dx, dy, dist, dist2, dist3, dist5;
    dx = (p[0] < 0.5) ? p[0] : p[0] - 1.0;
    dy = p[1];
    dist = std::pow(dx * dx + dy * dy, 0.5);

    v[0](0) = dx / dist;
    v[1](0) = dy / dist;

    if (order_ > 1) {
      dist2 = dist * dist;
      dist3 = dist2 * dist;

      v[0](1) = dy * dy / dist3;
      v[0](2) = -dx * dy / dist3;

      v[1](1) = -dx * dy / dist3;
      v[1](2) = dx * dx / dist3;
    }

    if (order_ > 2) {
      dist5 = dist3 * dist2;
      v[0](3) = -3 * dx * dy * dy / (2 * dist5);
      v[0](4) = dy * (2 * dx * dx - dy * dy) / dist5;
      v[0](5) = dx * (2 * dy * dy - dx * dx) / (2 * dist5);

      v[1](3) = dy * (2 * dx * dx - dy * dy) / (2 * dist5);
      v[1](4) = dx * (2 * dy * dy - dx * dx) / dist5;
      v[1](5) = -3 * dx * dx * dy / (2 * dist5);
    }

    if (order_ > 3) AMANZI_ASSERT(false);
  }

  // -- reaction
  virtual void ReactionTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::Polynomial& r) override
  {
    r.Reshape(d_, 0, true);
    r.set_origin(p);
  }

  // -- source term
  virtual void SourceTaylor(const Amanzi::AmanziGeometry::Point& p,
                            double t,
                            Amanzi::WhetStone::Polynomial& src) override
  {
    src.Reshape(d_, 0, true);
    src.set_origin(p);
  }
};

#endif
