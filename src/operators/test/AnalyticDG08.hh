/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Analytic solution is distance function for the notched circle.

  Solution: u = distance_func(x, y; r=0.2, w=0.1)
  Diffusion: K = 1
  Accumulation: a = 1
  Reaction: r = 0
  Velocity: v = [0.5 - y, x - 0.5]
  Source: f = 0
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_08_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_08_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG08 : public AnalyticDGBase {
 public:
  AnalyticDG08(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection){};
  ~AnalyticDG08(){};

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

    double tol(1e-5), x0, y0, dx, dy, sn(std::sin(t)), cs(std::cos(t));
    double uxp, uxm, uyp, uym;
    x0 = 0.5 - 0.25 * sn;
    y0 = 0.5 + 0.25 * cs;
    dx = p[0] - x0;
    dy = p[1] - y0;

    sol(0) = DistanceNotchedCircle_(dx, dy, sn, cs);

    if (order_ > 0) {
      uxm = DistanceNotchedCircle_(dx - tol, dy, sn, cs);
      uxp = DistanceNotchedCircle_(dx + tol, dy, sn, cs);

      uym = DistanceNotchedCircle_(dx, dy - tol, sn, cs);
      uyp = DistanceNotchedCircle_(dx, dy + tol, sn, cs);

      sol(1) = (uxp - uxm) / (2 * tol);
      sol(2) = (uyp - uym) / (2 * tol);
    }

    if (order_ > 1) {
      double umm, ump, upm, upp;
      umm = DistanceNotchedCircle_(dx - tol, dy - tol, sn, cs);
      ump = DistanceNotchedCircle_(dx - tol, dy + tol, sn, cs);
      upm = DistanceNotchedCircle_(dx + tol, dy - tol, sn, cs);
      upp = DistanceNotchedCircle_(dx + tol, dy + tol, sn, cs);

      sol(3) = (uxm - 2 * sol(0) + uxp) / (tol * tol);
      sol(4) = (upp + umm - upm - ump) / (4 * tol * tol);
      sol(5) = (uym - 2 * sol(0) + uyp) / (tol * tol);
    }

    if (order_ > 2) AMANZI_ASSERT(false);
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

  // -- velocity
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override
  {
    v.resize(d_);
    for (int i = 0; i < d_; ++i) {
      v[i].Reshape(d_, 1, true);
      v[i].set_origin(p);
    }
    v[0](0, 0) = 0.5 - p[1];
    v[1](0, 0) = p[0] - 0.5;

    v[0](1, 1) = -1.0;
    v[1](1, 0) = 1.0;
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

 private:
  double DistanceNotchedCircle_(double x0, double y0, double sn, double cs);
};


/* ******************************************************************
* Distance function for counter clockwise rotated notched circle.
* Assumptions: (a) circle center is at point (0, 0); (b) distance is
* positive inside notched circle.
****************************************************************** */
double
AnalyticDG08::DistanceNotchedCircle_(double x0, double y0, double sn, double cs)
{
  // rotate the input point clockwise
  double x = x0 * cs + y0 * sn;
  double y = -x0 * sn + y0 * cs;

  double R(0.2), W(0.05), H;
  H = std::pow(R * R - W * W, 0.5);

  // initial estimate of disdtance - to the perfect circle
  double dx, dy, tmp = std::pow(x * x + y * y, 0.5);
  double dist = R - tmp;

  // cone above notch
  if (y >= H && H * x >= -W * y && H * x <= W * y) {
    dy = y - H;
    if (x >= 0.0)
      dist = -std::pow((x - W) * (x - W) + dy * dy, 0.5);
    else
      dist = -std::pow((x + W) * (x + W) + dy * dy, 0.5);
  }
  // inside notch
  else if (y >= 0.0 && std::fabs(x) <= W) {
    dist = -std::min({ y, W - x, W + x });
  }
  // inside circle
  else if (dist > 0.0) {
    // right-top quarter
    if (x >= W && y >= 0.0) {
      dist = std::min(dist, x - W);
    }
    // left-top quarter
    else if (x <= -W && y >= 0.0) {
      dist = std::min(dist, -W - x);
    }
    // right-bottom quarter
    else if (x >= W && y < 0.0) {
      dx = x - W;
      dist = std::min(dist, std::pow(dx * dx + y * y, 0.5));
    }
    // left-bottom quarter
    else if (x <= -W && y < 0.0) {
      dx = x + W;
      dist = std::min(dist, std::pow(dx * dx + y * y, 0.5));
    }
    // below notch
    else {
      dist = std::min(dist, -y);
    }
  }

  return dist;
}

#endif
