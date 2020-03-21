/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution: u = exp(-a [(x - x0)^2 + (y - y0)^2]), a = 20
            default center: x0 = 0.75, y0 = 0.5
  Diffusion: K = 1
  Accumulation: a = 0
  Reaction: r = 0
  Velocity: v = cos(\pi t) [sin(\pi x) cos(\pi y), -cos(\pi x) sin(\pi y)]
  Source: f = 0
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_06_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_06_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG06 : public AnalyticDGBase {
 public:
  const double a = 20.0;

 public:
  AnalyticDG06(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection) {};
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
    sol.set_origin(p);

    double dx = p[0] - 0.75;
    double dy = p[1] - 0.5;
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
      v[i].Reshape(d_, order_, true); 
      v[i].set_origin(p);
    }

    double snx, sny, csx, csy;
    snx = std::sin(M_PI * p[0]);
    sny = std::sin(M_PI * p[1]);

    csx = std::cos(M_PI * p[0]);
    csy = std::cos(M_PI * p[1]);

    v[0](0, 0) = snx * csy;
    v[1](0, 0) =-csx * sny;

    if (order_ > 0) {
      v[0](1, 0) = M_PI * csx * csy;
      v[0](1, 1) =-M_PI * snx * sny;

      v[1](1, 0) = M_PI * snx * sny;
      v[1](1, 1) =-M_PI * csx * csy;
    }

    if (order_ > 1) {
      double factor = M_PI * M_PI;
      v[0](2, 0) =-factor * snx * csy / 2;
      v[0](2, 1) =-factor * csx * sny;
      v[0](2, 2) =-factor * snx * csy / 2;

      v[1](2, 0) = factor * csx * sny / 2;
      v[1](2, 1) = factor * snx * csy;
      v[1](2, 2) = factor * csx * sny / 2;
    }

    if (order_ > 2) {
      double factor = M_PI * M_PI * M_PI;
      v[0](3, 0) =-factor * csx * csy / 6;
      v[0](3, 1) = factor * snx * sny / 2;
      v[0](3, 2) =-factor * csx * csy / 2;
      v[0](3, 3) = factor * snx * sny / 6;

      v[1](3, 0) =-factor * snx * sny / 6;
      v[1](3, 1) = factor * csx * csy / 2;
      v[1](3, 2) =-factor * snx * sny / 2;
      v[1](3, 3) = factor * csx * csy / 6;
    }

    if (order_ > 3) AMANZI_ASSERT(false);

    v *= std::cos(M_PI * t);
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

