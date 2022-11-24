/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  This is a 3D extension of the test 06.

  Solution: u = exp(-a [(x - x0)^2 + (y - y0)^2] + (z - z0)^2), a = 20
            default center: x0 = 0.75, y0 = 0.5, z0 = 0.5
  Diffusion: K = 1
  Accumulation: a = 0
  Reaction: r = 0
  Velocity: v = cos(\pi t) [sin(\pi x) cos(\pi y), -cos(\pi x) sin(\pi y)]
  Source: f = 0
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_06C_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_06C_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG06c : public AnalyticDGBase {
 public:
  const double a = 2.0;

 public:
  AnalyticDG06c(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection){};
  ~AnalyticDG06c(){};

  // analytic data in conventional Taylor basis
  // -- diffusion tensor
  virtual Amanzi::WhetStone::Tensor
  Tensor(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    Amanzi::WhetStone::Tensor K(3, 1);
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

    double dx = p[0] - 0.75;
    double dy = p[1] - 0.5;
    double dz = p[2] - 0.5;
    double dist2 = dx * dx + dy * dy + dz * dz;
    double u = std::exp(-a * dist2);

    dx *= -2 * a;
    dy *= -2 * a;
    dz *= -2 * a;

    sol(0) = u;

    if (order_ > 0) {
      sol(1) = u * dx;
      sol(2) = u * dy;
      sol(3) = u * dz;
    }

    double dx2, dy2, dz2;
    if (order_ > 1) {
      dx2 = dx * dx;
      dy2 = dy * dy;
      dz2 = dz * dz;
      sol(2, 0) = u * (dx2 - 2 * a) / 2;
      sol(2, 1) = u * dx * dy;
      sol(2, 2) = u * dx * dz;
      sol(2, 3) = u * (dy2 - 2 * a) / 2;
      sol(2, 4) = u * dy * dz;
      sol(2, 5) = u * (dz2 - 2 * a) / 2;
    }

    if (order_ > 2) AMANZI_ASSERT(false);
  }

  // -- accumulation
  virtual void AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p,
                                  double t,
                                  Amanzi::WhetStone::Polynomial& acc) override
  {
    acc.Reshape(d_, 0, true);
    acc.set_origin(p);
  }

  // -- velocity
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override
  {
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
    v[1](0, 0) = -csx * sny;

    if (order_ > 0) {
      v[0](1, 0) = M_PI * csx * csy;
      v[0](1, 1) = -M_PI * snx * sny;

      v[1](1, 0) = M_PI * snx * sny;
      v[1](1, 1) = -M_PI * csx * csy;
    }

    if (order_ > 1) {
      double factor = M_PI * M_PI;
      v[0](2, 0) = -factor * snx * csy / 2;
      v[0](2, 1) = -factor * csx * sny;
      v[0](2, 3) = -factor * snx * csy / 2;

      v[1](2, 0) = factor * csx * sny / 2;
      v[1](2, 1) = factor * snx * csy;
      v[1](2, 3) = factor * csx * sny / 2;
    }

    if (order_ > 2) AMANZI_ASSERT(false);

    v *= std::cos(M_PI * t);
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
