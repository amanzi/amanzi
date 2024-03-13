/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  3D version of test AnalyticDG04.

  Solution: u = t sin(3x) sin(6y) sin(4z)
  Diffusion: K = 1
  Accumulation: a = 0
  Reaction: r = 0
  Velocity: v = 0
  Source: f = u / t + v . \grad u
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_04B_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_04B_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG04b : public AnalyticDGBase {
 public:
  AnalyticDG04b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection){};
  ~AnalyticDG04b(){};

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
    double snx(std::sin(3 * p[0])), csx(std::cos(3 * p[0]));
    double sny(std::sin(6 * p[1])), csy(std::cos(6 * p[1]));
    double snz(std::sin(4 * p[2])), csz(std::cos(4 * p[2]));

    sol(0, 0) = snx * sny * snz;

    if (order_ > 0) {
      sol(1, 0) = 3 * csx * sny * snz;
      sol(1, 1) = 6 * snx * csy * snz;
      sol(1, 2) = 4 * snx * sny * csz;
    }

    if (order_ > 1) {
      sol(2, 0) = -4.5 * snx * sny * snz;
      sol(2, 1) = 18.0 * csx * csy * snz;
      sol(2, 2) = 12.0 * csx * sny * csz;
      sol(2, 3) = -18.0 * snx * sny * snz;
      sol(2, 4) = 24.0 * snx * csy * csz;
      sol(2, 5) = -8.0 * snx * sny * csz;
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
    v.resize(d_);
    v.set_origin(p);

    for (int i = 0; i < d_; ++i) { v[i].Reshape(d_, 0, true); }
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
