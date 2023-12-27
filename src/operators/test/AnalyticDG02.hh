/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution: u = 1 + x^2 + y^2 + z^2
  Diffusion: K = [1 0.5 0; 0.5 2 0; 0 0 1]
  Accumulation: a = 0
  Reaction: r = 0
  Velocity: v = [0.1 + x - x^2, y - y^2, 0]
  Source: exact, e.g. f = -6 for 2D diffusion
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_02_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_02_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG02 : public AnalyticDGBase {
 public:
  AnalyticDG02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection){};
  ~AnalyticDG02(){};

  // diffusion tensor
  virtual Amanzi::WhetStone::Tensor
  Tensor(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    Amanzi::WhetStone::Tensor K(d_, 2);
    if (d_ == 3) {
      K.PutScalar(0.0);
      K(2, 2) = 1.0;
    }
    K(0, 0) = 1.0;
    K(1, 1) = 2.0;
    K(1, 0) = K(0, 1) = 0.5;
    return K;
  }

  // analytic data in conventional Taylor basis
  // -- solution
  virtual void SolutionTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::Polynomial& sol) override
  {
    sol.Reshape(d_, order_, true);
    sol.set_origin(p);

    sol(0, 0) = 1.0 + p * p;

    for (int i = 0; i < d_; ++i) { sol(1, i) = 2.0 * p[i]; }

    sol(2, 0) = 1.0;
    if (d_ == 2) {
      sol(2, 2) = 1.0;
    } else {
      sol(2, 3) = 1.0;
      sol(2, 5) = 1.0;
    }
  }

  // -- accumulation
  virtual void AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p,
                                  double t,
                                  Amanzi::WhetStone::Polynomial& a) override
  {
    a.Reshape(d_, 0, true);
    a.set_origin(p);
  }

  // -- velocity
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override
  {
    double x(p[0]), y(p[1]);

    v.resize(d_);
    if (!advection_) {
      for (int i = 0; i < 2; ++i) { v[i].Reshape(d_, 0, true); }
    } else {
      for (int i = 0; i < 2; ++i) { v[i].Reshape(d_, 2, true); }
      v[0](0, 0) = 0.1 + x - x * x;
      v[0](1, 0) = 1.0 - 2 * x;
      v[0](2, 0) = -1.0;

      v[1](0, 0) = y - y * y;
      v[1](1, 1) = 1.0 - 2 * y;

      if (d_ == 2) {
        v[1](2, 2) = -1.0;
      } else {
        v[2].Reshape(d_, 0, true);
        v[1](2, 3) = -1.0;
      }
    }

    v.set_origin(p);
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
    double x(p[0]), y(p[1]);

    if (!advection_) {
      src.Reshape(d_, 0, true);
      src(0, 0) = -6.0;
    } else {
      src.Reshape(d_, 3, true);
      src(0, 0) = -6.0 + 0.2 * x + 2 * x * x - 2 * x * x * x + 2 * y * y * (1.0 - y);
      src(1, 0) = 0.2 + 4 * x - 6 * x * x;
      src(1, 1) = 4 * y - 6 * y * y;
      src(2, 0) = 2.0 - 6 * x;
      src(2, 2) = 2.0 - 6 * y;
      src(3, 0) = -2.0;
      src(3, 3) = -2.0;
    }

    if (d_ == 3) src(0, 0) -= 2.0;

    src.set_origin(p);
  }
};

#endif
