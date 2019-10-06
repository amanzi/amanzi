/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution: u = t (1 + x^2 + y^2)
  Diffusion: K = 0
  Accumulation: a = 1
  Reaction: r = 0
  Velocity: v = [0.1 + x - x^2, y - y^2]
  Source: exact
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_02B_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_02B_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG02b : public AnalyticDGBase {
 public:
  AnalyticDG02b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection) {};
  ~AnalyticDG02b() {};

  // diffusion tensor
  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) override {
    Amanzi::WhetStone::Tensor K(d_, 1);
    K(0, 0) = 0.0;
    return K;
  }

  // analytic data in conventional Taylor basis
  // -- solution
  virtual void SolutionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& sol) override {
    AMANZI_ASSERT(order_ > 1);
    sol.Reshape(d_, order_, true); 
    sol.set_origin(p);

    sol(0) = 1.0 + p * p;
    sol(1) = 2.0 * p[0];
    sol(2) = 2.0 * p[1];
    sol(3) = 1.0;
    sol(5) = 1.0;

    sol *= t;
  }

  // -- accumulation
  virtual void AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                                  Amanzi::WhetStone::Polynomial& a) override {
    a.Reshape(d_, 0, true); 
    a(0) = 1.0;
    a.set_origin(p);
  }

  // -- velocity
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override {
    double x(p[0]), y(p[1]);

    v.resize(d_);
    if (! advection_) {
      for (int i = 0; i < 2; ++i) {
        v[i].Reshape(d_, 0, true); 
      }
    } else {
      for (int i = 0; i < 2; ++i) {
        v[i].Reshape(d_, 2, true); 
      }
      v[0](0, 0) = 0.1 + x - x * x;
      v[0](1, 0) = 1.0 - 2 * x;
      v[0](2, 0) =-1.0;

      v[1](0, 0) = y - y * y;
      v[1](1, 1) = 1.0 - 2 * y;
      v[1](2, 2) =-1.0;
    }

    v.set_origin(p);
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
    double x(p[0]), y(p[1]);

    if (!advection_) {
      src.Reshape(d_, 0, true);
    } else {
      src.Reshape(d_, 3, true);
      src(0, 0) = 2 * x * (0.1 + x - x * x) + 2 * y * y * (1.0 - y);
      src(1, 0) = 0.2 + 4 * x - 6 * x * x;
      src(1, 1) = 4 * y - 6 * y * y;
      src(2, 0) = 2.0 - 6 * x;
      src(2, 2) = 2.0 - 6 * y;
      src(3, 0) = -2.0;
      src(3, 3) = -2.0;
    }

    src *= t;
    src.set_origin(p);

    // add accumulation term whicth equals solution at time 1
    Amanzi::WhetStone::Polynomial sol;
    SolutionTaylor(p, 1.0, sol);
    src += sol;
  }
};

#endif

