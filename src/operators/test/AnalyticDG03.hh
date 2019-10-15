/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_03_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_03_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG03 : public AnalyticDGBase {
 public:
  AnalyticDG03(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order,
               bool advection)
    : AnalyticDGBase(mesh, order, advection){};
  ~AnalyticDG03(){};

  // diffusion tensor
  virtual Amanzi::WhetStone::Tensor
  Tensor(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    Amanzi::WhetStone::Tensor K(d_, 2);
    if (d_ == 3) {
      K.putScalar(0.0);
      K(2, 2) = 1.0;
    }
    K(0, 0) = 1.0;
    K(1, 1) = 2.0;
    K(1, 0) = K(0, 1) = 0.5;

    return K;
  }

  // analytic data in conventional Taylor basis
  // -- solution
  virtual void SolutionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& sol) override
  {
    double x(p[0]), y(p[1]);
    sol.Reshape(d_, order_, true);
    sol(0, 0) = 1.0 + x * x * x + y * y * y + x * y * y;

    sol(1, 0) = 3 * x * x + y * y;
    sol(1, 1) = 3 * y * y + 2 * x * y;

    sol(2, 0) = 3 * x;
    sol(2, 1) = 2 * y;
    sol(2, 2) = 3 * y + x;

    sol(3, 0) = 1.0;
    sol(3, 2) = 1.0;
    sol(3, 3) = 1.0;

    sol.set_origin(p);
  }

  // -- accumulation
  virtual void
  AccumulationTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                     Amanzi::WhetStone::Polynomial& a) override
  {
    a.Reshape(d_, 0, true);
  }

  // -- velocity
  virtual void VelocityTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::VectorPolynomial& v) override
  {
    v.resize(d_);
    for (int i = 0; i < 2; ++i) {
      v[i].Reshape(d_, 2, true);
      v[i](1, i) = 1.0;
      v[i](2, 2 * i) = -1.0;
    }
    v[0](0, 0) = 0.1;

    if (d_ == 3) { v[2].Reshape(d_, 0, true); }
  }

  // -- reaction
  virtual void ReactionTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                              Amanzi::WhetStone::Polynomial& r) override
  {
    r.Reshape(d_, 0, true);
  }

  // -- source term
  virtual void SourceTaylor(const Amanzi::AmanziGeometry::Point& p, double t,
                            Amanzi::WhetStone::Polynomial& src) override
  {
    src.Reshape(d_, 0, true);
    src(0, 0) = 0.0;
    src.set_origin(p);
  }
};

#endif
