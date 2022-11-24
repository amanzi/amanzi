/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  A simple step function moving with constant velocity V. Initial
  location of step is at x = X0:
     u(x,t) = 1 if x < X0 + t V
     u(x,t) = 0 otherwise.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_05_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_05_BASE_HH_

#include "AnalyticDGBase.hh"

class AnalyticDG05 : public AnalyticDGBase {
 public:
  const double X0 = 0.0;
  const double VEL = 0.0;

 public:
  AnalyticDG05(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order, bool advection)
    : AnalyticDGBase(mesh, order, advection){};
  ~AnalyticDG05(){};

  // analytic data in conventional Taylor basis
  // -- diffusion tensor
  virtual Amanzi::WhetStone::Tensor
  Tensor(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 0.0;
    return K;
  }

  // -- solution
  virtual void SolutionTaylor(const Amanzi::AmanziGeometry::Point& p,
                              double t,
                              Amanzi::WhetStone::Polynomial& sol) override
  {
    sol.Reshape(d_, order_, true);
    sol(0) = (p[0] < X0 + t * VEL) ? 1.0 : 0.0;
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
    v.Reshape(d_, d_, 0, true);
    v[0](0) = VEL;
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
    src.Reshape(d_, 0, true);
  }
};

#endif
