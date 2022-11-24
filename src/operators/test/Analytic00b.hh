/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  3D version of Analytic00.
  
  Polynomial solution and constant coefficient is defined by
  the user-provided gradient and polynomial order:
  Solution: p = 1  order=0
            p = 1 + gx x + gy y + gz z order=1
  Diffusion: K = 1
  Velocity: v = [vx, vy, vz]
  Source: f = -Laplacian(p)
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_00B_HH_
#define AMANZI_OPERATOR_ANALYTIC_00B_HH_

#include "Polynomial.hh"
#include "VectorObjects.hh"

#include "AnalyticBase.hh"

class Analytic00b : public AnalyticBase {
 public:
  Analytic00b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
              double gx,
              double gy,
              double gz,
              int order,
              const Amanzi::AmanziGeometry::Point v = Amanzi::AmanziGeometry::Point(3))
    : AnalyticBase(mesh), v_(v), poly_(3, order)
  {
    poly_(0, 0) = 1.0;

    if (order > 0) {
      poly_(1, 0) = gx;
      poly_(1, 1) = gy;
      poly_(1, 2) = gz;
    }

    grad_ = Gradient(poly_);

    Amanzi::WhetStone::VectorPolynomial tmp(3, 3);
    for (int i = 0; i < 3; ++i) { tmp[i] = v_[i] * poly_; }
    rhs_ = Amanzi::WhetStone::Divergence(tmp) - poly_.Laplacian();
  }
  ~Analytic00b(){};

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(3, 1);
    K(0, 0) = 1.0;
    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const
  {
    return poly_.Value(p);
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::AmanziGeometry::Point v(3);
    v[0] = -grad_[0].Value(p);
    v[1] = -grad_[1].Value(p);
    v[2] = -grad_[2].Value(p);
    return v;
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return -velocity_exact(p, t);
  }

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return v_;
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { return rhs_.Value(p); }

 private:
  Amanzi::AmanziGeometry::Point v_;
  Amanzi::WhetStone::Polynomial poly_, rhs_;
  Amanzi::WhetStone::VectorPolynomial grad_;
};

#endif
