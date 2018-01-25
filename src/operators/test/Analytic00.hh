/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Polynomial solution and constant coefficient is defined by
  the user-provided gradient and polynomial order.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_00_HH_
#define AMANZI_OPERATOR_ANALYTIC_00_HH_

#include "Polynomial.hh"
#include "VectorPolynomial.hh"

#include "AnalyticBase.hh"

class Analytic00 : public AnalyticBase {
 public:
  Analytic00(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, double gx, double gy, int order) :
      AnalyticBase(mesh),
      poly_(2, order) {
    poly_(1, 0) = gx;
    poly_(1, 1) = gy;

    if (order > 1) {
      poly_(2, 0) = 3.0;
      poly_(2, 1) = 4.0;
      poly_(2, 2) =-3.0;
    }

    if (order > 2) {
      poly_(3, 0) = 1.0;
      poly_(3, 1) = 6.0;
      poly_(3, 2) =-3.0;
      poly_(3, 3) =-2.0;
    }

    grad_.Gradient(poly_);
    rhs_ = poly_.Laplacian();
  }
  ~Analytic00() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return poly_.Value(p);
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    Amanzi::AmanziGeometry::Point v(2);
    v[0] = -grad_[0].Value(p);
    v[1] = -grad_[1].Value(p);
    return v;
  }
 
  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return -velocity_exact(p, t);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return -rhs_.Value(p);
  }

 private:
  Amanzi::WhetStone::Polynomial poly_, rhs_;
  Amanzi::WhetStone::VectorPolynomial grad_;
};

#endif

