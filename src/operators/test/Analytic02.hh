/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Linear solution for problem with constant tensorial coefficient
  working in 2D and 3D
  Solution: p = x + 2y - gy y
  Diffusion: K = [3  1]
                 [1  1]
  Velocity: v = [0, 0]
  Source: f = -div(K grad(p))
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_02_HH_
#define AMANZI_OPERATOR_ANALYTIC_02_HH_

#include "AnalyticBase.hh"

class Analytic02 : public AnalyticBase {
 public:
  Analytic02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticBase(mesh),
      g_(0.0) {};
  Analytic02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, double g) :
      AnalyticBase(mesh),
      g_(g) {};
  ~Analytic02() {};

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(d_, 2);
    K(0, 0) = 1.0;
    K(1, 1) = 3.0;
    K(0, 1) = 0.1;
    K(1, 0) = 0.1;
    if (d_ == 3) K(2, 2) = 1.0;

    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];

    return x + 2 * y - g_ * y;
  }

  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    Amanzi::AmanziGeometry::Point v(d_);
    v[0] = 1.0;
    v[1] = 2.0;
    return v;
  }

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) {
    return Amanzi::AmanziGeometry::Point(d_);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { return 0.0; }

 private:
  double g_;
};

#endif

