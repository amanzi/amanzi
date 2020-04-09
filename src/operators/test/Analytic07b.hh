/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Polynomial solution and constant coefficient to test zero boundary
  conditions.

  Solution: p = x (1 - x) y (1 - y) z (1 - z) 
  Diffusion: K = 1
  Velocity: v = 0
  Source: f = -Laplacian(p)
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_07B_HH_
#define AMANZI_OPERATOR_ANALYTIC_07B_HH_

#include "Polynomial.hh"
#include "VectorObjects.hh"

#include "AnalyticBase.hh"

class Analytic07b : public AnalyticBase {
 public:
  Analytic07b(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticBase(mesh) {};
  ~Analytic07b() {};

  Amanzi::WhetStone::Tensor TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K(3, 1);
    K(0, 0) = 1.0;
    return K;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const { 
    double x = p[0];
    double y = p[1];
    double z = p[2];
    return x * (1.0 - x) * y * (1.0 - y) * z * (1.0 - z);
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    double z = p[2];

    Amanzi::AmanziGeometry::Point v(3);
    v[0] = (2 * x - 1.0) * y * (1.0 - y) * z * (1.0 - z);
    v[1] = x * (1.0 - x) * (2 * y - 1.0) * z * (1.0 - z);
    v[2] = x * (1.0 - x) * y * (1.0 - y) * (2 * z - 1.0);
    return v;
  }
 
  Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return -velocity_exact(p, t);
  }

  Amanzi::AmanziGeometry::Point advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) {
    return Amanzi::AmanziGeometry::Point(3);
  }

  double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    double x = p[0];
    double y = p[1];
    double z = p[2];
    return 2 * (y * (1.0 - y) * z * (1.0 - z) + x * (1.0 - x) * z * (1.0 - z) + x * (1.0 - x) * y * (1.0 - y));
  }
};

#endif

