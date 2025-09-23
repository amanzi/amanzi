/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  The trigonometric solution.
  Solution: p = x - x0(y) & 1 / K (x - x0(y)), x0 = 0.5 + 0.04 sin(4 pi y)
  Diffusion: K = 1 & 100
  Source: f = -Laplacian(p)
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_08_HH_
#define AMANZI_OPERATOR_ANALYTIC_08_HH_

#include "AnalyticBase.hh"

const double AMP = 0.04;

class Analytic08 : public AnalyticBase {
 public:
  Analytic08(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order = 1, double K = 100.0)
    : AnalyticBase(mesh), K_(K) {};
  ~Analytic08(){};

  // diffusivity
  virtual Amanzi::WhetStone::Tensor
  TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    double x = p[0];
    double y = p[1];
    double xref = x / (1.0 + 2 * AMP * std::sin(4 * M_PI * y));

    Amanzi::WhetStone::Tensor K(d_, 1);
    K(0, 0) = (xref < 0.5) ? 1.0 : K_;
    return K;
  }

  // exact solution
  virtual double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const override
  { 
    double x = p[0];
    double y = p[1];
    double xref = x / (1.0 + 2 * AMP * std::sin(4 * M_PI * y));

    double f = x - (0.5 + AMP * std::sin(4 * M_PI * y));
    return (xref < 0.5) ? f : f / K_;
  }

  virtual Amanzi::AmanziGeometry::Point
  gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    double x = p[0];
    double y = p[1];
    double xref = x / (1.0 + 2 * AMP * std::sin(4 * M_PI * y));

    Amanzi::AmanziGeometry::Point grad(1.0, -4 * M_PI * AMP * std::cos(4 * M_PI * y));
    return (xref < 0.5) ? grad : grad / K_;
  }

  Amanzi::AmanziGeometry::Point
  advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    return Amanzi::AmanziGeometry::Point(2);
  }

  virtual double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) override
  {
    double y = p[1];
    return -AMP * 16 * M_PI * M_PI * std::sin(4 * M_PI * y);
  }

 private:
  double K_;
};

#endif
