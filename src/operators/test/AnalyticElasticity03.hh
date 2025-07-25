/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Quartic displacement which is zero on the boundary.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELASTICITY_03_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELASTICITY_03_HH_

#include "AnalyticElasticityBase.hh"

class AnalyticElasticity03 : public AnalyticElasticityBase {
 public:
  AnalyticElasticity03(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                       double mu = 1.0,
                       double lambda = 0.0,
                       bool flag = false)
    : AnalyticElasticityBase(mesh), mu_(mu), lambda_(lambda), a_(2.0) {};
  ~AnalyticElasticity03() {};

  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 4);
    K(0, 0) = K(1, 1) = 2 * mu_ + lambda_;
    K(1, 0) = K(0, 1) = lambda_;
    K(2, 2) = K(3, 3) = 2 * mu_;
    return K;
  }

  virtual Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p,
                                                       double t)
  {
    double x = p[0];
    double y = p[1];
    double tmp = x * (x - 1.0) * y * (y - 1.0);
    return Amanzi::AmanziGeometry::Point(tmp, a_ * tmp);
  }

  virtual double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { return 0.0; }

  virtual Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p,
                                                     double t)
  {
    double x = p[0];
    double y = p[1];
    double c1(2 * mu_ + lambda_);
    double c2 = (mu_ + lambda_) * (2 * x - 1) * (2 * y - 1);

    double fx = 2 * c1 * y * (y - 1.0) + 2 * mu_ * x * (x - 1.0) + a_ * c2;
    double fy = 2 * a_ * c1 * x * (x - 1.0) + 2 * a_ * mu_ * y * (y - 1.0) + c2;
    return Amanzi::AmanziGeometry::Point(-fx, -fy);
  }

  virtual Amanzi::WhetStone::Tensor stress_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor T(2, 2);
    T(0, 0) = 2 * mu_ - lambda_; // FIXME
    T(0, 1) = T(1, 0) = 2 * mu_;
    T(1, 1) = -4 * mu_ - lambda_;
    return T;
  }

  virtual double volumetric_strain_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    return (2 * x - 1.0) * y * (y - 1.0) + a_ * x * (x - 1.0) * (2 * y - 1.0);
  }

 private:
  double mu_, lambda_, a_;
};

#endif
