/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Linear vector field and constant tensor.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELASTICITY_01_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELASTICITY_01_HH_

#include "AnalyticElasticityBase.hh"

class AnalyticElasticity01 : public AnalyticElasticityBase {
 public:
  AnalyticElasticity01(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                       double mu = 1.0, double lambda = 0.0, bool flag = false)
    : AnalyticElasticityBase(mesh),
      mu_(mu), lambda_(lambda), flag_(flag) {};
  ~AnalyticElasticity01(){};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    if (lambda_ == 0.0 && !flag_) {
      Amanzi::WhetStone::Tensor K(2, 1);
      K(0, 0) = 2 * mu_;
      return K;
    } else {
      Amanzi::WhetStone::Tensor K(2, 4);
      K(0, 0) = K(1, 1) = 2 * mu_ + lambda_;
      K(1, 0) = K(0, 1) = lambda_;
      K(2, 2) = K(3, 3) = 2 * mu_;
      return K;
    }
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    return Amanzi::AmanziGeometry::Point(x + y, x - 2 * y);
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { return 0.0; }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(0.0, 0.0);
  }

  virtual Amanzi::WhetStone::Tensor stress_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor T(2, 2);
    T(0, 0) = 2 * mu_ - lambda_;
    T(0, 1) = T(1, 0) = 2 * mu_;
    T(1, 1) = -4 * mu_ - lambda_;
    return T;
  }

 private:
  double mu_, lambda_;
  bool flag_;
};

#endif
