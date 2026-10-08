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
                       double mu = 1.0,
                       double lambda = 0.0)
    : AnalyticElasticityBase(mesh), mu_(mu), lambda_(lambda), d_(mesh_->getSpaceDimension()) {};
  ~AnalyticElasticity01() {};

  Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    if (lambda_ == 0.0) {
      Amanzi::WhetStone::Tensor K(d_, 1);
      K(0, 0) = 2 * mu_;
      return K;
    } else {
      Amanzi::WhetStone::Tensor K(d_, 4);
      for (int i = 0; i < d_ * d_; ++i) K(i, i) = 2 * mu_;
      for (int i = 0; i < d_; ++i) {
        for (int j = 0; j < d_; ++j) {
          K(i, j) += lambda_;
        }
      }
      return K;
    }
  }

  Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::AmanziGeometry::Point out(d_);

    double x = p[0];
    double y = p[1];
    out[0] = x + y;
    out[1] = x - 2 * y;

    if (d_ == 3) {
      // double z = p[2];
      // out[2] = x - 2 * y - 3 * z;
    }
    return out;
  }

  double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { return 0.0; }

  Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return Amanzi::AmanziGeometry::Point(d_);
  }

  virtual Amanzi::WhetStone::Tensor stress_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor T(d_, 2);
    if (d_ == 2) {
      T(0, 0) = 2 * mu_ - lambda_;
      T(1, 1) =-4 * mu_ - lambda_;
      T(0, 1) = T(1, 0) = 2 * mu_;
    } else {
      T(0, 0) = 2 * mu_ - 4 * lambda_;
      T(1, 1) =-4 * mu_ - 4 * lambda_;
      T(2, 2) =-6 * mu_ - 4 * lambda_;
      T(0, 1) = T(1, 0) = 2 * mu_;
      T(1, 2) = T(2, 1) =-2 * mu_;
      T(0, 2) = T(2, 0) = mu_;
    }
    return T;
  }

  virtual double volumetric_strain_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return (d_ == 2) ? -1.0 : -4.0;
  }

 private:
  int d_;
  double mu_, lambda_;
};

#endif
