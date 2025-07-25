/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Nonlinear fields.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELASTICITY_02_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELASTICITY_02_HH_

#include "AnalyticElasticityBase.hh"

class AnalyticElasticity02 : public AnalyticElasticityBase {
 public:
  AnalyticElasticity02(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : AnalyticElasticityBase(mesh) {};
  ~AnalyticElasticity02() {};

  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

  virtual Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p,
                                                       double t)
  {
    double x = p[0];
    double y = p[1];
    return Amanzi::AmanziGeometry::Point(x + y * y * y, x * x * x - y);
  }

  virtual double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    return -3 * x * y + 0.75;
  }

  virtual Amanzi::WhetStone::Tensor stress_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    Amanzi::WhetStone::Tensor T(2, 2);

    T(0, 0) = 1.0;
    T(0, 1) = T(1, 0) = 1.5 * (y * y + x * x);
    T(1, 1) = -1.0;
    return T;
  }

  virtual double volumetric_strain_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    return 0.0;
  }

  virtual Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p,
                                                     double t)
  {
    return Amanzi::AmanziGeometry::Point(0.0, 0.0);
  }
};

#endif
