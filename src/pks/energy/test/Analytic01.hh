/*
  This is the energy component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Exact solution is described in "Exact solutions of a nonlinear
  diffusion­convection equation" by V.Vanaja, Phys. Scr. 80 (2009).
  The case m = 2, n = 3.
*/

#include <cmath>

#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "Tensor.hh"

#include "AnalyticBase.hh"

const double vel = 1.0;

class Analytic01 : public AnalyticBase {
 public:
  Analytic01(Teuchos::RCP<const Amanzi::CompositeVector> temp,
             Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : AnalyticBase(mesh), temp_(temp){};
  ~Analytic01(){};

  // Conductivity is m T^{m-1}
  Amanzi::WhetStone::Tensor Conductivity(int c, const Amanzi::AmanziGeometry::Point& p, double t)
  {
    int d = mesh_->space_dimension();
    Amanzi::WhetStone::Tensor K(1, d);

    const Epetra_MultiVector& temp_c = *temp_->ViewComponent("cell");
    K(0, 0) = 2 * temp_c[c][0];
    return K;
  }

  // Fluid velocity is n T^{n-1}
  Amanzi::AmanziGeometry::Point
  FluidVelocity(int c, const Amanzi::AmanziGeometry::Point& p, double t)
  {
    int d = mesh_->space_dimension();
    Amanzi::AmanziGeometry::Point v(d);

    const Epetra_MultiVector& temp_c = *temp_->ViewComponent("cell");
    v[0] = 3 * temp_c[c][0] * temp_c[c][0];
    return v;
  }

  // exact temperature (m = 2, n = 3)
  double temperature_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];

    return (vel / tan(std::pow(vel, 0.5) * (x - vel * t) / 2));
  }

  Amanzi::AmanziGeometry::Point flux_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    int d = mesh_->space_dimension();
    Amanzi::AmanziGeometry::Point v(d);
    return v;
  }

 private:
  Teuchos::RCP<const Amanzi::CompositeVector> temp_;
};
