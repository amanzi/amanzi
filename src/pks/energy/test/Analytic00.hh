/*
  This is the energy component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)

*/

#include <cmath>

#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "Tensor.hh"

#include "AnalyticBase.hh"

const double vel = 0.0;

class Analytic00 : public AnalyticBase {
 public:
  Analytic00(Teuchos::RCP<const Amanzi::CompositeVector> temp,
             Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) :
      AnalyticBase(mesh),
      temp_(temp) {};
  ~Analytic00() {};

  // Conductivity 
  Amanzi::WhetStone::Tensor Conductivity(int c, const Amanzi::AmanziGeometry::Point& p, double t) {
    int d = mesh_->space_dimension();
    Amanzi::WhetStone::Tensor K(1, d);
    K(0, 0) = 4.0;
    return K;
  }

  // Fluid velocity
  Amanzi::AmanziGeometry::Point FluidVelocity(int c, const Amanzi::AmanziGeometry::Point& p, double t) {
    int d = mesh_->space_dimension();
    Amanzi::AmanziGeometry::Point v(d);
    v[0] = 1.0;
    return v;
  }

  // exact temperature 
  double temperature_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    return p[0] * p[0];
  }

  Amanzi::AmanziGeometry::Point flux_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
    int d = mesh_->space_dimension();
    Amanzi::AmanziGeometry::Point v(d);
    return v;
  }

 private:
  Teuchos::RCP<const Amanzi::CompositeVector> temp_;
};

