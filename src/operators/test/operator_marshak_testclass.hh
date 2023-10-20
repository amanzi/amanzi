/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef OPERATOR_MARSHAK_TESTCLASS_HH_
#define OPERATOR_MARSHAK_TESTCLASS_HH_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh_Algorithms.hh"

namespace Amanzi {

class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh, double T0)
    : TemperatureSource(100.0), TemperatureFloor(T0), mesh_(mesh)
  {
    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

    values_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
    derivatives_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
  }
  ~HeatConduction(){};

  // main members
  void UpdateValues(const CompositeVector& u,
                    const std::vector<int>& bc_model,
                    const std::vector<double>& bc_value)
  {
    const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);
    Epetra_MultiVector& values_c = *values_->ViewComponent("cell", true);

    int ncells =
      mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
    for (int c = 0; c < ncells; c++) { values_c[0][c] = std::pow(uc[0][c], 3.0); }

    // add boundary face component
    Epetra_MultiVector& values_f = *values_->ViewComponent("face", true);
    for (int f = 0; f != bc_model.size(); ++f) {
      if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
        values_f[0][f] = std::pow(bc_value[f], 3.0);
      }
    }

    derivatives_->PutScalar(1.0);
  }

  double Conduction(int c, double T) const { return T * T * T; }

  Teuchos::RCP<CompositeVector> values() { return values_; }
  Teuchos::RCP<CompositeVector> derivatives() { return derivatives_; }

  double exact(double t, const Amanzi::AmanziGeometry::Point& p)
  {
    double x = p[0], c = 0.4;
    double xi = c * t - x;
    return (xi > 0.0) ? std::pow(3 * c * (c * t - x), 1.0 / 3) : TemperatureFloor;
  }

 public:
  double TemperatureSource;
  double TemperatureFloor;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<CompositeVector> values_, derivatives_;
};

} // namespace Amanzi

#endif
