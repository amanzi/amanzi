#ifndef OPERATOR_MARSHAK_TESTCLASS_HH_
#define OPERATOR_MARSHAK_TESTCLASS_HH_


#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "tensor.hh"
#include "mfd3d_diffusion.hh"

namespace Amanzi{

class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : mesh_(mesh),
      TemperatureFloor(0.),
      TemperatureSource(100.) { 
    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh_);
    cvs.SetGhosted(true);
    cvs.SetComponent("cell", AmanziMesh::CELL, 1);
    cvs.SetOwned(false);
    cvs.AddComponent("face", AmanziMesh::FACE, 1);

    values_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
    derivatives_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
  }
  ~HeatConduction() {};

  // main members
  void UpdateValues(const CompositeVector& u) { 
    const Epetra_MultiVector& uc = *u.ViewComponent("cell", true); 
    const Epetra_MultiVector& values_c = *values_->ViewComponent("cell", true); 

    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    for (int c = 0; c < ncells; c++) {
      values_c[0][c] = std::pow(uc[0][c], 3.0);
    }

    derivatives_->PutScalar(1.0);
  }

  double Conduction(int c, double T) const {
    // ASSERT(T > 0.0);
    return T * T * T;
  }

  Teuchos::RCP<CompositeVector> values() { return values_; }
  Teuchos::RCP<CompositeVector> derivatives() { return derivatives_; }
 

  double exact(double t, const Amanzi::AmanziGeometry::Point& p) {
    double x = p[0], c = 0.4;
    double xi = c * t - x;
    return (xi > 0.0) ? std::pow(3 * c * (c * t - x), 1.0 / 3)  : TemperatureFloor;
  }


public:
  double TemperatureSource;
  double TemperatureFloor;

  
private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<CompositeVector> values_, derivatives_;

};

typedef double(HeatConduction::*ModelUpwindFn)(int c, double T) const; 


}  // namespace Amanzi

  

#endif
