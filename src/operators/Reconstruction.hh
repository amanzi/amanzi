/*
  This is the operator component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_RECONSTRUCTION_HH_
#define AMANZI_RECONSTRUCTION_HH_

#include <vector>

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"


namespace Amanzi {
namespace Operators {

class Reconstruction {  
 public:
  Reconstruction() {};
  Reconstruction(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : mesh_(mesh), field_(Teuchos::null) {};
  ~Reconstruction() {};

  // main members
  virtual void Init(Teuchos::RCP<const Epetra_MultiVector> field) { field_ = field; }
  virtual void Compute() = 0;
  virtual void ApplyLimiter(Teuchos::RCP<Epetra_MultiVector> limiter) {};

  virtual double getValue(int id, const AmanziGeometry::Point& p) = 0;
 
 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const Epetra_MultiVector> field_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
