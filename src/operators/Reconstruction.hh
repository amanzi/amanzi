/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
  Reconstruction(){};
  Reconstruction(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : mesh_(mesh), field_(Teuchos::null), component_(0){};
  virtual ~Reconstruction() = default;

  // main members
  virtual void Init(Teuchos::RCP<const Epetra_MultiVector> field,
                    Teuchos::ParameterList& plist)
  {
    field_ = field;
  }
  virtual void Init(Teuchos::RCP<const Epetra_MultiVector> field,
                    Teuchos::ParameterList& plist, int component)
  {
    field_ = field;
    component_ = component;
  }
  virtual void ComputeGradient() = 0;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const Epetra_MultiVector> field_;
  int component_;
};

} // namespace Operators
} // namespace Amanzi

#endif
