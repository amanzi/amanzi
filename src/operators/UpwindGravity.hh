/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Upwind a cell-centered field (e.g. rel perm) using a given
  constant velocity (e.g. gravity).
*/

#ifndef AMANZI_GRAVITY_HH_
#define AMANZI_GRAVITY_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Mesh.hh"
#include "Point.hh"

// Operators
#include "Upwind.hh"

namespace Amanzi {
namespace Operators {

class UpwindGravity : public Upwind {
 public:
  UpwindGravity(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Upwind(mesh), g_(mesh->space_dimension()){};
  ~UpwindGravity(){};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void
  Compute(const CompositeVector& flux, const std::vector<int>& bc_model, CompositeVector& field);

 private:
  int method_, order_;
  double tolerance_;
  AmanziGeometry::Point g_;
};

} // namespace Operators
} // namespace Amanzi

#endif
