/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Upwind a cell-centered field defined on a network of manifolds
  using a given face-based flux.
*/

#ifndef AMANZI_UPWIND_FLUX_MANIFOLDS_HH_
#define AMANZI_UPWIND_FLUX_MANIFOLDS_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Mesh.hh"
#include "Mesh_Algorithms.hh"

// Operators
#include "UniqueLocalIndex.hh"
#include "Upwind.hh"

namespace Amanzi {
namespace Operators {

class UpwindFluxManifolds : public Upwind {
 public:
  UpwindFluxManifolds(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : Upwind(mesh){};
  ~UpwindFluxManifolds(){};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void
  Compute(const CompositeVector& flux, const std::vector<int>& bc_model, CompositeVector& field);

 private:
  int method_;
  double tolerance_;
};

} // namespace Operators
} // namespace Amanzi

#endif
