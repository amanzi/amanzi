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
  face-based flux (e.g. Darcy flux).
*/

#ifndef AMANZI_UPWIND_FLUX_HH_
#define AMANZI_UPWIND_FLUX_HH_

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

class UpwindFlux : public Upwind {
 public:
  UpwindFlux(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : Upwind(mesh){};
  ~UpwindFlux(){};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void Compute(const CompositeVector& flux,
               const CompositeVector& solution,
               const std::vector<int>& bc_model,
               CompositeVector& field);

 private:
  int method_, order_;
  double tolerance_;
};

} // namespace Operators
} // namespace Amanzi

#endif
