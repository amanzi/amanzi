/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Process Kernels

  Miscalleneous collection of simple non-member functions.
*/

#ifndef AMANZI_PK_UTILS_HH_
#define AMANZI_PK_UTILS_HH_

#include <map>
#include <string>
#include <vector>

#include "CompositeVector.hh"
#include "State.hh"
#include "Tag.hh"
#include "VerboseObject.hh"

namespace Amanzi {

// Average permeability tensor in horizontal direction.
void
PKUtils_CalculatePermeabilityFactorInWell(const Teuchos::Ptr<State>& S,
                                          Teuchos::RCP<Epetra_MultiVector>& Kxy);

AmanziGeometry::Point
PKUtils_EntityCoordinates(int id, AmanziMesh::Entity_ID kind, const AmanziMesh::Mesh& mesh);

void
PKUtils_FluxToVector(const State& S, const CompositeVector& flux, CompositeVector& grad);

} // namespace Amanzi

#endif
