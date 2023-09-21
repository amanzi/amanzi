/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Helper routines for fracture insertion.
*/

#ifndef AMANZI_FRACTURE_INSERTION_HELPER_HH_
#define AMANZI_FRACTURE_INSERTION_HELPER_HH_

#include "PDE_CouplingFlux.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {

// Populate advective coupling fluxes
void UpdateEnthalpyCouplingFluxes(
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_matrix,
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_fracture,
  const std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux>>& adv_coupling);

} // namespace Amanzi

#endif
