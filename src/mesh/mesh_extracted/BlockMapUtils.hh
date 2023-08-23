/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Mesh

*/

#ifndef AMANZI_BLOCK_MAP_UTILS_HH_
#define AMANZI_BLOCK_MAP_UTILS_HH_

// TPLs
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"

// Amanzi
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

// Convert discontinuous map to continuous map
// -- parent_maps - pair of consituous master and ghosted maps
// -- subset_maps - pair of discontinuous master and ghosted maps
//                 this must be a subset of parent_maps
// Returns pair of continuous master and ghosted maps that corresponds to subset maps
std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map>>
createContinuousMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                     const std::pair<Teuchos::RCP<const Epetra_BlockMap>,
                                     Teuchos::RCP<const Epetra_BlockMap>>& parent_maps,
                     const std::pair<Teuchos::RCP<const Epetra_BlockMap>,
                                     Teuchos::RCP<const Epetra_BlockMap>>& subset_maps);

} // namespace AmanziMesh
} // namespace Amanzi

#endif
