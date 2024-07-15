/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Kontantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Helper functions for unique numbering of entries in mesh lists.

  Some mesh funtions return list of entities that are have no specific order.
  Helper function calculate unique position of entities using their global IDs.

*/

#ifndef UNIQUE_LOCAL_INDEX_HH_
#define UNIQUE_LOCAL_INDEX_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

namespace Amanzi {
namespace Operators {

// unique local index in array returned by getXxxxCells
int
UniqueIndexFaceToCells(const AmanziMesh::Mesh& mesh, int f, int c);
int
UniqueIndexNodeToCells(const AmanziMesh::Mesh& mesh,
                       const AmanziMesh::Mesh& fracture,
                       int v, int v0, int c);

} // namespace Operators
} // namespace Amanzi


#endif
