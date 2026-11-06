/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Supporting tools for building reconstruction patches.
*/

#ifndef AMANZI_RECONSTRUCTION_UTILS_HH_
#define AMANZI_RECONSTRUCTION_UTILS_HH_

#include <set>

#include "Mesh.hh"

namespace Amanzi {
namespace Operators {

// all cells sharing at least one node with cell c, excluding c
void
getCellAllAdjCells(const AmanziMesh::Mesh& mesh,
                   AmanziMesh::Entity_ID c,
                   std::set<AmanziMesh::Entity_ID>& cells);

} // namespace Operators
} // namespace Amanzi

#endif
