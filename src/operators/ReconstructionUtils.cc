/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include "ReconstructionUtils.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Cells sharing a node with cell c.
****************************************************************** */
void
getCellAllAdjCells(const AmanziMesh::Mesh& mesh,
                   AmanziMesh::Entity_ID c,
                   std::set<AmanziMesh::Entity_ID>& cells)
{
  cells.clear();

  auto nodes = mesh.getCellNodes(c);
  for (int i = 0; i < nodes.size(); i++) {
    auto vcells = mesh.getNodeCells(nodes[i]);

    for (int k = 0; k < vcells.size(); ++k) {
      int c1 = vcells[k];
      if (c != c1) cells.insert(c1);
    }
  }
}

} // namespace Operators
} // namespace Amanzi
