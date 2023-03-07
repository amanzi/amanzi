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

#include <iterator>
#include <set>

#include "Epetra_BlockMap.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Local index of cells for common internal face
****************************************************************** */
int
UniqueIndexFaceToCells(const AmanziMesh::Mesh& mesh, int f, int c)
{
  int pos = 0;
  AmanziMesh::Entity_ID_List cells;

  mesh.face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  int ncells = cells.size();
  if (ncells > 1) {
    std::set<int> gids;
    const Epetra_BlockMap& cmap = mesh.cell_map(true);

    for (int i = 0; i < ncells; ++i) gids.insert(cmap.GID(cells[i]));

    auto it = std::find(gids.begin(), gids.end(), cmap.GID(c));
    pos = (it != gids.end()) ? std::distance(gids.begin(), it) : -1;
  }

  return pos;
}


/* ******************************************************************
* Order cells by their global ids. Returns 1 if cells were swapped.
****************************************************************** */
int
OrderCellsByGlobalId(const AmanziMesh::Mesh& mesh,
                     const AmanziMesh::Entity_ID_List& cells, int& c1, int& c2)
{
  c1 = cells[0];
  c2 = -1;

  int ncells = cells.size();
  if (ncells == 1) return 0;

  c2 = cells[1];
  if (mesh.cell_map(true).GID(c1) > mesh.cell_map(true).GID(c2)) {
    int c(c1);
    c1 = c2;
    c2 = c;
    return 1;
  }

  return 0;
}


/* ******************************************************************
* Order cells by their global ids.
****************************************************************** */
void
OrderCellsByGlobalId(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID_List& cells)
{
  int ncells = cells.size();

  for (int i = 0; i < ncells; ++i) {
    int i0 = i;
    int gid0 = mesh.cell_map(true).GID(cells[i]);
    for (int j = i + 1; j < ncells; ++j) {
      int gid1 = mesh.cell_map(true).GID(cells[j]);
      if (gid1 < gid0) {
        i0 = j; 
        gid0 = gid1;
      }
    }
    if (i != i0) std::swap(cells[i], cells[i0]);
  }
}

} // namespace Operators
} // namespace Amanzi


#endif
