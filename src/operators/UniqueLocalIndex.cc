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

#include <set>

#include "Epetra_BlockMap.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "UniqueLocalIndex.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Local indix of cell c for a common internal face f
****************************************************************** */
int
UniqueIndexFaceToCells(const AmanziMesh::Mesh& mesh, int f, int c)
{
  int pos = 0;
  auto cells = mesh.getFaceCells(f);
  int ncells = cells.size();
  if (ncells > 1) {
    std::set<int> gids;
    const Epetra_BlockMap& cmap = mesh.getMap(AmanziMesh::Entity_kind::CELL, true);

    for (int i = 0; i < ncells; ++i) gids.insert(cmap.GID(cells[i]));

    auto it = std::find(gids.begin(), gids.end(), cmap.GID(c));
    pos = (it != gids.end()) ? std::distance(gids.begin(), it) : -1;
  }

  return pos;
}


/* ******************************************************************
* Counting of duplicates of node vm in a patch of cells around it
****************************************************************** */
std::vector<int>
UniqueIndexNodeToCells(const AmanziMesh::Mesh& mesh, const std::set<int>& faces_int, int vm, int cm)
{
  std::set<int> gids, faces_all;
  const Epetra_BlockMap& cmap = mesh.getMap(AmanziMesh::Entity_kind::CELL, true);

  auto cells = mesh.getNodeCells(vm, AmanziMesh::Parallel_kind::ALL);
  for (auto c : cells) gids.insert(cmap.GID(c));

  auto faces = mesh.getNodeFaces(vm);
  for (auto f : faces) faces_all.insert(f);

  int cg = (cm < 0) ? -1 : cmap.GID(cm);
  int pos(0);
  std::vector<int> position;

  while (!gids.empty()) {
    int c0 = gids.extract(gids.begin()).value();
    position.push_back(pos);
    if (c0 == cg) return position;

    std::set<int> short_list({ c0 });

    while (!short_list.empty()) {
      int c1 = short_list.extract(short_list.begin()).value();
      const auto& new_faces = mesh.getCellFaces(cmap.LID(c1));

      for (int f : new_faces) {
        if (faces_int.find(f) != faces_int.end()) continue;
        if (faces_all.find(f) != faces_all.end()) {
          auto new_cells = mesh.getFaceCells(f);
          if (new_cells.size() > 1) {
            int c2 = cmap.GID(new_cells[0]) + cmap.GID(new_cells[1]) - c1;
            if (gids.find(c2) != gids.end()) {
              position.push_back(pos);
              if (c2 == cg) return position;

              faces_all.erase(f);
              short_list.insert(c2);
              gids.erase(c2);
            }
          }
        }
      }
    }

    pos++;
  }

  return position;
}

} // namespace Operators
} // namespace Amanzi


#endif
