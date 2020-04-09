/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_WHETSTONE_TYPEDEFS_HH_
#define AMANZI_WHETSTONE_TYPEDEFS_HH_

#include <vector>

#include "MeshDefs.hh"

// This variable allows usage of WhetStone as a standalone library.
#define AMANZI_CODE

namespace Amanzi {
namespace WhetStone {

#ifdef AMANZI_CODE
typedef AmanziMesh::Entity_ID Entity_ID;
typedef std::vector<Entity_ID> Entity_ID_List;
typedef AmanziMesh::Parallel_type Parallel_type;
typedef AmanziMesh::Entity_kind Entity_kind;

const int NODE = AmanziMesh::NODE;
const int EDGE = AmanziMesh::EDGE;
const int FACE = AmanziMesh::FACE;
const int CELL = AmanziMesh::CELL;
const int BOUNDARY_FACE = AmanziMesh::BOUNDARY_FACE;

#else
typedef long long int Entity_ID;
typedef std::vector<Entity_ID> Entity_ID_List;

enum Entity_kind { NODE = 0, EDGE, FACE, CELL, BOUNDARY_FACE };

enum class Parallel_type {
  OWNED = 1; // Owned by this processor
  GHOST = 2; // Owned by another processor
  ALL = 3;   // OWNED + GHOST
};
#endif

} // namespace WhetStone
} // namespace Amanzi

#endif
