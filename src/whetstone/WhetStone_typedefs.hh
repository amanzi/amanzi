/*
  WhetStone, version 2.1
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

#include "GeometryDefs.hh"
#include "MeshDefs.hh"

// This variable allows usage of WhetStone as a standalone library.
#define AMANZI_CODE

namespace Amanzi {
namespace WhetStone {

#ifdef AMANZI_CODE
typedef AmanziGeometry::Entity_ID Entity_ID;
typedef std::vector<Entity_ID> Entity_ID_List;
typedef AmanziMesh::Parallel_type ParallelTypeCast;

const int OWNED = AmanziMesh::USED;   // Owned by this processor
const int GHOST = AmanziMesh::GHOST;  // Owned by another processor
const int USED  = AmanziMesh::USED;   // OWNED + GHOST

typedef AmanziMesh::Entity_kind Entity_kind;

#else
typedef long long int Entity_ID;
typedef std::vector<Entity_ID> Entity_ID_List;

enum ParallelTypeCast
{
  PTYPE_UNKNOWN = 0, // Initializer
  OWNED = 1,         // Owned by this processor
  GHOST = 2,         // Owned by another processor
  USED  = 3          // OWNED + GHOST
};

enum Entity_kind
{
  NODE = 0,
  EDGE,
  FACE,
  CELL,
  BOUNDARY_FACE
};
#endif

class Polynomial;
typedef std::vector<std::vector<Polynomial> > MatrixPolynomial;

}  // namespace WhetStone
}  // namespace Amanzi

#endif

