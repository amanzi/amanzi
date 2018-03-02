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
typedef AmanziMesh::Parallel_type Parallel_type;

#else
typedef long long int Entity_ID;
typedef std::vector<Entity_ID> Entity_ID_List;
typedef int ParallelTypeCast;

enum class Parallel_type {
  OWNED = 1;  // Owned by this processor
  GHOST = 2;  // Owned by another processor
  ALL = 3;    // OWNED + GHOST
};
#endif

class Polynomial;
typedef std::vector<std::vector<Polynomial> > MatrixPolynomial;

}  // namespace WhetStone
}  // namespace Amanzi

#endif

