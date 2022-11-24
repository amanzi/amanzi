/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Helper functions for unique numbering of entries in mesh lists.

  Some mesh funtions return list of entities that are have no specific order.
  Helper function calculate unique position of entities using their global IDs.

  Author: Kontantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef UNIQUE_LOCAL_INDEX_HH_
#define UNIQUE_LOCAL_INDEX_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

namespace Amanzi {
namespace Operators {

// unique index for array returned by face_get_cells
int
UniqueIndexFaceToCells(const AmanziMesh::Mesh& mesh, int f, int c);

} // namespace Operators
} // namespace Amanzi


#endif
