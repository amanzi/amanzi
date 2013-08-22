/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------


License: see $AMANZI_DIR/COPYRIGHT
Author: Ethan Coon

A MeshPartition is a collection of non-overlapping Regions which cover a mesh.

------------------------------------------------------------------------- */

#ifndef AMANZI_MESH_FUNCTIONS_PARTITION_
#define AMANZI_MESH_FUNCTIONS_PARTITION_

#include "Teuchos_RCP.hpp"
#include "Epetra_IntVector.h"
#include "MeshDefs.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace Functions {

class MeshPartition {
 public:
  MeshPartition(AmanziMesh::Entity_kind kind,
                const std::vector<std::string>& regions);

  bool initialized() { return initialized_; }
  void Initialize(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  int operator[](AmanziMesh::Entity_ID id) { return (*map_)[id]; }

 protected:
  AmanziMesh::Entity_kind kind_;
  std::vector<std::string> regions_;
  bool initialized_;
  Teuchos::RCP<Epetra_IntVector> map_;

};

} // namespace
} // namespace

#endif
