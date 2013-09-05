/* -------------------------------------------------------------------------

License: see $AMANZI_DIR/COPYRIGHT
Author: Ethan Coon

A MeshPartition is a collection of non-overlapping Regions which cover 
(otionally) a mesh.

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

  void Initialize(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                  const int default_value);
  void Verify() const;

  // access
  int operator[](AmanziMesh::Entity_ID id) const { return (*map_)[id]; }
  bool initialized() { return initialized_; }

 protected:
  AmanziMesh::Entity_kind kind_;
  std::vector<std::string> regions_;
  bool initialized_;
  int default_value_;
  Teuchos::RCP<Epetra_IntVector> map_;
};

}  // namespace Functions
}  // namespace Amanzi

#endif
