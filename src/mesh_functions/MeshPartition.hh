/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

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
  MeshPartition() : initialized_(false){};

  MeshPartition(AmanziMesh::Entity_kind kind,
                const std::vector<std::string>& regions);

  // this routine could be used to optimize mesh coloring
  void Initialize(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                  AmanziMesh::Entity_kind kind,
                  const std::vector<std::vector<std::string>>& regions,
                  const int default_value);

  void Initialize(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                  const int default_value);

  void Verify() const;

  // access
  int operator[](AmanziMesh::Entity_ID id) const { return (*map_)[id]; }
  bool initialized() const { return initialized_; }
  const std::vector<std::string>& regions() const { return regions_; }

 protected:
  AmanziMesh::Entity_kind kind_;
  std::vector<std::string> regions_;
  bool initialized_;
  int default_value_;
  Teuchos::RCP<Epetra_IntVector> map_;
};

} // namespace Functions
} // namespace Amanzi

#endif
