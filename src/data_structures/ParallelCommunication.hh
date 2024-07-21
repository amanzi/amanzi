/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Data structures

  Supports parallel communications for integer arrays. Eventually this
  functionality should be absorbed in CompositeVector. But currently, it
  supports only double.
*/

#ifndef AMANZI_PARALLEL_COMMUNICATION_HH_
#define AMANZI_PARALLEL_COMMUNICATION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_IntVector.h"

#include "Mesh.hh"

namespace Amanzi {

class ParallelCommunication {
 public:
  ParallelCommunication(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh){};
  ~ParallelCommunication(){};

  // data management
  void CopyMasterEntity2GhostEntity(const AmanziMesh::Entity_kind& kind, Epetra_IntVector& vhost);
  void CombineGhostEntity2MasterEntity(const AmanziMesh::Entity_kind& kind,
                                       Epetra_IntVector& vghost,
                                       Epetra_CombineMode mode);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  std::map<AmanziMesh::Entity_kind, Teuchos::RCP<Epetra_Import>> importer_;
  std::map<AmanziMesh::Entity_kind, bool> importer_initialized_;
};

} // namespace Amanzi

#endif
