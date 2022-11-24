/*
  Data structures

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

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
  ParallelCommunication(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh)
  {
    importer_cell_initialized_ = false;
    importer_face_initialized_ = false;
  }
  ~ParallelCommunication(){};

  // data management
  void CopyMasterCell2GhostCell(Epetra_IntVector& vhost);
  void CopyMasterCell2GhostCell(const Epetra_IntVector& v, Epetra_IntVector& vhost);
  void CopyMasterFace2GhostFace(Epetra_IntVector& vhost);
  void CopyMasterFace2GhostFace(const Epetra_IntVector& v, Epetra_IntVector& vhost);
  void CombineGhostFace2MasterFace(Epetra_IntVector& v, Epetra_CombineMode mode = Insert);
  void CombineGhostCell2MasterCell(Epetra_IntVector& v, Epetra_CombineMode mode = Insert);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  Teuchos::RCP<Epetra_Import> importer_cell_;
  bool importer_cell_initialized_;

  Teuchos::RCP<Epetra_Import> importer_face_;
  bool importer_face_initialized_;
};

} // namespace Amanzi

#endif
