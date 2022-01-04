/*
  Data structures

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  Supports parallel communications for integer arrays. Eventually this 
  functionality should be absorbed in CompositeVector. But currently, 
  CV supports only double.
*/

#include "ParallelCommunication.hh"

namespace Amanzi {

/* *******************************************************************
* Copy cell-based data from master to ghost positions.
* WARNING: vector vghost must contain ghost cells.
******************************************************************* */
void ParallelCommunication::CopyMasterCell2GhostCell(Epetra_IntVector& vghost)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  const Epetra_BlockMap& target_cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);

  if (!importer_cell_initialized_) {
    importer_cell_ = Teuchos::rcp(new Epetra_Import(target_cmap, source_cmap));
    importer_cell_initialized_ = true;
  }

  int* vdata;
  vghost.ExtractView(&vdata);
  Epetra_IntVector vv(View, source_cmap, vdata);

  vghost.Import(vv, *importer_cell_, Insert);
#endif
}


/* *******************************************************************
* Copy face-based data from master to ghost positions.              
* WARNING: vector vghost must contain ghost cells.              
******************************************************************* */
void ParallelCommunication::CopyMasterCell2GhostCell(const Epetra_IntVector& v, Epetra_IntVector& vghost)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  const Epetra_BlockMap& target_cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);

  if (!importer_cell_initialized_) {
    importer_cell_ = Teuchos::rcp(new Epetra_Import(target_cmap, source_cmap));
    importer_cell_initialized_ = true;
  }

  vghost.Import(v, *importer_cell_, Insert);
#else
  vghost = v;
#endif
}


/* *******************************************************************
* Copy face-based data from master to ghost positions.              
* WARNING: vector vghost must contain ghost cells.              
******************************************************************* */
void ParallelCommunication::CopyMasterFace2GhostFace(Epetra_IntVector& vghost)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE, false);
  const Epetra_BlockMap& target_fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);

  if (!importer_face_initialized_) {
    importer_face_ = Teuchos::rcp(new Epetra_Import(target_fmap, source_fmap));
    importer_face_initialized_ = true;
  }

  int* vdata;
  vghost.ExtractView(&vdata);
  Epetra_IntVector vv(View, source_fmap, vdata);

  vghost.Import(vv, *importer_face_, Insert);
#endif
}


/* *******************************************************************
* Copy face-based data from master to ghost positions.              
* WARNING: vector vghost must contain ghost cells.              
******************************************************************* */
void ParallelCommunication::CopyMasterFace2GhostFace(const Epetra_IntVector& v, Epetra_IntVector& vghost)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE, false);
  const Epetra_BlockMap& target_fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);

  if (!importer_face_initialized_) {
    importer_face_ = Teuchos::rcp(new Epetra_Import(target_fmap, source_fmap));
    importer_face_initialized_ = true;
  }

  vghost.Import(v, *importer_face_, Insert);
#else
  vghost = v;
#endif
}


/* *******************************************************************
* Transfers face-based data from ghost to master positions and 
* performs the operation 'mode' there. 
* WARNING: Vector vghost must contain ghost faces.              
******************************************************************* */
void ParallelCommunication::CombineGhostFace2MasterFace(Epetra_IntVector& vghost, Epetra_CombineMode mode)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE, false);
  const Epetra_BlockMap& target_fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);

  if (!importer_face_initialized_) {
    importer_face_ = Teuchos::rcp(new Epetra_Import(target_fmap, source_fmap));
    importer_face_initialized_ = true;
  }

  int* vdata;
  vghost.ExtractView(&vdata);
  Epetra_IntVector vv(View, source_fmap, vdata);

  vv.Export(vghost, *importer_face_, mode);
#endif
}


/* *******************************************************************
* Transfers cell-based data from ghost to master positions and 
* performs the operation 'mode' there. 
* WARNING: Vector vghost must contain ghost cells.              
******************************************************************* */
void ParallelCommunication::CombineGhostCell2MasterCell(Epetra_IntVector& vghost, Epetra_CombineMode mode)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  const Epetra_BlockMap& target_cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);

  if (!importer_face_initialized_) {
    importer_cell_ = Teuchos::rcp(new Epetra_Import(target_cmap, source_cmap));
    importer_cell_initialized_ = true;
  }

  int* vdata;
  vghost.ExtractView(&vdata);
  Epetra_IntVector vv(View, source_cmap, vdata);

  vv.Export(vghost, *importer_cell_, mode);
#endif
}

}  // namespace Amanzi
