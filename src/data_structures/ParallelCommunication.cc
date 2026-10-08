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
  functionality should be absorbed in CompositeVector. But currently,
  CV supports only double.
*/

#include "ParallelCommunication.hh"

namespace Amanzi {

/* *******************************************************************
* Copy data from master to ghost positions.
* WARNING: vector vghost must contain ghosts.
******************************************************************* */
void
ParallelCommunication::CopyMasterEntity2GhostEntity(const AmanziMesh::Entity_kind& kind,
                                                    Epetra_IntVector& vghost)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_map = mesh_->getMap(kind, false);
  const Epetra_BlockMap& target_map = mesh_->getMap(kind, true);

  if (importer_initialized_.find(kind) == importer_initialized_.end()) {
    importer_[kind] = Teuchos::rcp(new Epetra_Import(target_map, source_map));
    importer_initialized_[kind] = true;
  }

  int* vdata;
  vghost.ExtractView(&vdata);
  Epetra_IntVector vv(View, source_map, vdata);

  vghost.Import(vv, *importer_[kind], Insert);
#endif
}


/* *******************************************************************
* Transfer data from ghost to master positions and perform operation 'mode'.
* WARNING: Vector vghost must contain ghost faces.
******************************************************************* */
void
ParallelCommunication::CombineGhostEntity2MasterEntity(const AmanziMesh::Entity_kind& kind,
                                                       Epetra_IntVector& vghost,
                                                       Epetra_CombineMode mode)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_map = mesh_->getMap(kind, false);
  const Epetra_BlockMap& target_map = mesh_->getMap(kind, true);

  if (importer_initialized_.find(kind) == importer_initialized_.end()) {
    importer_[kind] = Teuchos::rcp(new Epetra_Import(target_map, source_map));
    importer_initialized_[kind] = true;
  }

  int* vdata;
  vghost.ExtractView(&vdata);
  Epetra_IntVector vv(View, source_map, vdata);

  vv.Export(vghost, *importer_[kind], mode);
#endif
}

} // namespace Amanzi
