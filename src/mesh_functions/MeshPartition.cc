/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------


License: see $AMANZI_DIR/COPYRIGHT
Author: Ethan Coon

A MeshPartition is a collection of non-overlapping Regions which cover a mesh.

------------------------------------------------------------------------- */

#include "errors.hh"
#include "MeshPartition.hh"


namespace Amanzi {
namespace Functions {

MeshPartition::MeshPartition(AmanziMesh::Entity_kind kind,
                             const std::vector<std::string>& regions) :
    kind_(kind),
    regions_(regions),
    initialized_(false) {}

void MeshPartition::Initialize(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
  Teuchos::Ptr<const Epetra_BlockMap> mmap, mmap_ghost;

  // Create the data
  switch(kind_) {
    case AmanziMesh::CELL:
      mmap = Teuchos::ptr(&mesh->cell_map(false));
      mmap_ghost = Teuchos::ptr(&mesh->cell_map(true));
      break;
    case AmanziMesh::FACE:
      mmap = Teuchos::ptr(&mesh->face_map(false));
      mmap_ghost = Teuchos::ptr(&mesh->face_map(true));
      break;
    case AmanziMesh::NODE:
      mmap = Teuchos::ptr(&mesh->node_map(false));
      mmap_ghost = Teuchos::ptr(&mesh->node_map(true));
      break;
    default:
      Errors::Message message("Invalid mesh Entity_kind in MeshPartition");
      Exceptions::amanzi_throw(message);
  }

  map_ = Teuchos::rcp(new Epetra_IntVector(*mmap_ghost, false));

  // Initialize the data
  map_->PutValue(-1);
  for (int lcv=0; lcv!=regions_.size(); ++lcv) {
    AmanziMesh::Entity_ID_List block;
    mesh->get_set_entities(regions_[lcv], kind_, AmanziMesh::OWNED, &block);

    for (AmanziMesh::Entity_ID_List::iterator id = block.begin();
         id != block.end(); ++id) {
      // Check regions are non-overlapping
      if ((*map_)[*id] >= 0) {
        Errors::Message message("MeshPartition regions are overlapping");
        Exceptions::amanzi_throw(message);
      }
      (*map_)[*id] = lcv;
    }
  }

  // Scatter to ghost cells
#ifdef HAVE_MPI
  Epetra_Import importer(*mmap_ghost, *mmap);

  int* vdata;
  map_->ExtractView(&vdata);
  Epetra_IntVector vv(View, *mmap, vdata);

  map_->Import(vv, importer, Insert);
#endif

  // Ensure regions are a cover
  if (map_->MinValue() < 0) {
    Errors::Message message("MeshPartition regions do not cover the mesh.");
    Exceptions::amanzi_throw(message);
  }

  initialized_ = true;
}


} // namespace
} // namespace
