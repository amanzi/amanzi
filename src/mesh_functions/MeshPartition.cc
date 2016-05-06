/*
  Mesh Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A MeshPartition is a collection of non-overlapping regions which cover 
  (optionally) a mesh.
*/

#include "errors.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace Functions {

/* ******************************************************************
* Simple constructor.
****************************************************************** */
MeshPartition::MeshPartition(AmanziMesh::Entity_kind kind,
                             const std::vector<std::string>& regions) :
    kind_(kind),
    regions_(regions),
    initialized_(false) 
{
}


/* ******************************************************************
* Populate the map entity -> number in the list of names.
****************************************************************** */
void MeshPartition::Initialize(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                               const int default_value) 
{
  default_value_ = default_value;
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
      Errors::Message msg("Invalid mesh Entity_kind in MeshPartition");
      Exceptions::amanzi_throw(msg);
  }

  map_ = Teuchos::rcp(new Epetra_IntVector(*mmap_ghost, false));

  // Initialize the data
  map_->PutValue(default_value);

  for (int lcv = 0; lcv != regions_.size(); ++lcv) {
    AmanziMesh::Entity_ID_List block;
    mesh->get_set_entities(regions_[lcv], kind_, AmanziMesh::OWNED, &block);

    for (AmanziMesh::Entity_ID_List::iterator id = block.begin(); id != block.end(); ++id) {
      // Check regions are non-overlapping
      if ((*map_)[*id] >= 0) {
        Errors::Message msg("MeshPartition regions are overlapping");
        Exceptions::amanzi_throw(msg);
      }
      (*map_)[*id] = lcv;
    }
  }

#ifdef HAVE_MPI
  // Scatter to ghost cells
  Epetra_Import importer(*mmap_ghost, *mmap);

  int* vdata;
  map_->ExtractView(&vdata);
  Epetra_IntVector vv(View, *mmap, vdata);

  map_->Import(vv, importer, Insert);
#endif

  initialized_ = true;
}


/* ******************************************************************
* Populate the map entity -> number in the list of names.
****************************************************************** */
void MeshPartition::Initialize(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                               AmanziMesh::Entity_kind kind,
                               const std::vector<std::vector<std::string> >& regions,
                               const int default_value) 
{
  kind_ = kind;
  default_value_ = default_value;
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
      Errors::Message msg("Invalid mesh Entity_kind in MeshPartition");
      Exceptions::amanzi_throw(msg);
  }

  map_ = Teuchos::rcp(new Epetra_IntVector(*mmap_ghost, false));

  // Initialize the data
  map_->PutValue(default_value);

  for (int lcv = 0; lcv != regions.size(); ++lcv) {
    const std::vector<std::string>& regs = regions[lcv]; 
    for (int r = 0; r < regs.size(); ++r) {
      AmanziMesh::Entity_ID_List block;
      mesh->get_set_entities(regs[r], kind_, AmanziMesh::OWNED, &block);

      for (AmanziMesh::Entity_ID_List::iterator id = block.begin(); id != block.end(); ++id) {
        if ((*map_)[*id] >= 0) {
          Errors::Message msg("MeshPartition regions are overlapping");
          Exceptions::amanzi_throw(msg);
        }
        (*map_)[*id] = lcv;
      }
    }
  }

#ifdef HAVE_MPI
  // Scatter to ghost cells
  Epetra_Import importer(*mmap_ghost, *mmap);

  int* vdata;
  map_->ExtractView(&vdata);
  Epetra_IntVector vv(View, *mmap, vdata);

  map_->Import(vv, importer, Insert);
#endif

  initialized_ = true;
}


/* ******************************************************************
* In general, we allow incomplete coverage of the mesh. This 
* routine verifies that ther are no wholes in the map.
****************************************************************** */
void MeshPartition::Verify() const
{
  if (!initialized_) {
    Errors::Message msg("MeshPartition was not initialzied.");
    Exceptions::amanzi_throw(msg);
  }

  for (AmanziMesh::Entity_ID id = 0; id < map_->MyLength(); id++) {
    if ((*map_)[id] == default_value_) {
      Errors::Message msg("MeshPartition regions do not cover the mesh.");
      Exceptions::amanzi_throw(msg);
    }
  }
}

}  // namespace Functions
}  // namespace Amanzi
