/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Collects maps and importers for Mesh objects.

#pragma once

#include <set>
#include "AmanziTypes.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

// returns used, owned
template<class Mesh_type>
std::pair<Map_ptr_type, Map_ptr_type>
createMapsFromMeshGIDs(const Mesh_type& mesh, const Entity_kind kind)
{
  Entity_ID num_owned = mesh.getNumEntities(kind, Parallel_type::OWNED);
  Entity_GID_List gids = mesh.getEntityGIDs(kind, Parallel_type::ALL);
  return std::make_pair(
    Teuchos::rcp(new Epetra_Map(-1, gids.size(), gids.data(), 0, *mesh.getComm())),
    Teuchos::rcp(new Epetra_Map(-1, num_owned, gids.data(), 0, *mesh.getComm())));
}

// returns used, owned
template<class Mesh_type>
std::pair<Map_ptr_type, Map_ptr_type>
createMapsFromContiguousGIDs(const Mesh_type& mesh, const Entity_kind kind)
{
  Entity_ID num_owned = mesh.getNumEntities(kind, Parallel_type::OWNED);
  Entity_ID num_all = mesh.getNumEntities(kind, Parallel_type::ALL);
  return std::make_pair(
    Teuchos::rcp(new Epetra_Map(-1, num_all, 0, *mesh.getComm())),
    Teuchos::rcp(new Epetra_Map(-1, num_owned, 0, *mesh.getComm())));
}


template<class Mesh_type>
void MeshMaps::initialize(const Mesh_type& mesh, bool renumber)
{
  std::vector<Entity_kind> to_construct{Entity_kind::CELL,
    Entity_kind::FACE, Entity_kind::NODE};
  if (mesh.hasEdges()) to_construct.push_back(Entity_kind::EDGE);

  for (const auto& kind : to_construct) {
    std::pair<Teuchos::RCP<Epetra_Map>, Teuchos::RCP<Epetra_Map>> maps;
    if (renumber) {
      maps = createMapsFromContiguousGIDs(mesh, kind);
    } else {
      maps = createMapsFromMeshGIDs(mesh, kind);
    }
    all_[kind] = maps.first;
    owned_[kind] = maps.second;
    importer_[kind] = Teuchos::rcp(new Epetra_Import(*maps.first, *maps.second));
  }

  // boundary faces
  // -- get a list of all owned face LIDs with 1 cell
  std::size_t nfaces_owned = mesh.getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  std::size_t nfaces_all = mesh.getNumEntities(Entity_kind::FACE, Parallel_type::ALL);
  boundary_faces_.resize(nfaces_all, -1);

  std::size_t nbf_owned = 0;
  std::size_t nbf_all = 0;
  for (Entity_ID f=0; f!=nfaces_all; ++f) {
    Entity_ID_List fcells;
    mesh.getFaceCells(f, Parallel_type::ALL, fcells);
    if (fcells.size() == 1) {
      boundary_faces_[nbf_all++] = f;
      if (f < nfaces_owned) nbf_owned = nbf_all;
    }
  }

  // -- convert to GID
  Entity_GID_List boundary_face_GIDs(boundary_faces_.size());
  const auto& fmap_all = getMap(Entity_kind::FACE, true);
  for (std::size_t i=0; i!=nbf_all; ++i) {
    boundary_face_GIDs[i] = fmap_all.GID(boundary_faces_[i]);
  }

  // -- construct owned map
  owned_[Entity_kind::BOUNDARY_FACE] =
    Teuchos::rcp(new Epetra_Map(-1, nbf_owned, boundary_face_GIDs.data(), 0, *mesh.getComm()));

  // Get the local IDs (lc_id) of copies of owned boundary faces on remote
  // processors (pr_id).  We must check if a ghost face that has only 1 cell is
  // an exterior face, or on a process boundary.  This is done by seeing if it
  // appears in the owned boundary face list on another processor (pr_id >= 0)
  // -- if it does, it is a true exterior face.  If it does not, then it is
  // interior on the other process, so is a processor boundary.
  int nbf_notowned = nbf_all - nbf_owned;
  std::vector<int> pr_id(nbf_notowned), lc_id(nbf_notowned);
  owned_[Entity_kind::BOUNDARY_FACE]->RemoteIDList(nbf_notowned, &boundary_face_GIDs[nbf_owned],
          pr_id.data(), lc_id.data());

  int nbf_notowned_tight = 0;
  for (int i=0; i!=nbf_notowned; ++i) {
    if (pr_id[i] >= 0) {
      boundary_faces_[nbf_owned+nbf_notowned_tight] = boundary_faces_[nbf_owned+i];
      boundary_face_GIDs[nbf_owned+nbf_notowned_tight] = boundary_face_GIDs[nbf_owned+i];
      ++nbf_notowned_tight;
    }
  }
  nbf_all = nbf_owned + nbf_notowned_tight;
  boundary_faces_.resize(nbf_all);

  // -- construct the all map
  all_[Entity_kind::BOUNDARY_FACE] =
    Teuchos::rcp(new Epetra_Map(-1, nbf_all, boundary_face_GIDs.data(), 0, *mesh.getComm()));
  importer_[Entity_kind::BOUNDARY_FACE] =
    Teuchos::rcp(new Epetra_Import(*all_[Entity_kind::BOUNDARY_FACE], *owned_[Entity_kind::BOUNDARY_FACE]));
  // -- additional importer from face --> boundary_face
  boundary_face_importer_ = Teuchos::rcp(new Epetra_Import(*owned_[Entity_kind::BOUNDARY_FACE],
          *owned_[Entity_kind::FACE]));

  // boundary nodes
  // -- form an ordered set of all boundary nodes whose faces are boundary faces
  if (mesh.hasNodes()) {
    std::set<Entity_ID> bnodes_lid;
    for (std::size_t i=0; i!=nbf_all; ++i) {
      Entity_ID_List fnodes;
      mesh.getFaceNodes(boundary_faces_[i], fnodes);
      for (const auto& n : fnodes) bnodes_lid.insert(n);
    }

    // -- convert to GID
    boundary_nodes_.resize(bnodes_lid.size());
    Entity_GID_List boundary_node_GIDs(bnodes_lid.size());
    std::size_t nnodes_owned = mesh.getNumEntities(Entity_kind::NODE, Parallel_type::OWNED);
    std::size_t nbn_owned = 0;
    std::size_t nbn_all = 0;
    const auto& nmap_all = getMap(Entity_kind::NODE, true);
    for (const auto& bn : bnodes_lid) {
      boundary_nodes_[nbn_all] = bn;
      boundary_node_GIDs[nbn_all++] = nmap_all.GID(bn);
      if (bn < nnodes_owned) nbn_owned = nbn_all;
    }

    // -- construct map, importer
    all_[Entity_kind::BOUNDARY_NODE] =
      Teuchos::rcp(new Epetra_Map(-1, nbn_all, boundary_node_GIDs.data(), 0, *mesh.getComm()));
    owned_[Entity_kind::BOUNDARY_NODE] =
      Teuchos::rcp(new Epetra_Map(-1, nbn_owned, boundary_node_GIDs.data(), 0, *mesh.getComm()));
    importer_[Entity_kind::BOUNDARY_NODE] =
      Teuchos::rcp(new Epetra_Import(*all_[Entity_kind::BOUNDARY_NODE], *owned_[Entity_kind::BOUNDARY_NODE]));
    // -- additional importer from node --> boundary_node
    boundary_node_importer_ = Teuchos::rcp(new Epetra_Import(*owned_[Entity_kind::BOUNDARY_NODE],
            *owned_[Entity_kind::NODE]));
  }
}

const Map_type&
MeshMaps::getMap(Entity_kind kind, bool include_ghost) const
{
  if (all_.count(kind) == 0) {
    Errors::Message msg;
    msg << "This mesh does not support entity kind " << to_string(kind);
    Exceptions::amanzi_throw(msg);
  }
  if (include_ghost) {
    return *all_.at(kind);
  } else {
    return *owned_.at(kind);
  }
}

const Epetra_Import&
MeshMaps::getImporter(Entity_kind kind) const
{
  return *importer_.at(kind);
}


} // namespace AmanziMesh
} // namespace Amanzi