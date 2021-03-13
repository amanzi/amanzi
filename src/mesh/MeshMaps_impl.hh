/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Collects maps and importers for Mesh objects.

#pragma once

#include "AmanziTypes.hh"

namespace Amanzi {
namespace AmanziMesh {

// returns used, owned
template<class Mesh>
std::pair<Map, Epetra_Map>
createMapsFromMeshGIDs(const Mesh& mesh, const Entity_kind kind)
{
  Entity_ID num_owned = mesh.getNumEntities(kind, Parallel_type::OWNED);
  Entity_GID_List gids = mesh.getEntityGIDs(kind, Parallel_type::ALL);
  return std::make_pair(
    Teuchos::rcp(new Epetra_Map(-1, gids.size(), gids.data(), 0, *mesh.get_comm())),
    Teuchos::rcp(new Epetra_Map(-1, num_owned, gids.data(), 0, *mesh.get_comm())));
}

// returns used, owned
template<class Mesh>
std::pair<Map, Epetra_Map>
createMapsFromNaturalGIDs(const Mesh& mesh, const Entity_kind kind)
{
  Entity_ID num_owned = mesh.getNumEntities(kind, Parallel_type::OWNED);
  Entity_ID num_all = mesh.getNumEntities(kind, Parallel_type::ALL);
  return std::make_pair(
    Teuchos::rcp(new Epetra_Map(-1, num_all, 0, *mesh.get_comm())),
    Teuchos::rcp(new Epetra_Map(-1, num_owned, 0, *mesh.get_comm())));
}


template<class Mesh>
MeshMap::MeshMap(const Mesh& mesh, bool natural=false)
{
  std::vector<Entity_kind> to_construct{Entity_kind::CELL,
    Entity_kind::FACE, Entity_kind::NODE};
  if (mesh.has_edges()) to_construct.push_back(Entity_kind::EDGE);

  for (const auto& kind : to_construct) {
    std::pair<Teuchos::RCP<Epetra_Map>, Teuchos::RCP<Epetra_Map>> maps;
    if (natural) {
      maps = createMapsFromNaturalGIDs(mesh, kind);
    } else {
      maps = createMapsFromMeshGIDs(mesh, kind);
    }
    all_[kind] = maps.first;
    owned_[kind] = maps.second;
    importer_[kind] = Teuchos::rcp(new Epetra_Import(*maps.first, *maps.second));
  }

  // boundary faces
  // -- get a list of all face LIDs with 1 cell
  std::size_t nfaces_owned = mesh.getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  std::size_t nfaces_all = mesh.getNumEntities(Entity_kind::FACE, Parallel_type::ALL);
  Entity_ID_List bfaces(nfaces_all, -1);
  std::size_t nbf_all = 0;
  std::size_t nbf_owned = 0;
  for (Entity_ID f=0; f!=nfaces_all; ++f) {
    Entity_ID_List fcells;
    mesh.getFaceCells(f, Parallel_type::USED, fcells);
    if (fcells.size() == 1) {
      bfaces[nbf_all] = f;
      ++nbf_all;
      if (f < nfaces_owned) nbf_owned = nbf_all;
    }
  }

  // -- convert to GID
  const auto& fmap_all = *map(Entity_kind::FACE, true);
  for (std::size_t i=0; i!=nbf_all; ++i) {
    bfaces[i] = fmap_all.GID(bfaces[i]);
  }

  // -- construct map, importer
  all_[Entity_kind::BOUNDARY_FACE] =
    Teuchos::rcp(new Epetra_Map(-1, nbf_all, bfaces.data(), 0, *mesh.get_comm()));
  owned_[Entity_kind::BOUNDARY_FACE] =
    Teuchos::rcp(new Epetra_Map(-1, nbf_owned, bfaces.data(), 0, *mesh.get_comm()));
  importer_[Entity_kind::BOUNDARY_FACE] =
    Teuchos::rcp(new Epetra_Import(*all_[Entity_kind::BOUNDARY_FACE], *owned_[Entity_kind::BOUNDARY_FACE]));
  // -- additional importer from face --> boundary_face
  ext_face_importer_ = Teuchos::rcp(new Epetra_Import(*owned_[Entity_kind::BOUNDARY_FACE],
          *owned_[Entity_kind::FACE]));

  // boundary nodes
  // -- form an ordered set of all boundary nodes whose faces are boundary faces
  std::set<Entity_ID> bnodes_lid;
  for (std::size_t i=0; i!=nbf_all; ++i) {
    Entity_ID_List fnodes;
    mesh.getFaceNodes(f, fnodes);
    for (const auto& n : fnodes) bnodes_lid.insert(n);
  }

  // -- convert to GID
  Entity_ID_List bnodes(bnodes_lid.size());
  std::size_t nnodes_owned = mesh.getNumEntities(Entity_kind::NODE, Parallel_type::ALL);
  std::size_t nbn_owned = 0;
  std::size_t nbn_all = 0;
  const auto& nmap_all = *map(Entity_kind::NODE, true);
  for (const auto& bn : bnodes_lid) {
    bnodes[i] = nmap_all.GID(bn);
    ++nbn_all;
    if (bnodes[i] < nnodes_owned) nbn_owned = nbn_all;
  }

  // -- construct map, importer
  all_[Entity_kind::BOUNDARY_NODE] =
    Teuchos::rcp(new Epetra_Map(-1, nbn_all, bnodes.data(), 0, *mesh.get_comm()));
  owned_[Entity_kind::BOUNDARY_NODE] =
    Teuchos::rcp(new Epetra_Map(-1, nbn_owned, bnodes.data(), 0, *mesh.get_comm()));
  importer_[Entity_kind::BOUNDARY_NODE] =
    Teuchos::rcp(new Epetra_Import(*all_[Entity_kind::BOUNDARY_NODE], *owned_[Entity_kind::BOUNDARY_NODE]));
  // -- additional importer from node --> boundary_node
  ext_node_importer_ = Teuchos::rcp(new Epetra_Import(*owned_[Entity_kind::BOUNDARY_NODE],
          *owned_[Entity_kind::NODE]));
}

template<class Mesh>
const Epetra_Map&
map(Entity_kind kind, bool include_ghost) const
{
  if (include_ghost) {
    return *all_.at(kind);
  } else {
    return *owned_.at(kind);
  }
}

template<class Mesh>
const Epetra_Import&
importer(Entity_kind kind) const
{
  return *importer_.at(kind);
}


} // namespace AmanziMesh
} // namespace Amanzi
