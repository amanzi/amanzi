/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Collects maps and importers for Mesh objects.

#include <set>

#include "AmanziTypes.hh"
#include "AmanziMap.hh"
#include "AmanziVector.hh"
#include "MeshDefs.hh"
#include "MeshAlgorithms.hh"
#include "MeshMaps.hh"

namespace Amanzi {
namespace AmanziMesh {


//
// Create an owned, ghosted pair of maps, using a provided pair to communicate
// the new ghost entities.
//
std::pair<Map_ptr_type, Map_ptr_type>
createContiguousMaps(const Map_ptr_type& ghosted, const Map_ptr_type& owned)
{
  // create the owned contiguous map.
  auto owned_contiguous = Teuchos::rcp(
    new Map_type(owned->getGlobalNumElements(), owned->getLocalNumElements(), 0, owned->getComm()));

  // communicated owned to ghosted using the mesh's maps
  Vector_type_<GO> owned_contiguous_vec(owned);
  {
    auto data = owned_contiguous_vec.getLocalViewHost(Tpetra::Access::ReadWrite);
    for (int i = 0; i != data.size(); ++i) { data(i, 0) = owned_contiguous->getGlobalElement(i); }
  }
  Import_type importer(owned, ghosted);
  Vector_type_<GO> all_contiguous_vec(ghosted);
  all_contiguous_vec.doImport(owned_contiguous_vec, importer, Tpetra::CombineMode::INSERT);

  auto ghosted_contiguous =
    Teuchos::rcp(new Map_type(-1, all_contiguous_vec.getData()(), 0, ghosted->getComm()));
  return std::make_pair(ghosted_contiguous, owned_contiguous);
}


void
MeshMaps::initialize(const MeshFramework& mesh, bool renumber)
{
  std::vector<Entity_kind> to_construct{ Entity_kind::CELL, Entity_kind::FACE };
  if (mesh.hasEdges()) to_construct.push_back(Entity_kind::EDGE);
  if (mesh.hasNodes()) to_construct.push_back(Entity_kind::NODE);

  for (const auto& kind : to_construct) {
    std::pair<Map_ptr_type, Map_ptr_type> maps;
    if (renumber) {
      maps = createMapsFromContiguousGIDs(mesh, kind);
    } else {
      maps = createMapsFromMeshGIDs(mesh, kind);
    }
    all_[kind] = maps.first;
    owned_[kind] = maps.second;
    importer_[kind] =
      Teuchos::rcp(new Import_type(maps.second, maps.first)); // NOTE: Tpetra swapped order!
  }

  // boundary faces
  // -- get a list of all face LIDs with 1 cell
  std::size_t nfaces_owned = mesh.getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  std::size_t nfaces_all = mesh.getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);
  MeshFramework::Entity_ID_View boundary_faces("boundary_faces", nfaces_all);
  Kokkos::deep_copy(boundary_faces, -1);

  std::size_t nbf_owned = 0;
  std::size_t nbf_all = 0;
  // This used to loop over nfaces_all but logically that includes internal ghosted faces
  // that have one cell on this rank but that more than one cell on the other rank.
  //
  // On the other hand, if a boundary face is owned and only has one cell, then
  // it is a true boundary face.
  for (Entity_ID f = 0; f != nfaces_all; ++f) {
    MeshFramework::cEntity_ID_View fcells;
    mesh.getFaceCells(f, fcells);
    if (fcells.size() == 1) {
      boundary_faces[nbf_all++] = f;
      if (f < nfaces_owned) nbf_owned = nbf_all;
    }
  }
  Kokkos::resize(boundary_faces, nbf_all);

  // -- convert to GID
  MeshFramework::Entity_GID_View boundary_face_GIDs(
    "boundary_face_GIDs", boundary_faces.size()); // GID of the face corresponding to the bf
  MeshFramework::Entity_GID_View boundary_face_internal_cell_GIDs(
    "boundary_face_internal_cell_GIDs",
    boundary_faces.size()); // GID of the cell internal to the bf
  const auto& fmap_all = getMap(Entity_kind::FACE, true);
  const auto& cmap_all = getMap(Entity_kind::CELL, true);
  for (std::size_t i = 0; i != nbf_all; ++i) {
    boundary_face_GIDs[i] = fmap_all->getGlobalElement(boundary_faces[i]);
    MeshFramework::cEntity_ID_View fcells;
    mesh.getFaceCells(boundary_faces[i], fcells);
    boundary_face_internal_cell_GIDs[i] = cmap_all->getGlobalElement(fcells[0]);
  }

  // -- construct owned maps
  auto boundary_face_GIDs_owned =
    Kokkos::subview(boundary_face_GIDs, Kokkos::make_pair((std::size_t)0, nbf_owned));
  auto boundary_face_internal_cell_GIDs_owned =
    Kokkos::subview(boundary_face_internal_cell_GIDs, Kokkos::make_pair((std::size_t)0, nbf_owned));
  owned_[Entity_kind::BOUNDARY_FACE] =
    Teuchos::rcp(new Map_type(-1, boundary_face_GIDs_owned, 0, mesh.getComm()));
  auto bf_internal_cells_map =
    Teuchos::rcp(new Map_type(-1, boundary_face_internal_cell_GIDs_owned, 0, mesh.getComm()));

  // Get the local IDs (lc_id) of copies of owned boundary faces on remote
  // processors (pr_id).  We must check if a ghost face that has only 1 cell is
  // an exterior face, or on a process boundary.  This is done by seeing if it
  // appears in the owned boundary face list on another processor (pr_id >= 0)
  // -- if it does, it is a true exterior face.  If it does not, then it is
  // interior on the other process, so is a processor boundary.
  int nbf_notowned = nbf_all - nbf_owned;

  Teuchos::Array<int> pr_id(nbf_notowned), lc_id(nbf_notowned);
  Teuchos::ArrayView<int> bf_gids_notowned(&boundary_face_GIDs[nbf_owned], nbf_notowned);
  owned_[Entity_kind::BOUNDARY_FACE]->getRemoteIndexList(bf_gids_notowned, pr_id, lc_id);

  int nbf_notowned_tight = 0;
  for (int i = 0; i != nbf_notowned; ++i) {
    if (pr_id[i] >= 0) {
      boundary_faces[nbf_owned + nbf_notowned_tight] = boundary_faces[nbf_owned + i];
      boundary_face_GIDs[nbf_owned + nbf_notowned_tight] = boundary_face_GIDs[nbf_owned + i];
      ++nbf_notowned_tight;
    }
  }
  nbf_all = nbf_owned + nbf_notowned_tight;
  Kokkos::resize(boundary_faces, nbf_all);
  Kokkos::resize(boundary_face_GIDs, nbf_all);

  // create the dual view
  boundary_faces_ = asDualView(boundary_faces);

  // -- construct the all map
  all_[Entity_kind::BOUNDARY_FACE] =
    Teuchos::rcp(new Map_type(-1, boundary_face_GIDs, 0, mesh.getComm()));
  // NOTE: Tpetra swapped order!
  importer_[Entity_kind::BOUNDARY_FACE] = Teuchos::rcp(
    new Import_type(owned_[Entity_kind::BOUNDARY_FACE], all_[Entity_kind::BOUNDARY_FACE]));
  // -- additional importer from face --> boundary_face
  boundary_face_importer_ =
    Teuchos::rcp(new Import_type(owned_[Entity_kind::FACE], owned_[Entity_kind::BOUNDARY_FACE]));

  // -- importer from internal cell --> boundary_face
  boundary_face_internal_cell_importer_ =
    Teuchos::rcp(new Import_type(owned_[Entity_kind::CELL], bf_internal_cells_map));

  // boundary nodes
  // -- form an ordered set of all boundary nodes whose faces are boundary faces
  if (mesh.hasNodes()) {
    std::set<Entity_ID> bnodes_lid;
    for (std::size_t i = 0; i != nbf_all; ++i) {
      MeshFramework::cEntity_ID_View fnodes;
      mesh.getFaceNodes(boundary_faces[i], fnodes);
      for (const auto& n : fnodes) bnodes_lid.insert(n);
    }

    // -- convert to GID
    MeshFramework::Entity_ID_View boundary_nodes("boundary_nodes", bnodes_lid.size());
    Kokkos::deep_copy(boundary_nodes, -1);
    MeshFramework::Entity_GID_View boundary_node_GIDs("boundary_node_GIDs", bnodes_lid.size());
    std::size_t nnodes_owned = mesh.getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
    std::size_t nbn_owned = 0;
    std::size_t nbn_all = 0;
    const auto& nmap_all = getMap(Entity_kind::NODE, true);
    for (const auto& bn : bnodes_lid) {
      boundary_nodes[nbn_all] = bn;
      boundary_node_GIDs[nbn_all++] = nmap_all->getGlobalElement(bn);
      if (bn < nnodes_owned) nbn_owned = nbn_all;
    }
    Kokkos::resize(boundary_nodes, nbn_all);
    Kokkos::resize(boundary_node_GIDs, nbn_all);

    // create dual view
    boundary_nodes_ = asDualView(boundary_nodes);

    // -- construct map, importer
    all_[Entity_kind::BOUNDARY_NODE] =
      Teuchos::rcp(new Map_type(-1, boundary_node_GIDs, 0, mesh.getComm()));

    auto boundary_node_GIDs_owned =
      Kokkos::subview(boundary_node_GIDs, Kokkos::make_pair((std::size_t)0, nbn_owned));
    owned_[Entity_kind::BOUNDARY_NODE] =
      Teuchos::rcp(new Map_type(-1, boundary_node_GIDs_owned, 0, mesh.getComm()));
    // NOTE: Tpetra swapped order!
    importer_[Entity_kind::BOUNDARY_NODE] = Teuchos::rcp(
      new Import_type(owned_[Entity_kind::BOUNDARY_NODE], all_[Entity_kind::BOUNDARY_NODE]));
    // -- additional importer from node --> boundary_node
    boundary_node_importer_ =
      Teuchos::rcp(new Import_type(owned_[Entity_kind::NODE], owned_[Entity_kind::BOUNDARY_NODE]));
  }
}


const Map_ptr_type&
MeshMaps::getMap(Entity_kind kind, bool include_ghost) const
{
  if (all_.count(kind) == 0) {
    Errors::Message msg;
    msg << "This mesh does not support entity kind " << to_string(kind);
    Exceptions::amanzi_throw(msg);
  }
  if (include_ghost) {
    return all_.at(kind);
  } else {
    return owned_.at(kind);
  }
}

const Import_type&
MeshMaps::getImporter(Entity_kind kind) const
{
  return *importer_.at(kind);
}

std::size_t
MeshMaps::getNBoundaryFaces(Parallel_kind ptype) const
{
  if (Parallel_kind::OWNED == ptype)
    return owned_.at(Entity_kind::BOUNDARY_FACE)->getLocalNumElements();
  else
    return all_.at(Entity_kind::BOUNDARY_FACE)->// getLocalNumElements();
// }

// std::size_t
// MeshMaps::getNBoundaryNodes(Parallel_kind ptype) const
// {
//   if (Parallel_kind::OWNED == ptype)
//     return owned_.at(Entity_kind::BOUNDARY_NODE)->getLocalNumElements();
//   else
//     return all_.at(Entity_kind::BOUNDARY_NODE)->getLocalNumElements();
}


} // namespace AmanziMesh
} // namespace Amanzi
