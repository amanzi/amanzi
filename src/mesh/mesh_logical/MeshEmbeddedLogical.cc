/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "Epetra_IntVector.h"

#include "RegionEnumerated.hh"
#include "MeshLogical.hh"
#include "MeshEmbeddedLogical.hh"

namespace Amanzi {
namespace AmanziMesh {

/*
MeshEmbeddedLogical Constructor

Combines two meshes, each representing their own domain, into a
single mesh via set of specified connections.  The result is a
logical mesh in the sense that it provides a limited interface.

 - face_cell_list          : length nfaces array of length 2 arrays
                             defining the topology
 - face_cell_lengths       : length of the cell-to-face connection
 - face_area_normals       : length nfaces array of normals of the
                             face, points from cell 1 to 2 in
                             face_cell_list topology, magnitude
                             is area

Note that entities are ordered in the following way:

1. logical mesh entities
2. "extra" entities (valid only for faces, these are
   the connections between the two meshes.
3. background entities

*/
MeshEmbeddedLogical::MeshEmbeddedLogical(
  const Comm_ptr_type& comm,
  Teuchos::RCP<MeshFramework> bg_mesh,
  Teuchos::RCP<MeshFramework> log_mesh,
  const std::vector<Entity_ID_List>& face_cell_ids,
  const std::vector<std::vector<double>>& face_cell_lengths,
  const std::vector<AmanziGeometry::Point>& face_area_normals,
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : MeshFramework(comm, Teuchos::null, plist),
    bg_mesh_(bg_mesh),
    log_mesh_(log_mesh),
    extra_face_cell_ids_(face_cell_ids),
    extra_face_cell_lengths_(face_cell_lengths),
    extra_face_area_normals_(face_area_normals)
{
  setSpaceDimension(3);
  setManifoldDimension(3);
  setAlgorithms(Teuchos::rcp(new MeshLogicalAlgorithms()));

  // merge and remap to get new global IDs
  int ncells_bg_owned = bg_mesh->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  int ncells_bg_all = bg_mesh->getNumEntities(Entity_kind::CELL, Parallel_type::ALL);
  int ncells_log = log_mesh->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  int ncells_my_all = ncells_bg_all + ncells_log;

  int nfaces_bg_owned = bg_mesh->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  int nfaces_bg_all = bg_mesh->getNumEntities(Entity_kind::FACE, Parallel_type::ALL);
  int nfaces_log = log_mesh->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  int nfaces_extra = extra_face_cell_ids_.size<MemSpace_type::HOST>();

  // compute "extra" bisectors
  extra_face_cell_bisectors_.resize(nfaces_extra, 2);

  // now loop over new faces, adding the updates to the cell-ordered versions
  for (int fi = 0; fi != nfaces_extra; ++fi) {
    auto normal0 = extra_face_area_normals_[fi];
    extra_face_cell_bisectors_.getRow<MemSpace_type::HOST>(fi)[0] =
      face_cell_lengths[fi][0] / AmanziGeometry::norm(normal0) * normal0;

    auto normal1 = extra_face_area_normals_[fi];
    extra_face_cell_bisectors_.getRow<MemSpace_type::HOST>(fi)[1] =
      face_cell_lengths[fi][1] / AmanziGeometry::norm(normal1) * normal1;
  }
  extra_face_cell_bisectors_.update<MemSpace_type::DEVICE>();

  // Need to renumber the global IDs.  Use the corresponding mesh to communicate.
  // First cells
  // -- create a map of owned CELLs
  int ncells_total_owned = getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  int ncells_total_all = getNumEntities(Entity_kind::CELL, Parallel_type::ALL);
  Epetra_Map cell_owned_map(-1, ncells_total_owned, 0, *getComm());

  // -- create a map on the background mesh of all CELLs
  Entity_GID_List bg_cell_gids = bg_mesh_->getEntityGIDs(Entity_kind::CELL, Parallel_type::ALL);
  Epetra_Map bg_cell_all_map(-1, ncells_bg_all, bg_cell_gids.data(), 0, *bg_mesh_->getComm());
  Epetra_Map bg_cell_owned_map(-1, ncells_bg_owned, bg_cell_gids.data(), 0, *bg_mesh_->getComm());

  // -- create a vector, fill the owned entries, and scatter
  Epetra_IntVector bg_cell_owned_vec(bg_cell_owned_map);
  for (int c = 0; c != ncells_bg_owned; ++c)
    bg_cell_owned_vec[c] = cell_owned_map.GID(c + ncells_log);
  Epetra_IntVector bg_cell_all_vec(bg_cell_all_map);
  Epetra_Import cell_import(bg_cell_all_map, bg_cell_owned_map);
  bg_cell_all_vec.Import(bg_cell_owned_vec, cell_import, Insert);

  // -- copy ghost GIDs into the renumbered vector
  Entity_GID_List cell_gids(ncells_bg_all + ncells_log);
  for (int c = 0; c != ncells_total_owned; ++c) cell_gids[c] = cell_owned_map.GID(c);
  for (int c = ncells_total_owned; c != ncells_total_all; ++c)
    cell_gids[c] = bg_cell_all_vec[c - ncells_log + ncells_bg_owned];

  // -- create the map
  cell_map_ = Teuchos::rcp(new Epetra_Map(-1, cell_gids.size(), cell_gids.data(), 0, *getComm()));

  // Next faces
  // -- create a map of owned FACES
  int nfaces_total_owned = getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  int nfaces_total_all = getNumEntities(Entity_kind::FACE, Parallel_type::ALL);
  Epetra_Map face_owned_map(-1, nfaces_total_owned, 0, *getComm());

  // -- create a map on the background mesh of all FACEs
  Entity_GID_List bg_face_gids = bg_mesh_->getEntityGIDs(Entity_kind::FACE, Parallel_type::ALL);
  Epetra_Map bg_face_all_map(-1, nfaces_bg_all, bg_face_gids.data(), 0, *bg_mesh_->getComm());
  Epetra_Map bg_face_owned_map(-1, nfaces_bg_owned, bg_face_gids.data(), 0, *bg_mesh_->getComm());

  // -- create a vector, fill the owned entries, and scatter
  Epetra_IntVector bg_face_owned_vec(bg_face_owned_map);
  for (int f = 0; f != nfaces_bg_owned; ++f)
    bg_face_owned_vec[f] = face_owned_map.GID(f + nfaces_log + nfaces_extra);
  Epetra_IntVector bg_face_all_vec(bg_face_all_map);
  Epetra_Import face_import(bg_face_all_map, bg_face_owned_map);
  bg_face_all_vec.Import(bg_face_owned_vec, face_import, Insert);

  // -- copy ghost GIDs into the renumbered vector
  Entity_GID_List face_gids(nfaces_total_owned + nfaces_total_all);
  for (int f = 0; f != nfaces_total_owned; ++f) face_gids[f] = face_owned_map.GID(f);
  for (int f = nfaces_total_owned; f != nfaces_total_all; ++f)
    face_gids[f] = bg_face_all_vec[f - nfaces_log - nfaces_extra + nfaces_bg_owned];

  // -- create the map
  face_map_ = Teuchos::rcp(new Epetra_Map(-1, face_gids.size(), face_gids.data(), 0, *getComm()));
}


// Parent entity in the source mesh if mesh was derived from another mesh
Entity_ID
MeshEmbeddedLogical::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
{
  if (kind == Entity_kind::CELL) {
    auto ncells_log = log_mesh_->getNumEntities(kind, Parallel_type::OWNED);
    if (entid < ncells_log)
      return entid;
    else
      return entid - ncells_log;
  } else if (kind == Entity_kind::FACE) {
    auto nfaces_log = log_mesh_->getNumEntities(kind, Parallel_type::OWNED);
    int nfaces_extra = extra_face_cell_ids_.size<MemSpace_type::HOST>();
    if (entid < nfaces_log)
      return entid;
    else if (entid >= (nfaces_log + nfaces_extra))
      return entid - nfaces_log - nfaces_extra;
    else
      return -1;
  } else {
    return -1;
  }
}

Cell_type
MeshEmbeddedLogical::getCellType(const Entity_ID cellid) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  if (cellid < ncells_log)
    return log_mesh_->getCellType(getEntityParent(Entity_kind::CELL, cellid));
  else
    return bg_mesh_->getCellType(getEntityParent(Entity_kind::CELL, cellid));
}

std::size_t
MeshEmbeddedLogical::getNumEntities(const Entity_kind kind, const Parallel_type ptype) const
{
  auto count = log_mesh_->getNumEntities(kind, ptype) + bg_mesh_->getNumEntities(kind, ptype);
  if (kind == Entity_kind::FACE) count += extra_face_cell_ids_.size<MemSpace_type::HOST>();
  return count;
}


Entity_GID
MeshEmbeddedLogical::getEntityGID(const Entity_kind kind, const Entity_ID lid) const
{
  if (getComm()->NumProc() == 1) return lid;
  if (kind == Entity_kind::CELL) {
    return cell_map_->GID(lid);
  } else if (kind == Entity_kind::FACE) {
    return face_map_->GID(lid);
  } else {
    Errors::Message msg;
    msg << "MeshEmbeddedLogical does not have entities of kind " << to_string(kind);
    Exceptions::amanzi_throw(msg);
    return -1;
  }
}


//
// Nodal methods
//
AmanziGeometry::Point
MeshEmbeddedLogical::getNodeCoordinate(const Entity_ID node) const
{
  Errors::Message mesg("There are no nodes in a MeshEmbeddedLogical.");
  Exceptions::amanzi_throw(mesg);
  return AmanziGeometry::Point();
}

void
MeshEmbeddedLogical::getFaceNodes(const Entity_ID f, Entity_ID_List& nodes) const
{
  Errors::Message mesg("There are no nodes in a MeshEmbeddedLogical.");
  Exceptions::amanzi_throw(mesg);
}

void
MeshEmbeddedLogical::getNodeFaces(const Entity_ID nodeid,
                                  const Parallel_type ptype,
                                  Entity_ID_List& faceids) const
{
  Errors::Message mesg("There are no nodes in a MeshEmbeddedLogical.");
  Exceptions::amanzi_throw(mesg);
}


//
// These are the important ones -- MeshLogical defines cell quantities
//
void
MeshEmbeddedLogical::getCellFacesAndDirs(const Entity_ID c,
                                         Entity_ID_List& faces,
                                         Entity_Direction_List* const dirs) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  if (c < ncells_log) {
    log_mesh_->getCellFacesAndDirs(c, faces, dirs);

    // check for extras
    auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
    for (int f = 0; f != extra_face_cell_ids_.size<MemSpace_type::HOST>(); ++f) {
      if (c == extra_face_cell_ids_.get<MemSpace_type::HOST>(f, 0)) {
        faces.push_back(f + nfaces_log);
        if (dirs) dirs->push_back(1); // direction always from log to bg
      }
    }

  } else {
    auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
    auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_type::HOST>();
    bg_mesh_->getCellFacesAndDirs(c - ncells_log, faces, dirs);
    for (int i = 0; i != faces.size(); ++i) faces[i] += (nfaces_log + nfaces_extra);

    // check for extras
    auto nfaces_bg = bg_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
    for (int f = 0; f != extra_face_cell_ids_.size<MemSpace_type::HOST>(); ++f) {
      if (c == (extra_face_cell_ids_.get<MemSpace_type::HOST>(f, 1) + ncells_log)) {
        faces.push_back(f + nfaces_log);
        if (dirs) dirs->push_back(-1); // direction always from log to bg
      }
    }
  }
}

// Get the bisectors, i.e. vectors from cell centroid to face centroids.
void
MeshEmbeddedLogical::getCellFacesAndBisectors(const Entity_ID c,
                                              Entity_ID_List& faces,
                                              Point_List* const bisectors) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  if (c < ncells_log) {
    log_mesh_->getCellFacesAndBisectors(c, faces, bisectors);

    // check for extras
    auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
    for (int f = 0; f != extra_face_cell_ids_.size<MemSpace_type::HOST>(); ++f) {
      if (c == extra_face_cell_ids_.get<MemSpace_type::HOST>(f, 0)) {
        faces.push_back(f + nfaces_log);
        if (bisectors)
          bisectors->push_back(extra_face_cell_bisectors_.get<MemSpace_type::HOST>(f, 0));
      }
    }

  } else {
    auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
    auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_type::HOST>();
    bg_mesh_->getCellFacesAndBisectors(c - ncells_log, faces, bisectors);
    for (int i = 0; i != faces.size(); ++i) faces[i] += (nfaces_log + nfaces_extra);

    // check for extras
    auto nfaces_bg = bg_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
    for (int f = 0; f != extra_face_cell_ids_.size<MemSpace_type::HOST>(); ++f) {
      if (c == (extra_face_cell_ids_.get<MemSpace_type::HOST>(f, 1) + ncells_log)) {
        faces.push_back(f + nfaces_log);
        if (bisectors)
          bisectors->push_back(extra_face_cell_bisectors_.get<MemSpace_type::HOST>(f, 1));
      }
    }
  }
}

double
MeshEmbeddedLogical::getCellVolume(const Entity_ID c) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  if (c < ncells_log)
    return log_mesh_->getCellVolume(c);
  else
    return bg_mesh_->getCellVolume(c - ncells_log);
}

AmanziGeometry::Point
MeshEmbeddedLogical::getCellCentroid(const Entity_ID c) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  if (c < ncells_log)
    return log_mesh_->getCellCentroid(c);
  else
    return bg_mesh_->getCellCentroid(c - ncells_log);
}


void
MeshEmbeddedLogical::getFaceCells(const Entity_ID f,
                                  const Parallel_type ptype,
                                  Entity_ID_List& cells) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_type::HOST>();
  if (f < nfaces_log) {
    log_mesh_->getFaceCells(f, ptype, cells);
  } else if (f >= nfaces_log + nfaces_extra) {
    bg_mesh_->getFaceCells(f - nfaces_log - nfaces_extra, ptype, cells);
    for (auto& c : cells) c += ncells_log;
  } else {
    cells = asVector(extra_face_cell_ids_.getRow<MemSpace_type::HOST>(f - nfaces_log));
    cells[1] += ncells_log;
  }
}

double
MeshEmbeddedLogical::getFaceArea(const Entity_ID f) const
{
  auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_type::HOST>();
  if (f < nfaces_log) {
    return log_mesh_->getFaceArea(f);
  } else if (f >= nfaces_log + nfaces_extra) {
    return bg_mesh_->getFaceArea(f - nfaces_log - nfaces_extra);
  } else {
    return AmanziGeometry::norm(extra_face_area_normals_[f - nfaces_log]);
  }
}

AmanziGeometry::Point
MeshEmbeddedLogical::getFaceCentroid(const Entity_ID f) const
{
  auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_type::HOST>();
  if (f < nfaces_log) {
    return log_mesh_->getFaceCentroid(f);
  } else if (f >= nfaces_log + nfaces_extra) {
    return bg_mesh_->getFaceCentroid(f - nfaces_log - nfaces_extra);
  } else {
    auto log_cc =
      log_mesh_->getCellCentroid(extra_face_cell_ids_.get<MemSpace_type::HOST>(f - nfaces_log, 0));
    auto bg_cc =
      bg_mesh_->getCellCentroid(extra_face_cell_ids_.get<MemSpace_type::HOST>(f - nfaces_log, 1));
    return (log_cc + bg_cc) / 2.;
  }
}

AmanziGeometry::Point
MeshEmbeddedLogical::getFaceNormal(const Entity_ID f,
                                   const Entity_ID c,
                                   int* const orientation) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_type::HOST>();

  if (f < nfaces_log) {
    return log_mesh_->getFaceNormal(f, c, orientation);
  } else if (f >= nfaces_log + nfaces_extra) {
    return bg_mesh_->getFaceNormal(f - nfaces_log - nfaces_extra, c - ncells_log, orientation);
  } else {
    AmanziGeometry::Point normal;
    if (c < ncells_log) {
      normal = extra_face_area_normals_[f - nfaces_log];
      if (orientation) *orientation = 1;
    } else {
      normal = -extra_face_area_normals_[f - nfaces_log];
      if (orientation) *orientation = -1;
    }
    return normal;
  }
}


} // namespace AmanziMesh
} // namespace Amanzi
