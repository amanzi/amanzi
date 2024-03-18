/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "RegionEnumerated.hh"
#include "MeshLogical.hh"
#include "MeshEmbeddedLogical.hh"

namespace Amanzi {
namespace AmanziMesh {


void
MeshEmbeddedLogicalAlgorithms::computeCellFacesAndBisectors(
  const MeshHost& mesh,
  const Entity_ID cellid,
  MeshHost::cEntity_ID_View& faceids,
  MeshHost::cPoint_View* const bisectors) const
{
  static_cast<const MeshEmbeddedLogical*>(mesh.getMeshFramework().get())
    ->getCellFacesAndBisectors(cellid, faceids, bisectors);
}

double
MeshEmbeddedLogicalAlgorithms::computeCellVolume(const MeshHost& mesh, const Entity_ID c) const
{
  return static_cast<const MeshEmbeddedLogical*>(mesh.getMeshFramework().get())->getCellVolume(c);
}

AmanziGeometry::Point
MeshEmbeddedLogicalAlgorithms::computeCellCentroid(const MeshHost& mesh, const Entity_ID c) const
{
  return static_cast<const MeshEmbeddedLogical*>(mesh.getMeshFramework().get())->getCellCentroid(c);
}

double
MeshEmbeddedLogicalAlgorithms::computeFaceArea(const MeshHost& mesh, const Entity_ID f) const
{
  return static_cast<const MeshEmbeddedLogical*>(mesh.getMeshFramework().get())->getFaceArea(f);
}

AmanziGeometry::Point
MeshEmbeddedLogicalAlgorithms::computeFaceCentroid(const MeshHost& mesh, const Entity_ID f) const
{
  return static_cast<const MeshEmbeddedLogical*>(mesh.getMeshFramework().get())->getFaceCentroid(f);
}

AmanziGeometry::Point
MeshEmbeddedLogicalAlgorithms::computeFaceNormal(const MeshHost& mesh,
                                                 const Entity_ID f,
                                                 const Entity_ID c,
                                                 int* const orientation) const
{
  return static_cast<const MeshEmbeddedLogical*>(mesh.getMeshFramework().get())
    ->getFaceNormal(f, c, orientation);
}


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
MeshEmbeddedLogical::MeshEmbeddedLogical(const Comm_ptr_type& comm,
                                         Teuchos::RCP<MeshHost> bg_mesh,
                                         Teuchos::RCP<MeshHost> log_mesh,
                                         const std::vector<Entity_ID_List>& face_cell_ids,
                                         const std::vector<Double_List>& face_cell_lengths,
                                         const Point_List& face_area_normals,
                                         const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : MeshFramework(comm, Teuchos::null, plist),
    extra_face_cell_ids_(face_cell_ids),
    extra_face_cell_lengths_(face_cell_lengths),
    bg_mesh_(bg_mesh),
    log_mesh_(log_mesh)
{
  vectorToView(extra_face_area_normals_, face_area_normals);
  setSpaceDimension(3);
  setManifoldDimension(3);

  // merge and remap to get new global IDs
  int ncells_bg_owned = bg_mesh->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int ncells_bg_all = bg_mesh->getNumEntities(Entity_kind::CELL, Parallel_kind::ALL);
  int ncells_log = log_mesh->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int ncells_total_owned = ncells_bg_owned + ncells_log;
  int ncells_total_all = ncells_bg_all + ncells_log;

  int nfaces_bg_owned = bg_mesh->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  int nfaces_bg_all = bg_mesh->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);
  int nfaces_log = log_mesh->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  int nfaces_extra = extra_face_cell_ids_.size<MemSpace_kind::HOST>();
  int nfaces_total_owned = nfaces_bg_owned + nfaces_log + nfaces_extra;
  int nfaces_total_all = nfaces_bg_all + nfaces_log + nfaces_extra;

  // compute "extra" bisectors
  extra_face_cell_bisectors_.resize(nfaces_extra, 2);

  // now loop over new faces, adding the updates to the cell-ordered versions
  for (int fi = 0; fi != nfaces_extra; ++fi) {
    auto normal0 = extra_face_area_normals_[fi];
    extra_face_cell_bisectors_.getRowUnmanaged<MemSpace_kind::HOST>(fi)[0] =
      face_cell_lengths[fi][0] / AmanziGeometry::norm(normal0) * normal0;

    auto normal1 = extra_face_area_normals_[fi];
    extra_face_cell_bisectors_.getRowUnmanaged<MemSpace_kind::HOST>(fi)[1] =
      face_cell_lengths[fi][1] / AmanziGeometry::norm(normal1) * normal1;
  }
  extra_face_cell_bisectors_.update<MemSpace_kind::DEVICE>();

  // Need to renumber the global IDs.  Use the corresponding mesh to communicate.
  // First cells
  // -- create a map of owned CELLs
  auto cell_owned_map = Teuchos::rcp(new Map_type(-1, ncells_total_owned, 0, getComm()));

  // -- create a map on the background mesh of all CELLs
  auto [bg_cell_all_map, bg_cell_owned_map] = createMapsFromMeshGIDs(*this, Entity_kind::CELL);

  // -- create a vector, fill the owned entries, and scatter
  IntVector_type bg_cell_owned_vec(bg_cell_owned_map);
  for (int c = 0; c != ncells_bg_owned; ++c)
    bg_cell_owned_vec.replaceLocalValue(c, cell_owned_map->getGlobalElement(c + ncells_log));
  IntVector_type bg_cell_all_vec(bg_cell_all_map);
  Import_type cell_import(bg_cell_all_map, bg_cell_owned_map);
  bg_cell_all_vec.doImport(bg_cell_owned_vec, cell_import, Tpetra::INSERT);

  // -- copy ghost GIDs into the renumbered vector
  MeshFramework::Entity_GID_View cell_gids("cell_gids", ncells_bg_all + ncells_log);
  for (int c = 0; c != ncells_total_owned; ++c) cell_gids[c] = cell_owned_map->getGlobalElement(c);
  {
    auto bg_cell_all_vec_data = bg_cell_all_vec.getData();
    for (int c = ncells_total_owned; c != ncells_total_all; ++c)
      cell_gids[c] = bg_cell_all_vec_data[c - ncells_log + ncells_bg_owned];
  }

  Kokkos::View<GO*, Amanzi::DefaultExecutionSpace> cell_gids_d("cell_gids on device", cell_gids.size());
  Kokkos::deep_copy(cell_gids_d, cell_gids);

  // -- create the map
  cell_map_ = Teuchos::rcp(new Map_type(-1, cell_gids_d, 0, getComm()));

  // Next faces
  // -- create a map of owned FACES
  auto face_owned_map = Teuchos::rcp(new Map_type(-1, nfaces_total_owned, 0, getComm()));

  // -- create a map on the background mesh of all FACEs
  auto [bg_face_all_map, bg_face_owned_map] = createMapsFromMeshGIDs(*this, Entity_kind::FACE);

  // -- create a vector, fill the owned entries, and scatter
  IntVector_type bg_face_owned_vec(bg_face_owned_map);
  for (int f = 0; f != nfaces_bg_owned; ++f)
    bg_face_owned_vec.replaceLocalValue(
      f, face_owned_map->getGlobalElement(f + nfaces_log + nfaces_extra));
  IntVector_type bg_face_all_vec(bg_face_all_map);
  Import_type face_import(bg_face_all_map, bg_face_owned_map);
  bg_face_all_vec.doImport(bg_face_owned_vec, face_import, Tpetra::INSERT);

  // -- copy ghost GIDs into the renumbered vector
  MeshFramework::Entity_GID_View face_gids("face_gids", nfaces_total_owned + nfaces_total_all);
  for (int f = 0; f != nfaces_total_owned; ++f) face_gids[f] = face_owned_map->getGlobalElement(f);
  {
    auto bg_face_all_vec_data = bg_face_all_vec.getData();
    for (int f = nfaces_total_owned; f != nfaces_total_all; ++f)
      face_gids[f] = bg_face_all_vec_data[f - nfaces_log - nfaces_extra + nfaces_bg_owned];
  }

  Kokkos::View<GO*, Amanzi::DefaultExecutionSpace> face_gids_d("face_gids on device", face_gids.size());
  Kokkos::deep_copy(face_gids_d, face_gids);

  // -- create the map
  face_map_ = Teuchos::rcp(new Map_type(-1, face_gids_d, 0, getComm()));
}


// Parent entity in the source mesh if mesh was derived from another mesh
Entity_ID
MeshEmbeddedLogical::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
{
  if (kind == Entity_kind::CELL) {
    auto ncells_log = log_mesh_->getNumEntities(kind, Parallel_kind::OWNED);
    if (entid < ncells_log)
      return entid;
    else
      return entid - ncells_log;
  } else if (kind == Entity_kind::FACE) {
    auto nfaces_log = log_mesh_->getNumEntities(kind, Parallel_kind::OWNED);
    int nfaces_extra = extra_face_cell_ids_.size<MemSpace_kind::HOST>();
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

Cell_kind
MeshEmbeddedLogical::getCellType(const Entity_ID cellid) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  if (cellid < ncells_log)
    return log_mesh_->getCellType(getEntityParent(Entity_kind::CELL, cellid));
  else
    return bg_mesh_->getCellType(getEntityParent(Entity_kind::CELL, cellid));
}

std::size_t
MeshEmbeddedLogical::getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const
{
  auto count = log_mesh_->getNumEntities(kind, ptype) + bg_mesh_->getNumEntities(kind, ptype);
  if (kind == Entity_kind::FACE) count += extra_face_cell_ids_.size<MemSpace_kind::HOST>();
  return count;
}


Entity_GID
MeshEmbeddedLogical::getEntityGID(const Entity_kind kind, const Entity_ID lid) const
{
  if (getComm()->getSize() == 1) return lid;
  if (kind == Entity_kind::CELL) {
    return cell_map_->getGlobalElement(lid);
  } else if (kind == Entity_kind::FACE) {
    return face_map_->getGlobalElement(lid);
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
MeshEmbeddedLogical::getFaceNodes(const Entity_ID f, MeshHost::cEntity_ID_View& nodes) const
{
  Errors::Message mesg("There are no nodes in a MeshEmbeddedLogical.");
  Exceptions::amanzi_throw(mesg);
}

void
MeshEmbeddedLogical::getNodeFaces(const Entity_ID nodeid, MeshHost::cEntity_ID_View& faceids) const
{
  Errors::Message mesg("There are no nodes in a MeshEmbeddedLogical.");
  Exceptions::amanzi_throw(mesg);
}


//
// These are the important ones -- MeshLogical defines cell quantities
//
void
MeshEmbeddedLogical::getCellFacesAndDirs(const Entity_ID c,
                                         MeshHost::cEntity_ID_View& faces,
                                         MeshHost::cDirection_View* const dirs) const
{
  MeshHost::Entity_ID_View lfaces;
  MeshHost::Direction_View ldirs;
  Entity_ID_List vfaces;
  std::vector<int> vdirs;
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  if (c < ncells_log) {
    log_mesh_->getCellFacesAndDirs(c, faces, dirs);
    lfaces.fromConst(faces);
    if (dirs) ldirs.fromConst(*dirs);
    // check for extras
    auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
    for (int f = 0; f != extra_face_cell_ids_.size<MemSpace_kind::HOST>(); ++f) {
      if (c == extra_face_cell_ids_.get<MemSpace_kind::HOST>(f, 0)) {
        vfaces.push_back(f + nfaces_log);
        if (dirs) vdirs.push_back(1); // direction always from log to bg
      }
    }
    int size = lfaces.size();
    Kokkos::resize(lfaces, size + vfaces.size());
    for (int i = size, j = 0; i < lfaces.size(); ++i, ++j) lfaces[i] = vfaces[j];
    if (dirs) {
      int lsize = ldirs.size();
      Kokkos::resize(ldirs, lsize + vdirs.size());
      for (int i = lsize, j = 0; i < ldirs.size(); ++i, ++j) ldirs[i] = vdirs[j];
      *dirs = ldirs;
    }

  } else {
    auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
    auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_kind::HOST>();
    bg_mesh_->getCellFacesAndDirs(c - ncells_log, faces, dirs);
    if (dirs) ldirs.fromConst(*dirs);
    lfaces.fromConst(faces);
    for (int i = 0; i != faces.size(); ++i) lfaces[i] += (nfaces_log + nfaces_extra);

    // check for extras
    for (int f = 0; f != extra_face_cell_ids_.size<MemSpace_kind::HOST>(); ++f) {
      if (c == (extra_face_cell_ids_.get<MemSpace_kind::HOST>(f, 1) + ncells_log)) {
        vfaces.push_back(f + nfaces_log);
        if (dirs) vdirs.push_back(-1); // direction always from log to bg
      }
    }
    int size = lfaces.size();
    Kokkos::resize(lfaces, size + vfaces.size());
    for (int i = size, j = 0; i < lfaces.size(); ++i, ++j) lfaces[i] = vfaces[j];
    if (dirs) {
      int lsize = ldirs.size();
      Kokkos::resize(ldirs, lsize + vdirs.size());
      for (int i = lsize, j = 0; i < ldirs.size(); ++i, ++j) ldirs[i] = vdirs[j];
      *dirs = ldirs;
    }
  }
  faces = lfaces;
}

// Get the bisectors, i.e. vectors from cell centroid to face centroids.
void
MeshEmbeddedLogical::getCellFacesAndBisectors(const Entity_ID c,
                                              MeshHost::cEntity_ID_View& faces,
                                              MeshHost::cPoint_View* const bisectors) const
{
  MeshHost::Entity_ID_View lfaces;
  MeshHost::Point_View lbis;

  Entity_ID_List vfaces;
  Point_List vbisectors;
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  if (c < ncells_log) {
    log_mesh_->getCellFacesAndBisectors(c, faces, bisectors);
    lbis.fromConst(*bisectors);
    lfaces.fromConst(faces);

    // check for extras
    auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
    for (int f = 0; f != extra_face_cell_ids_.size<MemSpace_kind::HOST>(); ++f) {
      if (c == extra_face_cell_ids_.get<MemSpace_kind::HOST>(f, 0)) {
        vfaces.push_back(f + nfaces_log);
        if (bisectors)
          vbisectors.push_back(extra_face_cell_bisectors_.get<MemSpace_kind::HOST>(f, 0));
      }
    }
    int size = lfaces.size();
    Kokkos::resize(lfaces, size + vfaces.size());
    for (int i = size, j = 0; i < lfaces.size(); ++i, ++j) lfaces[i] = vfaces[j];
    if (bisectors) {
      int lsize = lbis.size();
      Kokkos::resize(lbis, lsize + vbisectors.size());
      for (int i = lsize, j = 0; i < lbis.size(); ++i, ++j) lbis[i] = vbisectors[j];
      *bisectors = lbis;
    }
  } else {
    auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
    auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_kind::HOST>();
    bg_mesh_->getCellFacesAndBisectors(c - ncells_log, faces, bisectors);
    lbis.fromConst(*bisectors);
    lfaces.fromConst(faces);
    for (int i = 0; i != lfaces.size(); ++i) lfaces[i] += (nfaces_log + nfaces_extra);

    // check for extras
    for (int f = 0; f != extra_face_cell_ids_.size<MemSpace_kind::HOST>(); ++f) {
      if (c == (extra_face_cell_ids_.get<MemSpace_kind::HOST>(f, 1) + ncells_log)) {
        vfaces.push_back(f + nfaces_log);
        if (bisectors)
          vbisectors.push_back(extra_face_cell_bisectors_.get<MemSpace_kind::HOST>(f, 1));
      }
    }
    int size = lfaces.size();
    Kokkos::resize(lfaces, size + vfaces.size());
    for (int i = size, j = 0; i < lfaces.size(); ++i, ++j) lfaces[i] = vfaces[j];
    if (bisectors) {
      int lsize = lbis.size();
      Kokkos::resize(lbis, lsize + vbisectors.size());
      for (int i = lsize, j = 0; i < lbis.size(); ++i, ++j) lbis[i] = vbisectors[j];
      *bisectors = lbis;
    }
  }
  faces = lfaces;
}

double
MeshEmbeddedLogical::getCellVolume(const Entity_ID c) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  if (c < ncells_log)
    return log_mesh_->getCellVolume(c);
  else
    return bg_mesh_->getCellVolume(c - ncells_log);
}

AmanziGeometry::Point
MeshEmbeddedLogical::getCellCentroid(const Entity_ID c) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  if (c < ncells_log)
    return log_mesh_->getCellCentroid(c);
  else
    return bg_mesh_->getCellCentroid(c - ncells_log);
}


void
MeshEmbeddedLogical::getFaceCells(const Entity_ID f, MeshHost::cEntity_ID_View& cells) const
{
  MeshHost::Entity_ID_View lcells;
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_kind::HOST>();
  if (f < nfaces_log) {
    log_mesh_->getFaceCells(f, cells);
  } else if (f >= (nfaces_log + nfaces_extra)) {
    bg_mesh_->getFaceCells(f - nfaces_log - nfaces_extra, cells);
    lcells.fromConst(cells);
    for (auto& c : lcells) c += ncells_log;
    cells = lcells;
  } else {
    lcells.fromConst(extra_face_cell_ids_.getRow<MemSpace_kind::HOST>(f - nfaces_log));
    lcells[1] += ncells_log;
    cells = lcells;
  }
}

double
MeshEmbeddedLogical::getFaceArea(const Entity_ID f) const
{
  auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_kind::HOST>();
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
  auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_kind::HOST>();
  if (f < nfaces_log) {
    return log_mesh_->getFaceCentroid(f);
  } else if (f >= nfaces_log + nfaces_extra) {
    return bg_mesh_->getFaceCentroid(f - nfaces_log - nfaces_extra);
  } else {
    auto log_cc =
      log_mesh_->getCellCentroid(extra_face_cell_ids_.get<MemSpace_kind::HOST>(f - nfaces_log, 0));
    auto bg_cc =
      bg_mesh_->getCellCentroid(extra_face_cell_ids_.get<MemSpace_kind::HOST>(f - nfaces_log, 1));
    return (log_cc + bg_cc) / 2.;
  }
}

AmanziGeometry::Point
MeshEmbeddedLogical::getFaceNormal(const Entity_ID f,
                                   const Entity_ID c,
                                   int* const orientation) const
{
  auto ncells_log = log_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  auto nfaces_log = log_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  auto nfaces_extra = extra_face_cell_ids_.size<MemSpace_kind::HOST>();

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
