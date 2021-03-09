/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

// Default imlementations of MeshFramework.

#include "VerboseObject.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "MeshFramework.hh"
#include "MeshAlgorithms_impl.hh"

namespace Amanzi {
namespace AmanziMesh {

MeshFramework::MeshFramework(const Comm_ptr_type& comm,
                             const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
                             const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : comm_(comm),
    gm_(gm),
    plist_(plist),
    space_dim_(-1),
    manifold_dim_(-1),
    vis_mesh_(Teuchos::null)
{
  if (!plist_.get()) plist_ = Teuchos::rcp(new Teuchos::ParameterList("mesh"));
  vo_ = Teuchos::rcp(new VerboseObject(comm, "MeshFramework", *plist_));
}


Parallel_type
MeshFramework::getEntityPtype(const Entity_kind kind, const Entity_ID entid) const
{
  return entid >= getNumEntities(kind, Parallel_type::OWNED) ? Parallel_type::OWNED : Parallel_type::GHOST;
}

Entity_GID
MeshFramework::getEntityGID(const Entity_kind kind, const Entity_ID lid) const
{
  AMANZI_ASSERT(comm_->NumProc() == 1);
  return lid;
}

Entity_GID_List
MeshFramework::getEntityGIDs(const Entity_kind kind, const Parallel_type ptype) const
{
  Entity_GID_List gids(getNumEntities(kind, ptype));
  for (int lid=0; lid!=gids.size(); ++lid)
    gids[lid] = getEntityGID(kind, lid);
  return gids;
}

Entity_ID MeshFramework::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
{
  Errors::Message msg("Parent entities not implemented.");
  Exceptions::amanzi_throw(msg);
  return -1;
}

Cell_type
MeshFramework::getCellType(const Entity_ID c) const
{
  Entity_ID_List faces;
  getCellFaces(c, faces);
  return getCellType_(c, faces);
}

Cell_type
MeshFramework::getCellType_(const Entity_ID c,
                           const Entity_ID_List& faces) const
{
  if (get_manifold_dimension() == 2) {
    switch (faces.size()) {
    case 3:
      return Cell_type::TRI;
      break;
    case 4:
      return Cell_type::QUAD;
      break;
    default:
      return Cell_type::POLYGON;
    }
  } else if (get_manifold_dimension() == 3) {
    int nquads = 0;
    for (const auto& f : faces) {
      Entity_ID_List fnodes;
      getFaceNodes(f, fnodes);
      if (fnodes.size() == 4) nquads++;
    }

    switch (faces.size()) {
      case 4:
        if (nquads == 0) return Cell_type::TET;
        else return Cell_type::POLYHED;
        break;
      case 5:
        if (nquads == 1) return Cell_type::PYRAMID;
        else if (nquads == 3) return Cell_type::PRISM;
        else return Cell_type::POLYHED;
        break;
      case 6:
        if (nquads == 6) return Cell_type::HEX;
        else return Cell_type::POLYHED;
        break;
      default:
        return Cell_type::POLYHED;
    }
  } else {
    Errors::Message msg;
    msg << "Mesh of manifold_dimension = " << get_manifold_dimension() << " not supported";
    Exceptions::amanzi_throw(msg);
  }
  return Cell_type::UNKNOWN;
}


void
MeshFramework::setNodeCoordinate(const Entity_ID nodeid, const AmanziGeometry::Point& ncoord)
{
  Errors::Message msg("MeshFramework is not deformable.");
  Exceptions::amanzi_throw(msg);
}


Point_List
MeshFramework::getCellCoordinates(const Entity_ID c) const
{
  Entity_ID_List nodes;
  getCellNodes(c, nodes);

  Point_List coords;
  coords.reserve(nodes.size());

  for (const auto& n : nodes) {
    coords.emplace_back(getNodeCoordinate(n));
  }
  return coords;
}

Point_List
MeshFramework::getFaceCoordinates(const Entity_ID f) const
{
  Entity_ID_List nodes;
  getFaceNodes(f, nodes);

  Point_List coords;
  coords.reserve(nodes.size());

  for (const auto& n : nodes) {
    coords.emplace_back(getNodeCoordinate(n));
  }
  return coords;
}

Point_List
MeshFramework::getEdgeCoordinates(const Entity_ID e) const
{
  Entity_ID_List nodes;
  getEdgeNodes(e, nodes);

  Point_List coords;
  coords.reserve(nodes.size());

  for (const auto& n : nodes) {
    coords.emplace_back(getNodeCoordinate(n));
  }
  return coords;
}

//
// Cell geometry
//
std::pair<double, AmanziGeometry::Point>
MeshFramework::computeCellGeometry(const Entity_ID c) const
{
  return MeshAlgorithms::computeCellGeometry(*this, c);
}

double
MeshFramework::getCellVolume(const Entity_ID c) const
{
  return computeCellGeometry(c).first;
}

AmanziGeometry::Point
MeshFramework::getCellCentroid(const Entity_ID c) const
{
  return computeCellGeometry(c).second;
}


//
// face geometry
//
std::tuple<double, AmanziGeometry::Point, Point_List>
MeshFramework::computeFaceGeometry(const Entity_ID f) const
{
  return MeshAlgorithms::computeFaceGeometry(*this, f);
}

double
MeshFramework::getFaceArea(const Entity_ID f) const
{
  return std::get<0>(computeFaceGeometry(f));
}

AmanziGeometry::Point
MeshFramework::getFaceCentroid(const Entity_ID f) const
{
  return std::get<1>(computeFaceGeometry(f));
}

AmanziGeometry::Point
MeshFramework::getFaceNormal(const Entity_ID f, const Entity_ID c, int * const orientation) const
{
  auto geom = computeFaceGeometry(f);
  if (c >= 0) {
    Entity_ID_List cells;
    getFaceCells(f, Parallel_type::ALL, cells);

    for (int i=0; i!=cells.size(); ++i) {
      if (c == cells[i]) {
        int dir = MeshAlgorithms::getFaceDirectionInCell(*this, f, c);
        if (orientation) *orientation = dir;
        return dir * std::get<2>(geom)[i];
      }
    }
  }
  // called with a negative argument, returns the natural normal
  return std::get<2>(geom)[0];
}


//
// Edge geometry
//
std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
MeshFramework::computeEdgeGeometry(const Entity_ID c) const
{
  return MeshAlgorithms::computeEdgeGeometry(*this, c);
}

double
MeshFramework::getEdgeLength(const Entity_ID e) const
{
  return AmanziGeometry::norm(computeEdgeGeometry(e).first);
}

AmanziGeometry::Point
MeshFramework::getEdgeVector(const Entity_ID e, const Entity_ID n, int * const orientation) const
{
  auto geom = computeEdgeGeometry(e);
  if (n >= 0) {
    Entity_ID_List nodes;
    getEdgeNodes(e, nodes);
    if (n == nodes[0]) {
      if (orientation) *orientation = 1;
      return geom.first;
    } else if (n == nodes[1]) {
      if (orientation) *orientation = -1;
      return -geom.first;
    } else {
      AMANZI_ASSERT(0);
    }
  }
  return geom.first;
}

AmanziGeometry::Point
MeshFramework::getEdgeCentroid(const Entity_ID e) const
{
  return computeEdgeGeometry(e).second;
}


//
// bisectors
//
void
MeshFramework::getCellFacesAndBisectors(const Entity_ID cellid,
        Entity_ID_List& faceids, Point_List * const bisectors) const
{
  getCellFaces(cellid, faceids);
  if (bisectors)
    MeshAlgorithms::computeBisectors(*this, cellid, faceids, *bisectors);
}


//
// topology
//
void
MeshFramework::getCellEdges(const Entity_ID c, Entity_ID_List& edges) const
{
  hasEdgesOrThrow_();

  edges.resize(0);
  Entity_ID_List faces, fedges;
  getCellFaces(c, faces);
  for (const auto& f : faces) {
    getFaceEdges(f, fedges);
    for (const auto& e : fedges) {
      if (std::find(edges.begin(), edges.end(), e) == edges.end()) {
        edges.emplace_back(e);
      }
    }
  }
}

void
MeshFramework::getCellNodes(const Entity_ID c, Entity_ID_List& nodes) const
{
  nodes.resize(0);
  Entity_ID_List faces, fnodes;
  getCellFaces(c, faces);
  for (const auto& f : faces) {
    getFaceNodes(f, fnodes);
    for (const auto& n : fnodes) {
      if (std::find(nodes.begin(), nodes.end(), n) == nodes.end()) {
        nodes.emplace_back(n);
      }
    }
  }
}

void
MeshFramework::getFaceEdgesAndDirs(const Entity_ID f,
        Entity_ID_List& edges,
        Entity_Direction_List * const dirs) const {
  hasEdgesOrThrow_();
  throwNotImplemented_("getFaceEdgesAndDirs");
}


void
MeshFramework::getEdgeNodes(const Entity_ID e, Entity_ID_List& nodes) const
{
  hasEdgesOrThrow_();
  throwNotImplemented_("getEdgeNodes");
}


void
MeshFramework::getEdgeCells(const Entity_ID e, const Parallel_type ptype, Entity_ID_List& cells) const
{
  hasEdgesOrThrow_();
  throwNotImplemented_("getEdgeCells");
}

void
MeshFramework::getEdgeFaces(const Entity_ID e, const Parallel_type ptype, Entity_ID_List& faces) const
{
  hasEdgesOrThrow_();
  throwNotImplemented_("getEdgeFaces");
}

void
MeshFramework::getNodeCells(const Entity_ID n, const Parallel_type ptype, Entity_ID_List& cells) const
{
  cells.resize(0);
  Entity_ID_List faces, fcells;
  getNodeFaces(n, Parallel_type::ALL, faces);
  for (const auto& f : faces) {
    getFaceCells(f, ptype, fcells);
    for (const auto& c : fcells) {
      if (std::find(cells.begin(), cells.end(), c) == cells.end()) {
        cells.emplace_back(c);
      }
    }
  }
}

void
MeshFramework::getNodeEdges(const Entity_ID nodeid, const Parallel_type ptype, Entity_ID_List& edgeids) const
{
  hasEdgesOrThrow_();
  throwNotImplemented_("getNodeEdges");
}

void
MeshFramework::getSetEntities(const AmanziGeometry::RegionLabeledSet& region,
        const Entity_kind kind, const Parallel_type ptype, Entity_ID_List& entids) const
{
  throwNotImplemented_("getSetEntities");
}


void
MeshFramework::hasEdgesOrThrow_() const
{
  if (!has_edges()) {
    Errors::Message msg("MeshFramework does not include edges.");
    Exceptions::amanzi_throw(msg);
  }
}

void
MeshFramework::throwNotImplemented_(const std::string& methodname) const
{
  Errors::Message msg;
  msg << "MeshFramework does not yet implement " << methodname;
  Exceptions::amanzi_throw(msg);
}


} // namespace AmanziMesh
} // namespace Amanzi
