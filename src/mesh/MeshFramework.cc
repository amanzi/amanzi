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
#include "Mesh_Helpers_impl.hh"

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
  algorithms_ = Teuchos::rcp(new MeshFrameworkAlgorithms());
}


Parallel_kind
MeshFramework::getEntityPtype(const Entity_kind kind, const Entity_ID entid) const
{
  return entid >= getNumEntities(kind, Parallel_kind::OWNED) ? Parallel_kind::OWNED : Parallel_kind::GHOST;
}

Entity_GID
MeshFramework::getEntityGID(const Entity_kind kind, const Entity_ID lid) const
{
  AMANZI_ASSERT(comm_->NumProc() == 1);
  return lid;
}

Entity_GID_View
MeshFramework::getEntityGIDs(const Entity_kind kind, const Parallel_kind ptype) const
{
  Entity_GID_View gids("gids",getNumEntities(kind, ptype));
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

Cell_kind
MeshFramework::getCellType(const Entity_ID c) const
{
  Entity_ID_View faces;
  getCellFaces(c, faces);
  return getCellType_(c, faces);
}

Cell_kind
MeshFramework::getCellType_(const Entity_ID c,
                           const Entity_ID_View& faces) const
{
  if (getManifoldDimension() == 2) {
    switch (faces.size()) {
    case 3:
      return Cell_kind::TRI;
      break;
    case 4:
      return Cell_kind::QUAD;
      break;
    default:
      return Cell_kind::POLYGON;
    }
  } else if (getManifoldDimension() == 3) {
    int nquads = 0;
    for (const auto& f : faces) {
      Entity_ID_View fnodes;
      getFaceNodes(f, fnodes);
      if (fnodes.size() == 4) nquads++;
    }

    switch (faces.size()) {
      case 4:
        if (nquads == 0) return Cell_kind::TET;
        else return Cell_kind::POLYHED;
        break;
      case 5:
        if (nquads == 1) return Cell_kind::PYRAMID;
        else if (nquads == 3) return Cell_kind::PRISM;
        else return Cell_kind::POLYHED;
        break;
      case 6:
        if (nquads == 6) return Cell_kind::HEX;
        else return Cell_kind::POLYHED;
        break;
      default:
        return Cell_kind::POLYHED;
    }
  } else {
    Errors::Message msg;
    msg << "Mesh of manifold_dimension = " << getManifoldDimension() << " not supported";
    Exceptions::amanzi_throw(msg);
  }
  return Cell_kind::UNKNOWN;
}


void
MeshFramework::setNodeCoordinate(const Entity_ID nodeid, const AmanziGeometry::Point& ncoord)
{
  Errors::Message msg("MeshFramework is not deformable.");
  Exceptions::amanzi_throw(msg);
}


Point_View
MeshFramework::getCellCoordinates(const Entity_ID c) const
{
  return MeshAlgorithms::getCellCoordinates(*this, c);
}

Point_View
MeshFramework::getFaceCoordinates(const Entity_ID f) const
{
  return MeshAlgorithms::getFaceCoordinates(*this, f);
}

Point_View
MeshFramework::getEdgeCoordinates(const Entity_ID e) const
{
  return MeshAlgorithms::getEdgeCoordinates(*this, e);
}

//
// Cell geometry
//
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

  Entity_ID_View fcells;
  getFaceCells(f, Parallel_kind::ALL, fcells);
  if (orientation) *orientation = 0;
  Entity_ID cc = (c < 0) ? fcells[0] : c;

  int i = std::find(fcells.begin(), fcells.end(), cc) - fcells.begin();
  AmanziGeometry::Point normal = std::get<2>(geom)[i];

  if (getSpaceDimension() == getManifoldDimension()) {
    if (c < 0) {
      normal *= MeshAlgorithms::getFaceDirectionInCell(*this, f, cc);
    } else if (orientation) {
      *orientation = MeshAlgorithms::getFaceDirectionInCell(*this, f, cc);
    }
  } else if (c < 0) {
    Errors::Message msg("MeshFramework: asking for the natural normal of a submanifold mesh is not valid.");
    Exceptions::amanzi_throw(msg);
  }
  return normal;
}


//
// Edge geometry
//
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
    Entity_ID_View nodes;
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
        Entity_ID_View& faceids, Point_View * const bisectors) const
{
  getCellFaces(cellid, faceids);
  if (bisectors)
    *bisectors = MeshAlgorithms::computeBisectors(*this, cellid, faceids);
}


//
// topology
//
void
MeshFramework::getCellEdges(const Entity_ID c, Entity_ID_View& edges) const
{
  hasEdgesOrThrow();
  edges = MeshAlgorithms::computeCellEdges(*this, c);
}

void
MeshFramework::getCellNodes(const Entity_ID c, Entity_ID_View& nodes) const
{
  nodes = MeshAlgorithms::computeCellNodes(*this, c);
}

void
MeshFramework::getFaceEdgesAndDirs(const Entity_ID f,
        Entity_ID_View& edges,
        Entity_Direction_View * const dirs) const {
  hasEdgesOrThrow();
  throwNotImplemented_("getFaceEdgesAndDirs");
}


void
MeshFramework::getEdgeNodes(const Entity_ID e, Entity_ID_View& nodes) const
{
  hasEdgesOrThrow();
  throwNotImplemented_("getEdgeNodes");
}


void
MeshFramework::getEdgeCells(const Entity_ID e, const Parallel_kind ptype, Entity_ID_View& cells) const
{
  hasEdgesOrThrow();
  throwNotImplemented_("getEdgeCells");
}

void
MeshFramework::getEdgeFaces(const Entity_ID e, const Parallel_kind ptype, Entity_ID_View& faces) const
{
  hasEdgesOrThrow();
  throwNotImplemented_("getEdgeFaces");
}

void
MeshFramework::getNodeCells(const Entity_ID n, const Parallel_kind ptype, Entity_ID_View& cells) const
{
  Entity_ID_View faces, fcells;
  Entity_ID_List vcells; 
  getNodeFaces(n, Parallel_kind::ALL, faces);
  for (const auto& f : faces) {
    getFaceCells(f, ptype, fcells);
    for (const auto& c : fcells) {
      if (std::find(vcells.begin(), vcells.end(), c) == vcells.end()) {
        vcells.emplace_back(c);
      }
    }
  }
  vectorToView(cells,vcells); 
}

void
MeshFramework::getNodeEdges(const Entity_ID nodeid, const Parallel_kind ptype, Entity_ID_View& edgeids) const
{
  hasEdgesOrThrow();
  throwNotImplemented_("getNodeEdges");
}

void
MeshFramework::getSetEntities(const AmanziGeometry::RegionLabeledSet& region,
        const Entity_kind kind, const Parallel_kind ptype, Entity_ID_View& entids) const
{
  throwNotImplemented_("getSetEntities");
}


void
MeshFramework::hasEdgesOrThrow() const
{
  if (!hasEdges()) {
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


std::pair<double, AmanziGeometry::Point>
MeshFrameworkAlgorithms::computeCellGeometry(const MeshFramework& mesh, const Entity_ID c) const
{
  return MeshAlgorithms::computeCellGeometry(mesh, c);
}
// std::pair<double, AmanziGeometry::Point>
// MeshFrameworkAlgorithms::computeCellGeometry(const Mesh& mesh, const Entity_ID c) const
// {
//   return MeshAlgorithms::computeCellGeometry(mesh, c);
// }

std::tuple<double, AmanziGeometry::Point, Point_View>
MeshFrameworkAlgorithms::computeFaceGeometry(const MeshFramework& mesh, const Entity_ID f) const
{
  return MeshAlgorithms::computeFaceGeometry(mesh, f);
}
// std::tuple<double, AmanziGeometry::Point, Point_View>
// MeshFrameworkAlgorithms::computeFaceGeometry(const Mesh& mesh, const Entity_ID f) const
// {
//   return MeshAlgorithms::computeFaceGeometry(mesh, f);
// }

std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
MeshFrameworkAlgorithms::computeEdgeGeometry(const MeshFramework& mesh, const Entity_ID c) const
{
  return MeshAlgorithms::computeEdgeGeometry(mesh, c);
}
// std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
// MeshFrameworkAlgorithms::computeEdgeGeometry(const Mesh& mesh, const Entity_ID c) const
// {
//   return MeshAlgorithms::computeEdgeGeometry(mesh, c);
// }

// Point_View
// MeshFrameworkAlgorithms::computeBisectors(const Mesh& mesh, const Entity_ID c, const Entity_ID_View& faces) const
// {
//   return MeshAlgorithms::computeBisectors(mesh, c, faces);
// }




} // namespace AmanziMesh
} // namespace Amanzi
