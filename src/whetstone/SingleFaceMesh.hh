/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Single face mesh has one less topological dimension than the
  base mesh used in the constructor. Like in the base mesh, edge
  and face are equivalent.
*/

#ifndef AMANZI_WHETSTONE_SINGLE_FACE_MESH_HH_
#define AMANZI_WHETSTONE_SINGLE_FACE_MESH_HH_

#include "Point.hh"

#include "Mesh.hh"
#include "SurfaceCoordinateSystem.hh"

namespace Amanzi {
namespace WhetStone {

class SingleFaceMesh : public AmanziMesh::MeshFramework {
 public:
  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                 int f,
                 const SurfaceCoordinateSystem& coordsys)
  : MeshFramework(mesh->getComm(), Teuchos::null, Teuchos::null) 
  {
    BuildCache_(mesh, f, coordsys);
  }

  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int f)
    : MeshFramework(mesh->getComm(), Teuchos::null, Teuchos::null) 
  {
    const auto& xf = mesh->getFaceCentroid(f);
    const auto& normal = mesh->getFaceNormal(f);
    SurfaceCoordinateSystem coordsys(xf, normal);
    BuildCache_(mesh, f, coordsys);
  }

  ~SingleFaceMesh(){};

  // ---------------------
  // Downward connectivity
  // ---------------------
  void
  getCellNodes(const AmanziMesh::Entity_ID c, AmanziMesh::Entity_ID_List& nodes) const override
  {
    AMANZI_ASSERT(c == 0);
    nodes = cell_node_ids_;
  }

  void
  getFaceNodes(const AmanziMesh::Entity_ID f, AmanziMesh::Entity_ID_List& nodes) const override
  {
    AMANZI_ASSERT(f < nnodes_);
    nodes = face_node_ids_[f];
  }

  void getEdgeNodes(const AmanziMesh::Entity_ID e,
                              AmanziMesh::Entity_ID_List& nodes) const override
  {
    AMANZI_ASSERT(e < nnodes_);
    nodes[0] = face_node_ids_[e][0];
    nodes[1] = face_node_ids_[e][1];
  }

  // -------------------
  // Upward connectivity
  // -------------------
  void getNodeCells(const AmanziMesh::Entity_ID v,
                              const AmanziMesh::Parallel_type ptype,
                              AmanziMesh::Entity_ID_List& cells) const override
  {
    AMANZI_ASSERT(v < nnodes_);
    cells.resize(1, 0);
  }

  void getNodeFaces(const AmanziMesh::Entity_ID v,
                              const AmanziMesh::Parallel_type ptype,
                              AmanziMesh::Entity_ID_List& faces) const override
  {
    AMANZI_ASSERT(v < nnodes_);
    faces.resize(2);
    faces[0] = v;
    faces[1] = (v + 1) % nnodes_;
  }

  // --------
  // Geometry
  // -------- 
  AmanziGeometry::Point getNodeCoordinate(const AmanziMesh::Entity_ID v) const override
  {
    AMANZI_ASSERT(v < nnodes_);
    return cell_coords_[v];
  }

  std::size_t getNumEntities(const AmanziMesh::Entity_kind kind,
                                    const AmanziMesh::Parallel_type ptype) const override
  {
    return (kind == AmanziMesh::Entity_kind::CELL) ? 1 : nnodes_;
  }

  AmanziMesh::Parallel_type getEntityPtype(const AmanziMesh::Entity_kind kind,
                                                     const AmanziMesh::Entity_ID ent) const override
  {
    return AmanziMesh::Parallel_type::OWNED;
  }

  void getCellFacesAndDirs(
    const AmanziMesh::Entity_ID c,
    AmanziMesh::Entity_ID_List& faces,
    AmanziMesh::Entity_Direction_List * const dirs) const override {}

  void getFaceCells(const AmanziMesh::Entity_ID f,
                            const AmanziMesh::Parallel_type ptype,
                            AmanziMesh::Entity_ID_List& cells) const override {}

  bool hasEdges() const { return true; }

 private:
  void BuildCache_(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                   int f,
                   const SurfaceCoordinateSystem& coordsys);

 private:
  int nnodes_;
  AmanziMesh::Entity_ID_List cell_node_ids_;
  std::vector<AmanziGeometry::Point> cell_coords_;
  std::vector<AmanziMesh::Entity_ID_List> face_node_ids_;
};

/* ******************************************************************
* Construct cache data
****************************************************************** */
inline void
SingleFaceMesh::BuildCache_(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                            int f,
                            const SurfaceCoordinateSystem& coordsys)
{
  int d = mesh->getSpaceDimension() - 1;
  setSpaceDimension(d);
  setManifoldDimension(d);

  AmanziMesh::Entity_ID_List fedges, fnodes, enodes, cells;
  fnodes = mesh->getFaceNodes(f);
  nnodes_ = fnodes.size();

  std::vector<int> fdirs(nnodes_, 1);
  if (mesh->hasEdges()) { mesh->getFaceEdgesAndDirs(f, fedges, &fdirs); }

  // single surface cell
  fedges.resize(nnodes_);
  for (int i = 0; i < nnodes_; ++i) fedges[i] = i;

  AmanziGeometry::Point xyz(d);

  // cell nodes
  cell_node_ids_.resize(nnodes_);
  cell_coords_.resize(nnodes_);

  for (int i = 0; i < nnodes_; ++i) {
    cell_node_ids_[i] = i;
    xyz = mesh->getNodeCoordinate(fnodes[i]);
    cell_coords_[i] = coordsys.Project(xyz, true);
  }

  // cell faces
  enodes.resize(2);
  face_node_ids_.resize(nnodes_);

  for (int i = 0; i < nnodes_; ++i) {
    int n0(i), n1((i + 1) % nnodes_);
    enodes[0] = (fdirs[i] > 0) ? n0 : n1;
    enodes[1] = (fdirs[i] > 0) ? n1 : n0;
    face_node_ids_[i] = enodes;
  }
}

} // namespace WhetStone
} // namespace Amanzi

#endif
