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

#ifndef AMANZI_MESH_SINGLE_FACE_MESH_HH_
#define AMANZI_MESH_SINGLE_FACE_MESH_HH_

#include "Point.hh"

#include "Mesh.hh"
#include "SurfaceCoordinateSystem.hh"

namespace Amanzi {
namespace AmanziMesh {

class SingleFaceMesh : public AmanziMesh::MeshFramework {
 public:
  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                 int f,
                 const AmanziGeometry::SurfaceCoordinateSystem& coordsys)
  : MeshFramework(mesh->getComm(), Teuchos::null, Teuchos::null) 
  {
    BuildCache_(mesh, f, coordsys);
  }

  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int f)
    : MeshFramework(mesh->getComm(), Teuchos::null, Teuchos::null) 
  {
    const auto& xf = mesh->getFaceCentroid(f);
    const auto& normal = mesh->getFaceNormal(f);
    AmanziGeometry::SurfaceCoordinateSystem coordsys(xf, normal);
    BuildCache_(mesh, f, coordsys);
  }

  ~SingleFaceMesh(){};

  // ---------------------
  // Downward connectivity
  // ---------------------
  void
  getCellNodes(const AmanziMesh::Entity_ID c, AmanziMesh::cEntity_ID_View& nodes) const override
  {
    AMANZI_ASSERT(c == 0);
    nodes = cell_node_ids_;
  }

  void
  getFaceNodes(const AmanziMesh::Entity_ID f, AmanziMesh::cEntity_ID_View& nodes) const override
  {
    AMANZI_ASSERT(f < nnodes_);
    nodes = face_node_ids_[f];
  }

  void getEdgeNodes(const AmanziMesh::Entity_ID e,
                              AmanziMesh::cEntity_ID_View& nodes) const override
  {
    AMANZI_ASSERT(e < nnodes_);
    AmanziMesh::Entity_ID_View lnodes("nodes", 2); 
    lnodes[0] = face_node_ids_[e][0];
    lnodes[1] = face_node_ids_[e][1];
    nodes = lnodes; 
  }

  // -------------------
  // Upward connectivity
  // -------------------
  void getNodeCells(const AmanziMesh::Entity_ID v,
                              const AmanziMesh::Parallel_kind ptype,
                              AmanziMesh::cEntity_ID_View& cells) const override
  {
    AMANZI_ASSERT(v < nnodes_);
    AmanziMesh::Entity_ID_View lcells; 
    Kokkos::resize(lcells, 1);
    cells = lcells; 
  }

  void getNodeFaces(const AmanziMesh::Entity_ID v,
                              const AmanziMesh::Parallel_kind ptype,
                              AmanziMesh::cEntity_ID_View& faces) const override
  {
    AMANZI_ASSERT(v < nnodes_);
    AmanziMesh::Entity_ID_View lfaces("faces", 2); 
    lfaces[0] = v;
    lfaces[1] = (v + 1) % nnodes_;
    faces = lfaces;
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
                                    const AmanziMesh::Parallel_kind ptype) const override
  {
    return (kind == AmanziMesh::Entity_kind::CELL) ? 1 : nnodes_;
  }

  AmanziMesh::Parallel_kind getEntityPtype(const AmanziMesh::Entity_kind kind,
                                                     const AmanziMesh::Entity_ID ent) const override
  {
    return AmanziMesh::Parallel_kind::OWNED;
  }

  void getCellFacesAndDirs(
    const AmanziMesh::Entity_ID c,
    AmanziMesh::cEntity_ID_View& faces,
    AmanziMesh::cEntity_Direction_View * const dirs) const override {}

  void getFaceCells(const AmanziMesh::Entity_ID f,
                            const AmanziMesh::Parallel_kind ptype,
                            AmanziMesh::cEntity_ID_View& cells) const override {}

  bool hasEdges() const { return true; }

 private:
  void BuildCache_(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                   int f,
                   const AmanziGeometry::SurfaceCoordinateSystem& coordsys);

 private:
  int nnodes_;
  AmanziMesh::Entity_ID_View cell_node_ids_;
  Point_List cell_coords_;
  std::vector<AmanziMesh::Entity_ID_View> face_node_ids_;
};

/* ******************************************************************
* Construct cache data
****************************************************************** */
inline void
SingleFaceMesh::BuildCache_(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                            int f,
                            const AmanziGeometry::SurfaceCoordinateSystem& coordsys)
{
  int d = mesh->getSpaceDimension() - 1;
  setSpaceDimension(d);
  setManifoldDimension(d);

  AmanziMesh::cEntity_ID_View fnodes, fedges; 
  AmanziMesh::Entity_ID_View enodes, lfnodes;
  fnodes = mesh->getFaceNodes(f);
  nnodes_ = fnodes.size();

  AmanziMesh::Entity_Direction_View fdirs;
  if (mesh->hasEdges()) { 
    AmanziMesh::cEntity_Direction_View cfdirs;
    mesh->getFaceEdgesAndDirs(f, fedges, &cfdirs); 
    fdirs.fromConst(cfdirs); 
  } 
  else 
    Kokkos::resize(fdirs, 1); 

  // cell nodes
  Kokkos::resize(cell_node_ids_, nnodes_);
  cell_coords_.resize(nnodes_);

  for (int i = 0; i < nnodes_; ++i) {
    cell_node_ids_[i] = i;
    auto xyz = mesh->getNodeCoordinate(fnodes[i]);
    cell_coords_[i] = coordsys.Project(xyz, true);
  }

  // cell faces
  Kokkos::resize(enodes, 2);
  face_node_ids_.resize(nnodes_);

  for (int i = 0; i < nnodes_; ++i) {
    int n0(i), n1((i + 1) % nnodes_);
    enodes[0] = (fdirs[i] > 0) ? n0 : n1;
    enodes[1] = (fdirs[i] > 0) ? n1 : n0;
    face_node_ids_[i] = enodes;
  }
}

} // namespace AmanziMesh
} // namespace Amanzi

#endif
