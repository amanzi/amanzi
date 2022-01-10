/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Single face mesh has one less topological dimension than the 
  base mesh used in the constructor. Like in the base mesh, edge
  and face are equivalent.
*/

#ifndef AMANZI_WHETSTONE_SINGLE_FACE_MESH_HH_
#define AMANZI_WHETSTONE_SINGLE_FACE_MESH_HH_

#include "Point.hh"

#include "MeshFramework.hh"
#include "SurfaceCoordinateSystem.hh"

namespace Amanzi {
namespace WhetStone {

class SingleFaceMesh : public AmanziMesh::MeshFramework {
 public:
  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int f,
                 const SurfaceCoordinateSystem& coordsys)
      : MeshFramework(Amanzi::getDefaultComm(), Teuchos::null, Teuchos::null),
        mesh_(mesh),
        coordsys_(coordsys),
        f_(f) {
    Init_(mesh, f, coordsys);
  }

  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int f)
      : MeshFramework(Amanzi::getDefaultComm(), Teuchos::null, Teuchos::null),
        mesh_(mesh),
        coordsys_(SurfaceCoordinateSystem(mesh->getFaceCentroid(f), mesh->getFaceNormal(f))),
        f_(f) {
    Init_(mesh, f, coordsys_);
   }

  ~SingleFaceMesh() {};

  // ---------------------
  // Downward connectivity
  // ---------------------
  virtual
  void getCellNodes(const AmanziMesh::Entity_ID c,
                    AmanziMesh::Entity_ID_List& nodes) const override;

  virtual
  void getFaceNodes(const AmanziMesh::Entity_ID f,
                    AmanziMesh::Entity_ID_List& nodes) const override;

  virtual 
  void getCellFacesAndDirs(const AmanziMesh::Entity_ID c,
                           AmanziMesh::Entity_ID_List& faces,
                           AmanziMesh::Entity_Direction_List* dirs) const final;

  // -------------------
  // Upward connectivity
  // -------------------
  virtual
  void getFaceCells(const AmanziMesh::Entity_ID f,
                    const AmanziMesh::Parallel_type ptype,
                    AmanziMesh::Entity_ID_List& cells) const final;

  virtual
  void getNodeFaces(const AmanziMesh::Entity_ID v,
                    const AmanziMesh::Parallel_type ptype,
                    AmanziMesh::Entity_ID_List& faces) const override;

  // --------
  // Geometry
  // --------
  virtual
  AmanziGeometry::Point getNodeCoordinate(const AmanziGeometry::Entity_ID v) const final;

  virtual 
  std::size_t getNumEntities(const AmanziMesh::Entity_kind kind,
                             const AmanziMesh::Parallel_type ptype) const override;

  virtual
  AmanziMesh::Parallel_type getEntityPtype(const AmanziMesh::Entity_kind kind, 
                                           const AmanziMesh::Entity_ID entid) const final {
    return AmanziMesh::Parallel_type::OWNED;
  }

 private:
  void Init_(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int f,
    const SurfaceCoordinateSystem& coordsys);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  SurfaceCoordinateSystem coordsys_;
  const AmanziMesh::Entity_ID f_;

  int nnodes_;
};


inline
void SingleFaceMesh::getCellNodes(
    const AmanziMesh::Entity_ID c,
    AmanziMesh::Entity_ID_List& nodes) const
{
  AMANZI_ASSERT(c == 0);
  nodes.resize(nnodes_);
  for (int i = 0; i < nnodes_; ++i) nodes[i] = i;
}


inline
void SingleFaceMesh::getCellFacesAndDirs(
    const AmanziMesh::Entity_ID c,
    AmanziMesh::Entity_ID_List& faces,
    AmanziMesh::Entity_Direction_List* const dirs) const
{
  AmanziMesh::Entity_ID_List fedges;
  mesh_->getFaceEdgesAndDirs(f_, fedges, dirs);

  faces.resize(nnodes_);
  for (int i = 0; i < nnodes_; ++i) faces[i] = i;
}


inline
void SingleFaceMesh::getFaceNodes(
    const AmanziMesh::Entity_ID f,
    AmanziMesh::Entity_ID_List& nodes) const 
{
  AMANZI_ASSERT(f < nnodes_);
  AmanziMesh::Entity_ID_List fedges;
  AmanziMesh::Entity_Direction_List fdirs;
  mesh_->getFaceEdgesAndDirs(f_, fedges, &fdirs);

  nodes.resize(2);
  int n0(f), n1((f + 1) % nnodes_);

  nodes[0] = (fdirs[f] > 0) ? n0 : n1;
  nodes[1] = (fdirs[f] > 0) ? n1 : n0;
}


inline
void SingleFaceMesh::getFaceCells(
    const AmanziMesh::Entity_ID f,
    const AmanziMesh::Parallel_type ptype,
    AmanziMesh::Entity_ID_List& cells) const
{
  cells.resize(1, 0);
}


inline
AmanziGeometry::Point SingleFaceMesh::getNodeCoordinate(
    const AmanziGeometry::Entity_ID v) const
{
  AMANZI_ASSERT(v < nnodes_);
  const auto& nodes = mesh_->getFaceNodes(f_);
  const auto& xyz = mesh_->getNodeCoordinate(nodes[v]);
  return coordsys_.Project(xyz, true);
}


inline
void SingleFaceMesh::getNodeFaces(
    const AmanziMesh::Entity_ID v,
    const AmanziMesh::Parallel_type ptype,
   AmanziMesh::Entity_ID_List& faces) const 
{
  AMANZI_ASSERT(v < nnodes_);
  faces.resize(2);
  faces[0] = v; 
  faces[1] = (v + 1) % nnodes_;
}


inline
std::size_t SingleFaceMesh::getNumEntities(
    const AmanziMesh::Entity_kind kind,
    const AmanziMesh::Parallel_type ptype) const
{
  if (kind == AmanziMesh::Entity_kind::CELL) {
    return 1;
  } else if (kind == AmanziMesh::Entity_kind::NODE || 
             kind == AmanziMesh::Entity_kind::EDGE || 
             kind == AmanziMesh::Entity_kind::FACE) {
    return nnodes_;
  }
  return 0;
}


/* ******************************************************************
* Initialize framework from a cache
****************************************************************** */
inline
void SingleFaceMesh::Init_(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int f,
    const SurfaceCoordinateSystem& coordsys)
{
  int d = mesh->getSpaceDimension() - 1;
  nnodes_ = mesh_->getFaceNodes(f).size();

  setSpaceDimension(d);
  setManifoldDimension(d);
}

} // namespace WhetStone
} // namespace Amanzi

#endif
