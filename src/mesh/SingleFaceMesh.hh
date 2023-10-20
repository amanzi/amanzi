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

class SingleFaceMesh : public AmanziMesh::Mesh {
 public:
  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                 int f,
                 const AmanziGeometry::SurfaceCoordinateSystem& coordsys)
  {
    cacheSFM(mesh, f, coordsys);
    isSFM_ = true;
  }

  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int f)
  {
    AmanziGeometry::SurfaceCoordinateSystem coordsys(mesh->getFaceCentroid(f),
                                                     mesh->getFaceNormal(f));
    cacheSFM(mesh, f, coordsys);
    isSFM_ = true;
  }

  ~SingleFaceMesh(){};

  void cacheSFM(const Teuchos::RCP<const AmanziMesh::Mesh>& mc,
                int f,
                const AmanziGeometry::SurfaceCoordinateSystem& coordsys)
  {
    algorithms_ = mc->getAlgorithms();
    has_edges_ = mc->hasEdges();

    space_dim_ = mc->getSpaceDimension() - 1;
    manifold_dim_ = space_dim_;

    ncells_owned = 1;
    ncells_all = 1;
    nfaces_owned = 1;
    nfaces_all = 1;
    if (has_edges_) {
      cEntity_ID_View fn;
      mc->getFaceNodes(f, fn);
      nedges_owned = fn.size();
      nedges_all = nedges_owned;
    }
    AmanziMesh::cEntity_ID_View cfedges, cfnodes, cenodes;
    AmanziMesh::Entity_ID_View fedges, enodes, cells;
    mc->getFaceNodes(f, cfnodes);
    nnodes_owned = cfnodes.size();
    nnodes_all = nnodes_owned;
    cEntity_Direction_View fdirs;

    if (has_edges_) {
      mc->getFaceEdgesAndDirs(f, cfedges, &fdirs);
    } else {
      Entity_Direction_View tfdirs("", nnodes_owned);
      Kokkos::deep_copy(tfdirs, 1);
      fdirs = tfdirs;
    }


    // single surface cell
    Kokkos::resize(fedges, nnodes_owned);
    for (int i = 0; i < nnodes_owned; ++i) fedges[i] = i;

    data_.cell_faces.resize(1, nnodes_owned);
    Kokkos::deep_copy(data_.cell_faces.getRowUnmanaged<MemSpace_kind::HOST>(0), fedges);

    data_.cell_face_directions.resize(1, fdirs.size());
    Kokkos::deep_copy(data_.cell_face_directions.getRowUnmanaged<MemSpace_kind::HOST>(0), fdirs);

    const auto& xf = mc->getFaceCentroid(f);
    data_.cell_centroids.resize(1);
    data_.cell_centroids.h_view[0] = coordsys.Project(xf, true);

    double volume = mc->getFaceArea(f);
    data_.cell_volumes.resize(1);
    data_.cell_volumes.h_view[0] = volume;
    // cell nodes
    data_.cell_nodes.resize(1, nnodes_owned);
    data_.cell_coordinates.resize(1, nnodes_owned);
    data_.node_coordinates.resize(nnodes_owned);

    for (int i = 0; i < nnodes_owned; ++i) {
      data_.cell_nodes.getRowUnmanaged<MemSpace_kind::HOST>(0)[i] = i;
      auto xyz = mc->getNodeCoordinate(cfnodes[i]);
      auto res = coordsys.Project(xyz, true);
      data_.cell_coordinates.getRowUnmanaged<MemSpace_kind::HOST>(0)[i] = res;
      data_.node_coordinates.h_view[i] = res;
    }

    // cell faces
    data_.face_nodes.resize(nnodes_owned, 2);
    data_.face_areas.resize(nnodes_owned);
    data_.face_centroids.resize(nnodes_owned);
    std::vector<std::vector<AmanziGeometry::Point>> normals_v(nnodes_owned);
    for (int i = 0; i < nnodes_owned; ++i) {
      int n0(i), n1((i + 1) % nnodes_owned);
      auto cn1 = data_.cell_coordinates.getRowUnmanaged<MemSpace_kind::HOST>(0)[n1];
      auto cn0 = data_.cell_coordinates.getRowUnmanaged<MemSpace_kind::HOST>(0)[n0];

      data_.face_nodes.getRowUnmanaged<MemSpace_kind::HOST>(i)[0] = (fdirs[i] > 0) ? n0 : n1;
      data_.face_nodes.getRowUnmanaged<MemSpace_kind::HOST>(i)[1] = (fdirs[i] > 0) ? n1 : n0;

      auto tau = (cn1 - cn0) * fdirs[i];

      data_.face_areas.h_view[i] = norm(tau);
      normals_v[i].resize(1, AmanziGeometry::Point(tau[1], -tau[0]));
      data_.face_centroids.h_view[i] = (cn1 + cn0) / 2;
    }

    data_.face_normals = RaggedArray_DualView(normals_v);

    // cell edges
    data_.edge_vectors.resize(nnodes_owned);
    data_.edge_lengths.resize(nnodes_owned);
    Kokkos::deep_copy(data_.edge_lengths.h_view, data_.face_areas.h_view);

    for (int i = 0; i < nnodes_owned; ++i) {
      int n0(i), n1((i + 1) % nnodes_owned);
      auto cn1 = data_.cell_coordinates.getRowUnmanaged<MemSpace_kind::HOST>(0)[n1];
      auto cn0 = data_.cell_coordinates.getRowUnmanaged<MemSpace_kind::HOST>(0)[n0];

      data_.edge_vectors.h_view[i] = (cn1 + cn0);
    }

    // backward connectivity
    data_.face_cells.resize(nnodes_owned, 1);

    data_.face_edges.resize(nnodes_owned, 1);
    data_.face_edge_directions.resize(nnodes_owned, 1);

    for (int i = 0; i < nnodes_owned; ++i) {
      data_.face_cells.getRowUnmanaged<MemSpace_kind::HOST>(i)[0] = 0;
      data_.face_edge_directions.getRowUnmanaged<MemSpace_kind::HOST>(i)[0] = 1;
      data_.face_edges.getRowUnmanaged<MemSpace_kind::HOST>(i)[0] = i;
    }

    data_.edge_nodes.resize(nnodes_owned, 2);
    for (int i = 0; i < nnodes_owned; ++i) {
      data_.edge_nodes.getRowUnmanaged<MemSpace_kind::HOST>(i)[0] =
        data_.face_nodes.getRowUnmanaged<MemSpace_kind::HOST>(i)[0];
      data_.edge_nodes.getRowUnmanaged<MemSpace_kind::HOST>(i)[1] =
        data_.face_nodes.getRowUnmanaged<MemSpace_kind::HOST>(i)[1];
    }


    data_.cell_geometry_cached = true;
    data_.cell_faces_cached = true;
    data_.cell_edges_cached = true;
    data_.cell_nodes_cached = true;
    data_.cell_coordinates_cached = true;

    data_.face_geometry_cached = true;
    data_.face_cells_cached = true;
    data_.face_edges_cached = true;
    data_.face_nodes_cached = true;
    data_.face_coordinates_cached = true;

    data_.edge_geometry_cached = true;
    data_.edge_cells_cached = true;
    data_.edge_faces_cached = true;
    data_.edge_nodes_cached = true;
    data_.edge_coordinates_cached = true;
    data_.edge_lengths_cached = true;

    data_.node_cells_cached = true;
    data_.node_faces_cached = true;
    data_.node_edges_cached = true;
    data_.node_coordinates_cached = true;
    data_.parent_entities_cached = true;
    data_.cell_cellbelow_cached = true;
  }

 private:
  AmanziMesh::Entity_ID_View cell_node_ids_;
  Point_List cell_coords_;
  std::vector<AmanziMesh::Entity_ID_View> face_node_ids_;
};

} // namespace AmanziMesh
} // namespace Amanzi

#endif
