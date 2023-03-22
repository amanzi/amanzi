/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <algorithm>

#include <Teuchos_RCP.hpp>

#include "dbc.hh"
#include "errors.hh"
#include "MeshCache.hh"

#include "MeshFrameworkColumn.hh"

namespace Amanzi {
namespace AmanziMesh {

std::pair<double, AmanziGeometry::Point>
MeshFrameworkColumnAlgorithms::computeCellGeometry(const Mesh& mesh, const Entity_ID c) const
{
  return MeshAlgorithms::computeMeshColumnCellGeometry(mesh, c);
}

// -----------------------------------------------------------------------------
// Constructor: instantiates base mesh, generates new nodal coordinates,
//              fixes faces, and makes maps.
// -----------------------------------------------------------------------------
MeshFrameworkColumn::MeshFrameworkColumn(const Teuchos::RCP<MeshFramework>& col3D_mesh,
                       const Teuchos::RCP<Teuchos::ParameterList>& plist) :
  MeshFramework(col3D_mesh->getComm(), col3D_mesh->getGeometricModel(), plist),
  col3D_mesh_(col3D_mesh)
{
  AMANZI_ASSERT(col3D_mesh_->getSpaceDimension() == 3);
  AMANZI_ASSERT(col3D_mesh_->getManifoldDimension() == 3);

  // set supporting subclasses
  setSpaceDimension(3);
  setManifoldDimension(3);
  setVisMesh(col3D_mesh);

  // compute special geometric quantities for column entities (node
  // coordinates, face centroids, cell centroids, face areas)
  computeSpecialNodeCoordinates_();
}


// -----------------------------------------------------------------------------
// Compute special coordinates for the nodes - all the other
// quantities will follow suit
// -----------------------------------------------------------------------------
void MeshFrameworkColumn::computeSpecialNodeCoordinates_()
{
  // Assume that the column is vertical - nodes are stacked vertically
  // above each other. Assume that the base face is perfectly horizontal
  //
  // Base the X-Y coordinates of all the other nodes based on the X-Y
  // coordinates of the nodes of the bottom face. The Z-coordinate of
  // a node is the same as the Z-coordinate of the centroid of the
  // inter-cell (non-lateral) face they belong to. If these faces are
  // planar, then the volume of the cells will incur a second order
  // error due to "small" rotations in general cases. In 2D, when the
  // sides are vertical the error will be zero. In 3D, if the sides
  // are vertical and the face is symmetrical about the centroid, the
  // error will be zero.

  // Create a cached object, so that we can use columns.
  MeshCache<MemSpace_kind::HOST> col3D_mesh(col3D_mesh_, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null);
  col3D_mesh.buildColumns();

  // Get the ordered face indexes of the column
  auto colfaces = col3D_mesh.columns.getFaces<MemSpace_kind::HOST>(0);
  column_faces_ = colfaces;

  // mask for face index in the column of faces
  Kokkos::resize(face_in_column_,col3D_mesh.getNumEntities(Entity_kind::FACE,
          Parallel_kind::ALL));
  Kokkos::deep_copy(face_in_column_, -1); 

  // How many nodes each "horizontal" face has in the column
  cEntity_ID_View face_nodes;
  col3D_mesh.getFaceNodes(column_faces_[0], face_nodes);
  nfnodes_ = face_nodes.size();

  // Set up the new node coordinates This is done in two passes, which may be
  // unnecessary, but I'm not sure if face_centroid() would break if done in
  // one.
  int space_dim = getSpaceDimension(); // from parent mesh
  int nfaces = column_faces_.size();
  int nnodes = nfaces*nfnodes_;
  AmanziGeometry::Point p(space_dim);
  Point_View node_coordinates("node_coordinates", nnodes);
  Kokkos::deep_copy(node_coordinates, p); 

  for (int j=0; j!=nfaces; ++j) {
    // set the mask
    face_in_column_[column_faces_[j]] = j;

    // calculate node coordinates
    col3D_mesh.getFaceNodes(column_faces_[j], face_nodes);
    auto face_coordinates = col3D_mesh.getFaceCoordinates(column_faces_[j]);
    auto fcen = col3D_mesh.getFaceCentroid(column_faces_[j]);

    for (int i=0; i!=nfnodes_; ++i) {
      AmanziGeometry::Point coords(space_dim);

      // last coordinate is z-coordinate of face centroid
      coords[space_dim-1] = fcen[space_dim-1];

      // remain coordinates are coordinates of the corresponding node on
      // the bottom face
      for (int d=0; d!=space_dim-1; ++d) coords[d] = face_coordinates[i][d];
      node_coordinates[face_nodes[i]] = coords;
    }
  }

  // Set the mesh coordinates
  for (int n=0; n!=nnodes; ++n)
    setNodeCoordinate(n, node_coordinates[n]);
}


}  // namespace AmanziMesh
}  // namespace Amanzi
