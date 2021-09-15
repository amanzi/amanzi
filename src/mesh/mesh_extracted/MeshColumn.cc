/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <algorithm>

#include <Teuchos_RCP.hpp>

#include "MeshColumn.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziMesh {

// -----------------------------------------------------------------------------
// Constructor: instantiates base mesh, generates new nodal coordinates,
//              fixes faces, and makes maps.
// -----------------------------------------------------------------------------
MeshColumn::MeshColumn(const Teuchos::RCP<Mesh>& col3D_mesh,
                       const Teuchos::RCP<const Teuchos::ParameterList>& plist) :
  Mesh(col3D_mesh->get_comm(), col3D_mesh->geometric_model(), plist, true, false),
  col3D_mesh_(col3D_mesh)
{
  AMANZI_ASSERT(col3D_mesh_->space_dimension() == 3);
  AMANZI_ASSERT(col3D_mesh_->manifold_dimension() == 3);

  // set supporting subclasses
  set_space_dimension(3);
  set_manifold_dimension(3);
  //  set_manifold_dimension(1); // ETC: this should be done, but it breaks overland flow

  // compute special geometric quantities for column entities (node
  // coordinates, face centroids, cell centroids, face areas)
  compute_special_node_coordinates_();

  // build epetra maps
  build_epetra_maps_();
}


MeshColumn::~MeshColumn() {
  if (face_map_) delete face_map_;
  if (exterior_face_map_) delete exterior_face_map_;
  if (exterior_face_importer_) delete exterior_face_importer_;
}



// Deform the mesh by moving given nodes to given coordinates
// If the flag keep_valid is true, then the nodes are moved
// only as much as possible without making the mesh invalid
// The final positions of the nodes is returned in final_positions
//
// Note, this is a bit wierd -- nodes are based on the 3D mesh, not this 1D
// mesh.  This should get fixed, but it would break physics code. --etc
int
MeshColumn::deform(const Entity_ID_List& nodeids,
                   const AmanziGeometry::Point_List& new_positions,
                   const bool keep_valid,
                   AmanziGeometry::Point_List *final_positions) {
  int ierr = col3D_mesh_->deform(nodeids, new_positions, keep_valid, final_positions);

  // recompute all geometric quantities
  compute_cell_geometric_quantities_();
  compute_face_geometric_quantities_();
  return ierr;
}

// Deform a mesh so that cell volumes conform as closely as possible
// to target volumes without dropping below the minimum volumes.  If
// move_vertical = true, nodes will be allowed to move only in the
// vertical direction (right now arbitrary node movement is not allowed)
// Nodes in any set in the fixed_sets will not be permitted to move.
int
MeshColumn::deform(const std::vector<double>& target_cell_volumes_in,
       const std::vector<double>& min_cell_volumes_in,
       const Entity_ID_List& fixed_nodes,
       const bool move_vertical) {
  int ierr = col3D_mesh_->deform(target_cell_volumes_in, min_cell_volumes_in,
          fixed_nodes, move_vertical);
  // recompute all geometric quantities
  compute_cell_geometric_quantities_();
  compute_face_geometric_quantities_();
  return ierr;
}


// -----------------------------------------------------------------------------
// Compute special coordinates for the nodes - all the other
// quantities will follow suit
// -----------------------------------------------------------------------------
void MeshColumn::compute_special_node_coordinates_() {
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

  // First build info about column topology
  col3D_mesh_->build_columns();

  // Get the ordered face indexes of the column
  const Entity_ID_List& colfaces = col3D_mesh_->faces_of_column(0);
  column_faces_ = colfaces;

  // mask for face index in the column of faces
  face_in_column_.resize(col3D_mesh_->num_entities(FACE, AmanziMesh::Parallel_type::ALL), -1);

  // How many nodes each "horizontal" face has in the column
  Entity_ID_List face_nodes;
  col3D_mesh_->face_get_nodes(column_faces_[0],&face_nodes);
  nfnodes_ = face_nodes.size();

  // Set up the new node coordinates This is done in two passes, which may be
  // unnecessary, but I'm not sure if face_centroid() would break if done in
  // one.
  int space_dim = space_dimension(); // from parent mesh
  int nfaces = column_faces_.size();
  int nnodes = nfaces*nfnodes_;
  AmanziGeometry::Point p(space_dim);
  std::vector<AmanziGeometry::Point> node_coordinates(nnodes, p);

  for (int j=0; j!=nfaces; ++j) {
    // set the mask
    face_in_column_[column_faces_[j]] = j;

    // calculate node coordinates
    std::vector<AmanziGeometry::Point> face_coordinates;
    col3D_mesh_->face_get_nodes(column_faces_[j], &face_nodes);
    col3D_mesh_->face_get_coordinates(column_faces_[j], &face_coordinates);

    AmanziGeometry::Point fcen = col3D_mesh_->face_centroid(column_faces_[j]);

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
    node_set_coordinates(n, node_coordinates[n]);
}


// -----------------------------------------------------------------------------
// Build Epetra_maps indicating the distribution of entities across
// processors and their dependencies (through global IDs).
//
// In this case since the columns are all on one processor, the map is
// just a contiguous sequence of numbers and the communicator is a serial
// communicator
// -----------------------------------------------------------------------------
void MeshColumn::build_epetra_maps_() {
  int indexBase = 0;

  int nfaces = column_faces_.size();
  face_map_ = new Epetra_Map(nfaces,indexBase,*get_comm());

  std::vector<int> ext_gids(2,-1);
  ext_gids[0] = 0;
  ext_gids[1] = nfaces-1;

  exterior_face_map_ = new Epetra_Map(-1, 2, &ext_gids[0], 0, *get_comm());
  exterior_face_importer_ = new Epetra_Import(*exterior_face_map_, *face_map_);
}

}  // namespace AmanziMesh
}  // namespace Amanzi
