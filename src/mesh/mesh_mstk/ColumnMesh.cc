/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <algorithm>

#include <Teuchos_RCP.hpp>

#include "ColumnMesh.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziMesh {


// -----------------------------------------------------------------------------
// Constructor: instantiates base MSTK mesh, generates new nodal coordiantes,
//              fixes faces, and makes maps.
// -----------------------------------------------------------------------------
ColumnMesh::ColumnMesh (const Mesh& inmesh,
                        const int column_id, 
                        const VerboseObject *verbosity_obj) :
    Mesh(verbosity_obj, true, false),
    parent_mesh_(inmesh),
    column_id_(column_id),
    extracted_(new Epetra_MpiComm(MPI_COMM_SELF),
              inmesh, inmesh.cells_of_column(column_id), CELL,
              false, false, true, false)
{
  ASSERT(column_id_ >= 0 && column_id_ < inmesh.num_columns());

  // set supporting subclasses
  set_comm(extracted_.get_comm());
  set_geometric_model(extracted_.geometric_model());
  
  // compute special geometric quantities for column entities (node
  // coordinates, face centroids, cell centroids, face areas)
  compute_special_node_coordinates_();
  
  // build epetra maps
  build_epetra_maps_();
}


ColumnMesh::~ColumnMesh () {
  if (face_map_) delete face_map_;
  if (exterior_face_map_) delete exterior_face_map_;
  if (exterior_face_importer_) delete exterior_face_importer_;  
}


// -----------------------------------------------------------------------------
// Compute special coordinates for the nodes - all the other 
// quantities will follow suit
// -----------------------------------------------------------------------------
void ColumnMesh::compute_special_node_coordinates_() {
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

  // Get the ordered face indexes of the column
  const Entity_ID_List& colfaces = extracted_.faces_of_column(0);
  column_faces_ = colfaces;

  // mask for face index in the column of faces
  face_in_column_.resize(extracted_.num_entities(FACE, AmanziMesh::USED), -1);
  
  // How many nodes each "horizontal" face has in the column
  Entity_ID_List face_nodes;
  extracted_.face_get_nodes(column_faces_[0],&face_nodes);
  nfnodes_ = face_nodes.size(); 

  // Set up the new node coordinates This is done in two passes, which may be
  // unnecessary, but I'm not sure if face_centroid() would break if done in
  // one.
  int spacedim = space_dimension(); // from parent mesh
  int nfaces = column_faces_.size();
  int nnodes = nfaces*nfnodes_;
  AmanziGeometry::Point p(spacedim);
  std::vector<AmanziGeometry::Point> node_coordinates(nnodes, p);
  
  for (int j=0; j!=nfaces; ++j) {
    // set the mask
    face_in_column_[column_faces_[j]] = j;

    // calculate node coordinates
    std::vector<AmanziGeometry::Point> face_coordinates;
    extracted_.face_get_nodes(column_faces_[j], &face_nodes);
    extracted_.face_get_coordinates(column_faces_[j], &face_coordinates);

    AmanziGeometry::Point fcen = extracted_.face_centroid(column_faces_[j]);
    
    for (int i=0; i!=nfnodes_; ++i) {
      AmanziGeometry::Point coords(spacedim);

      // last coordinate is z-coordinate of face centroid
      coords[spacedim-1] = fcen[spacedim-1];

      // remain coordinates are coordinates of the corresponding node on
      // the bottom face
      for (int d=0; d!=spacedim-1; ++d) coords[d] = face_coordinates[i][d];
        
      node_coordinates[face_nodes[i]] = coords;
    }
  }

  // Set the mesh coordinates
  for (int n=0; n!=nnodes; ++n)
    node_set_coordinates(n, node_coordinates[n]);
  
}


// Build Epetra_maps indicating the distribution of entities across
// processors and their dependencies (through global IDs).
//
// In this case since the columns are all on one processor, the map is
// just a contiguous sequence of numbers and the communicator is a serial
// communicator
void ColumnMesh::build_epetra_maps_() {
  Epetra_SerialComm epcomm;
  int indexBase = 0;

  int nfaces = column_faces_.size();
  face_map_ = new Epetra_Map(nfaces,indexBase,epcomm);

  std::vector<int> ext_gids(2,-1);
  ext_gids[0] = 0;
  ext_gids[1] = nfaces-1;

  exterior_face_map_ = new Epetra_Map(-1, 2, &ext_gids[0], 0, *get_comm());
  exterior_face_importer_ = new Epetra_Import(*exterior_face_map_, *face_map_);
  
}



} // close namespace AmanziMesh
} // close namespace Amanzi
