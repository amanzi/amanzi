/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <algorithm>

#include <Teuchos_RCP.hpp>

#include "Column_mesh.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi
{

namespace AmanziMesh
{

Column_mesh::Column_mesh (Mesh& inmesh,
                          const int column_id, 
                          const VerboseObject *verbosity_obj) :
  parent_mesh_(inmesh),
  column_id_(column_id),
  column_cells_(inmesh.cells_of_column(column_id)),
  column_faces_(inmesh.faces_of_column(column_id)) {

  // sanity check
  
  ASSERT(column_id_ >= 0 && column_id_ < inmesh.num_columns());
  
  // Figure out how many nodes each "horizontal" face has in the
  // column
  
  Entity_ID_List fnodes;
  parent_mesh_.face_get_nodes(column_faces_[0],&fnodes);
  nfnodes_ = fnodes.size(); 
  
  // compute special geometric quantities for column entities (node
  // coordinates, face centroids, cell centroids, face areas)
  
  compute_special_node_coordinates_();
  
  // build epetra maps
  
  build_epetra_maps_();
  
}

// Parent of entity

Entity_ID Column_mesh::entity_get_parent(const Entity_kind kind,
                                         const Entity_ID entid) const {
  
  switch (kind) {
  case CELL:
    return column_cells_[entid];

  case FACE:
    return column_faces_[entid];

  case EDGE:
    return -1;

  case NODE: {
    int faceid = entid%nfnodes_;
    int parent_face = column_faces_[faceid];
    std::vector<Entity_ID> fnodes;
    parent_mesh_.face_get_nodes(parent_face,&fnodes);
    return fnodes[entid - faceid*nfnodes_];
  }
  default:
    amanzi_throw(Errors::Message("Unknown entity kind"));
  }
  
}

// Compute special coordinates for the nodes - all the other 
// quantities will follow suit

void Column_mesh::compute_special_node_coordinates_() {

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

  int spacedim = space_dimension(); // from parent mesh

  int nfaces = column_faces_.size();
  int nnodes = nfaces*nfnodes_;

  int parent_face = column_faces_[0];
  Entity_ID_List fnodes;
  parent_mesh_.face_get_coordinates(parent_face,&node_coordinates_);

  for (int j = 1; j < nfaces; ++j) {
    AmanziGeometry::Point fcen(spacedim);

    parent_mesh_.face_centroid(column_faces_[j],&fcen);
    
    for (int i = 0; i < nfnodes_; ++i) {
      AmanziGeometry::Point coords(spacedim);

      // node id

      int k = j*nfnodes_ + i;

      // last coordinate is z-coordinate of face centroid

      coords[spacedim-1] = fcen[spacedim-1];

      // remain coordinates are coordinates of the corresponding node on
      // the bottom face

      for (int d = 0; d < spacedim-1; ++d)
        coords[d] = (node_coordinates_[i])[d]; // node_coordinates_[i][d] ?
        
      node_coordinates_.push_back(coords);
    }
  }
  
}


// Build Epetra_maps indicating the distribution of entities across
// processors and their dependencies (through global IDs).
//
// In this case since the columns are all on one processor, the map is
// just a contiguous sequence of numbers and the communicator is a serial
// communicator


void Column_mesh::build_epetra_maps_() {
  Epetra_SerialComm epcomm;
  int indexBase = 0;

  int ncells = column_cells_.size();
  cell_map_ = new Epetra_Map(ncells,indexBase,epcomm);

  int nfaces = column_faces_.size();
  face_map_ = new Epetra_Map(nfaces,indexBase,epcomm);

  int nnodes = nfnodes_*nfaces;
  node_map_ = new Epetra_Map(nnodes,indexBase,epcomm);
}

} // close namespace AmanziMesh
} // close namespace Amanzi
