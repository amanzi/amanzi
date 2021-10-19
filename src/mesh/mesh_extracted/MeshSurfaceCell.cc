/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//
// This is a mesh for a single surface cell.
//
// This exists solely because we need "surface meshes" extracted from
// MeshColumn.  This is really just 1 cell.  Really.
//

#include <vector>
#include <string>
#include <algorithm>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Epetra_SerialComm.h"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

#include "Region.hh"
#include "MeshSurfaceCell.hh"

namespace Amanzi {
namespace AmanziMesh {

MeshSurfaceCell::MeshSurfaceCell(const Teuchos::RCP<const Mesh>& parent_mesh,
        const std::string& setname, bool flatten)
    : Mesh(getCommSelf(), parent_mesh->geometric_model(), parent_mesh->parameter_list(), true, false)
{
  parent_ = parent_mesh;

  // set dimensions
  if (flatten) {
    setSpaceDimension(2);
  } else {
    setSpaceDimension(3);
  }
  //    setManifoldDimension(0); // ETC: this should be done, but it breaks overland flow
  setManifoldDimension(2);

  // set my cells
  Entity_ID_List my_cells;
  parent_->get_set_entities(setname, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED, &my_cells);
  AMANZI_ASSERT(my_cells.size() == 1);
  parent_face_ = my_cells[0];

  // set my nodes
  Entity_ID_List my_nodes;
  parent_->face_get_nodes(parent_face_, &my_nodes);
  nodes_.resize(my_nodes.size());
  if (flatten) {
    for (int i=0; i!=my_nodes.size(); ++i) {
      AmanziGeometry::Point parent_node;
      parent_->node_get_coordinates(my_nodes[i], &parent_node);
      AmanziGeometry::Point child_node(2);
      child_node[0] = parent_node[0];
      child_node[1] = parent_node[1];
      nodes_[i] = child_node;
    }
  } else {
    for (int i=0; i!=my_nodes.size(); ++i) {
      parent_->node_get_coordinates(my_nodes[i], &nodes_[i]);
    }
  }

  // set the maps
  cell_map_ = Teuchos::rcp(new Epetra_Map(1, 0, *get_comm()));
  face_map_ = Teuchos::rcp(new Epetra_Map((int)nodes_.size(), 0, *get_comm()));
  exterior_face_importer_ =
      Teuchos::rcp(new Epetra_Import(*face_map_,*face_map_));

  // set the geometric model and sets
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm = parent_->geometric_model();

  for (auto& r : *gm) {
    // set to false as default
    sets_[r->id()] = false;

    // determine if true
    if ((r->type() == AmanziGeometry::LABELEDSET ||
         r->type() == AmanziGeometry::ENUMERATED) &&
        parent_->valid_set_name(r->name(), FACE)) {
      // label pulled from parent
      Entity_ID_List faces_in_set;
      std::vector<double> vofs;
      parent_->get_set_entities_and_vofs(r->name(), FACE, Parallel_type::OWNED, &faces_in_set, &vofs);
      sets_[r->id()] = std::find(faces_in_set.begin(), faces_in_set.end(),
              parent_face_) != faces_in_set.end();

    } else if (r->is_geometric()) {
      // check containment
      if (r->space_dimension() == 3) {
        sets_[r->id()] = r->inside(parent_->face_centroid(parent_face_));

      } else if (r->space_dimension() == 2 && flatten) {
        sets_[r->id()] = r->inside(cell_centroid(0));
      }
    }
  }

  // set the cell type
  if (nodes_.size() == 3) {
    cell_type_ = TRI;
  } else if (nodes_.size() == 4) {
    cell_type_ = QUAD;
  } else {
    cell_type_ = POLYGON;
  }
}


//
// General mesh information
// -------------------------
//

// Number of entities of any kind (cell, face, node) and in a
// particular category (OWNED, GHOST, ALL)
unsigned int MeshSurfaceCell::num_entities(const Entity_kind kind,
        const Parallel_type ptype) const {
  int count;
  switch (kind) {
    case CELL:
      count = 1;
      break;

    default: // num_nodes == num_faces == num_boundary_faces
      count = nodes_.size();
      break;
  }
  return count;
}



//
// Mesh Entity Adjacencies
//-------------------------

// Downward Adjacencies
//---------------------
// Get nodes of a cell
void MeshSurfaceCell::cell_get_nodes(const Entity_ID cellid,
        Entity_ID_List *nodeids) const {
  AMANZI_ASSERT(cellid == 0);
  AMANZI_ASSERT(nodeids);
  nodeids->resize(nodes_.size());
  for (int i=0; i!=nodes_.size(); ++i) (*nodeids)[i] = i;
}


// Get nodes of face
// On a distributed mesh, all nodes (OWNED or GHOST) of the face
// are returned
// In 3D, the nodes of the face are returned in ccw order consistent
// with the face normal
// In 2D, nfnodes is 2
void MeshSurfaceCell::face_get_nodes(const Entity_ID faceid,
        Entity_ID_List *nodeids) const {
  AMANZI_ASSERT(faceid < nodes_.size());
  nodeids->resize(2);
  (*nodeids)[0] = faceid;
  (*nodeids)[1] = (faceid + 1) % nodes_.size();
}


// Get nodes of edge
void MeshSurfaceCell::edge_get_nodes(const Entity_ID edgeid,
        Entity_ID *nodeid0, Entity_ID *nodeid1) const {
  Errors::Message mesg("Not implemented");
  amanzi_throw(mesg);
}


// Upward adjacencies
//-------------------

// Cells of type 'ptype' connected to a node - The order of cells is
// not guaranteed to be the same for corresponding nodes on
// different processors
void MeshSurfaceCell::node_get_cells(const Entity_ID nodeid,
        const Parallel_type ptype,
        Entity_ID_List *cellids) const {
  cellids->resize(1);
  (*cellids)[0] = 0;
}


// Faces of type 'ptype' connected to a node - The order of faces is
// not guarnateed to be the same for corresponding nodes on
// different processors
void MeshSurfaceCell::node_get_faces(const Entity_ID nodeid,
        const Parallel_type ptype,
        Entity_ID_List *faceids) const {
  Errors::Message mesg("Not implemented");
  amanzi_throw(mesg);
}


// Cells of type 'ptype' connected to an edges
void MeshSurfaceCell::edge_get_cells(const Entity_ID edgeid,
        const Parallel_type ptype,
        Entity_ID_List *cellids) const {
  Errors::Message mesg("Not implemented");
  amanzi_throw(mesg);
}


// Same level adjacencies
//-----------------------

// Face connected neighboring cells of given cell of a particular ptype
// (e.g. a hex has 6 face neighbors)

// The order in which the cellids are returned cannot be
// guaranteed in general except when ptype = ALL, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces
void MeshSurfaceCell::cell_get_face_adj_cells(const Entity_ID cellid,
        const Parallel_type ptype,
        Entity_ID_List *fadj_cellids) const {
  fadj_cellids->resize(0);
}


//
// Mesh entity geometry
//--------------
//


// Deform a mesh so that cell volumes conform as closely as possible
// to target volumes without dropping below the minimum volumes.  If
// move_vertical = true, nodes will be allowed to move only in the
// vertical direction (right now arbitrary node movement is not allowed)
// Nodes in any set in the fixed_sets will not be permitted to move.
int MeshSurfaceCell::deform(const std::vector<double>& target_cell_volumes_in,
                            const std::vector<double>& min_cell_volumes_in,
                            const Entity_ID_List& fixed_nodes,
                            const bool move_vertical) {
  Errors::Message mesg("Not implemented");
  Exceptions::amanzi_throw(mesg);
  return -1;
}



//
// Mesh Sets for ICs, BCs, Material Properties and whatever else
//--------------------------------------------------------------
//
bool MeshSurfaceCell::valid_set_type(const AmanziGeometry::RegionType rtype,
        const Entity_kind kind) const {
  if (rtype == AmanziGeometry::LABELEDSET && kind == CELL) return true;
  if (rtype == AmanziGeometry::ENUMERATED && kind == CELL) return true;
  if (rtype == AmanziGeometry::BOX) return true;
  if (rtype == AmanziGeometry::PLANE) return true;
  if (rtype == AmanziGeometry::POINT) return true;
  return false;
}

// Get number of entities of type 'category' in set
unsigned int
MeshSurfaceCell::get_set_size(const std::string& setname,
        const Entity_kind kind,
        const Parallel_type ptype) const {
  auto setid = geometric_model()->FindRegion(setname)->id();
  if (sets_.at(setid)) {
    return kind == CELL ? 1 : nodes_.size();
  }
  return 0;
}

// Get list of entities of type 'category' in set
void MeshSurfaceCell::get_set_entities_and_vofs(const std::string& setname,
        const Entity_kind kind,
        const Parallel_type ptype,
        Entity_ID_List *entids,
        std::vector<double> *vofs) const {
  auto setid = geometric_model()->FindRegion(setname)->id();
  if (sets_.at(setid)) {
    if (kind == CELL) {
      entids->resize(1,0);
    } else {
      entids->resize(nodes_.size());
      for (int i=0; i!=nodes_.size(); ++i) (*entids)[i] = i;
    }
  } else {
    entids->resize(0);
  }
}


// Miscellaneous functions
void MeshSurfaceCell::write_to_exodus_file(const std::string filename) const {
  Errors::Message mesg("Not implemented");
  Exceptions::amanzi_throw(mesg);
}


// get faces and face dirs of a cell. This can be called by
// cell_get_faces_and_dirs method of the base class and the data
// cached or it can be called directly by the
// cell_get_faces_and_dirs method of this class
void MeshSurfaceCell::cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
        Entity_ID_List *faceids,
        std::vector<int> *face_dirs,
        const bool ordered) const {
  AMANZI_ASSERT(cellid == 0);
  faceids->resize(nodes_.size());
  for (int i=0; i!=nodes_.size(); ++i) (*faceids)[i] = i;
  face_dirs->resize(nodes_.size(),1);
}


// Cells connected to a face - this function is implemented in each
// mesh framework. The results are cached in the base class
void MeshSurfaceCell::face_get_cells_internal_(const Entity_ID faceid,
        const Parallel_type ptype,
        Entity_ID_List *cellids) const {
  cellids->resize(1,0);
}


// edges of a face - this function is implemented in each mesh
// framework. The results are cached in the base class
void MeshSurfaceCell::face_get_edges_and_dirs_internal_(const Entity_ID faceid,
        Entity_ID_List *edgeids,
        std::vector<int> *edge_dirs,
        const bool ordered) const {
  Errors::Message mesg("Not implemented");
  Exceptions::amanzi_throw(mesg);
}

// edges of a cell - this function is implemented in each mesh
// framework. The results are cached in the base class.
void MeshSurfaceCell::cell_get_edges_internal_(const Entity_ID cellid,
        Entity_ID_List *edgeids) const {
  Errors::Message mesg("Not implemented");
  Exceptions::amanzi_throw(mesg);
}

} // close namespace AmanziMesh
} // close namespace Amanzi




