/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

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
#include "AmanziComm.hh"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

#include "Region.hh"
#include "MeshSurfaceCell.hh"

namespace Amanzi {
namespace AmanziMesh {

MeshSurfaceCell::MeshSurfaceCell(const Teuchos::RCP<const Mesh>& parent_mesh,
                                 const std::string& setname, bool flatten)
  : Mesh(getCommSelf(), parent_mesh->geometric_model(),
         parent_mesh->parameter_list(), true, false),
    parent_mesh_(parent_mesh)
{
  std::cout << "Constructor" << std::endl;

  // set dimensions
  if (flatten) {
    set_space_dimension(2);
  } else {
    set_space_dimension(3);
  }
  //    set_manifold_dimension(0); // ETC: this should be done, but it breaks
  //    overland flow
  set_manifold_dimension(2);

  // set my cells
  Kokkos::View<Entity_ID*> my_cells;
  parent_mesh->get_set_entities(setname,
                                AmanziMesh::Entity_kind::FACE,
                                AmanziMesh::Parallel_type::OWNED,
                                my_cells);
  AMANZI_ASSERT(my_cells.extent(0) == 1);
  parent_face_ = my_cells(0);

  // set my nodes
  Kokkos::View<Entity_ID*> my_nodes;
  parent_mesh->face_get_nodes(parent_face_, my_nodes);
  nodes_.resize(my_nodes.extent(0));
  if (flatten) {
    for (int i = 0; i != my_nodes.extent(0); ++i) {
      AmanziGeometry::Point parent_node;
      parent_mesh->node_get_coordinates(my_nodes(i), &parent_node);
      AmanziGeometry::Point child_node(2);
      child_node[0] = parent_node[0];
      child_node[1] = parent_node[1];
      nodes_[i] = child_node;
    }
  } else {
    for (int i = 0; i != my_nodes.extent(0); ++i) {
      parent_mesh->node_get_coordinates(my_nodes(i), &nodes_[i]);
    }
  }

  // set the maps
  cell_map_ = Teuchos::rcp(new Map_type(1, 1, 0, get_comm()));
  face_map_ = Teuchos::rcp(
    new Map_type((int)nodes_.size(), (int)nodes_.size(), 0, get_comm()));
  exterior_face_importer_ = Teuchos::rcp(new Import_type(face_map_, face_map_));

  // set the geometric model and sets
  auto gm = parent_mesh->geometric_model();

  init_cache();

  for (AmanziGeometry::GeometricModel::RegionConstIterator r =
         gm->RegionBegin();
       r != gm->RegionEnd();
       ++r) {
    // set to false as default
    sets_[(*r)->id()] = false;

    // determine if true
    if ((*r)->type() == AmanziGeometry::LABELEDSET ||
        (*r)->type() == AmanziGeometry::ENUMERATED) {
      // label pulled from parent
      Kokkos::View<Entity_ID*> faces_in_set;
      Kokkos::View<double*> vofs;
      parent_mesh->get_set_entities_and_vofs(
        (*r)->name(), FACE, Parallel_type::OWNED, faces_in_set, &vofs);
      sets_[(*r)->id()] = 0;
      for (int i = 0; i < faces_in_set.extent(0); ++i) {
        if (faces_in_set(i) == parent_face_) {
          sets_[(*r)->id()] = 1;
          break;
        }
      }
      // sets_[(*r)->id()] = std::find(faces_in_set.begin(), faces_in_set.end(),
      //        parent_face_) != faces_in_set.end();

    } else if ((*r)->is_geometric()) {
      // check containment
      if ((*r)->space_dimension() == 3) {
        sets_[(*r)->id()] =
          (*r)->inside(parent_mesh->face_centroid(parent_face_));

      } else if ((*r)->space_dimension() == 2 && flatten) {
        sets_[(*r)->id()] = (*r)->inside(cell_centroid(0));
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

unsigned int
MeshSurfaceCell::num_entities(const Entity_kind kind,
                              const Parallel_type ptype) const
{
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
void
MeshSurfaceCell::cell_get_nodes(const Entity_ID cellid,
                                Kokkos::View<Entity_ID*>& nodeids) const
{
  AMANZI_ASSERT(cellid == 0);
  // AMANZI_ASSERT(nodeids);
  Kokkos::resize(nodeids, nodes_.size());
  for (int i = 0; i != nodes_.size(); ++i) nodeids(i) = i;
}


// Get nodes of face
// On a distributed mesh, all nodes (OWNED or GHOST) of the face
// are returned
// In 3D, the nodes of the face are returned in ccw order consistent
// with the face normal
// In 2D, nfnodes is 2
void
MeshSurfaceCell::face_get_nodes(const Entity_ID faceid,
                                Kokkos::View<Entity_ID*>& nodeids) const
{
  AMANZI_ASSERT(faceid < nodes_.size());
  Kokkos::resize(nodeids, 2);
  nodeids(0) = faceid;
  nodeids(1) = (faceid + 1) % nodes_.size();
}


// Get nodes of edge
void
MeshSurfaceCell::edge_get_nodes(const Entity_ID edgeid, Entity_ID* nodeid0,
                                Entity_ID* nodeid1) const
{
  Errors::Message mesg("Not implemented");
  amanzi_throw(mesg);
}


// Upward adjacencies
//-------------------

// Cells of type 'ptype' connected to a node - The order of cells is
// not guaranteed to be the same for corresponding nodes on
// different processors
void
MeshSurfaceCell::node_get_cells(const Entity_ID nodeid,
                                const Parallel_type ptype,
                                Kokkos::View<Entity_ID*>& cellids) const
{
  Kokkos::resize(cellids, 1);
  cellids(0) = 0;
}


// Faces of type 'ptype' connected to a node - The order of faces is
// not guarnateed to be the same for corresponding nodes on
// different processors
void
MeshSurfaceCell::node_get_faces(const Entity_ID nodeid,
                                const Parallel_type ptype,
                                Kokkos::View<Entity_ID*>& faceids) const
{
  Errors::Message mesg("Not implemented");
  amanzi_throw(mesg);
}


// Get faces of ptype of a particular cell that are connected to the
// given node - The order of faces is not guarnateed to be the same
// for corresponding nodes on different processors
void
MeshSurfaceCell::node_get_cell_faces(const Entity_ID nodeid,
                                     const Entity_ID cellid,
                                     const Parallel_type ptype,
                                     Kokkos::View<Entity_ID*>& faceids) const
{
  Errors::Message mesg("Not implemented");
  amanzi_throw(mesg);
}

// Cells of type 'ptype' connected to an edges
void
MeshSurfaceCell::edge_get_cells(const Entity_ID edgeid,
                                const Parallel_type ptype,
                                Kokkos::View<Entity_ID*>& cellids) const
{
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
void
MeshSurfaceCell::cell_get_face_adj_cells(
  const Entity_ID cellid, const Parallel_type ptype,
  Kokkos::View<Entity_ID*>& fadj_cellids) const
{
  Kokkos::resize(fadj_cellids, 0);
}


// Node connected neighboring cells of given cell
// (a hex in a structured mesh has 26 node connected neighbors)
// The cells are returned in no particular order
void
MeshSurfaceCell::cell_get_node_adj_cells(
  const Entity_ID cellid, const Parallel_type ptype,
  Kokkos::View<Entity_ID*>& nadj_cellids) const
{
  Kokkos::resize(nadj_cellids, 0);
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
int
MeshSurfaceCell::deform(const std::vector<double>& target_cell_volumes_in,
                        const std::vector<double>& min_cell_volumes_in,
                        const Kokkos::View<Entity_ID*>& fixed_nodes,
                        const bool move_vertical)
{
  Errors::Message mesg("Not implemented");
  Exceptions::amanzi_throw(mesg);
  return -1;
}


//
// Mesh Sets for ICs, BCs, Material Properties and whatever else
//--------------------------------------------------------------
//

// Get number of entities of type 'category' in set
unsigned int
MeshSurfaceCell::get_set_size(const Set_ID setid, const Entity_kind kind,
                              const Parallel_type ptype) const
{
  if (sets_.at(setid)) { return kind == CELL ? 1 : nodes_.size(); }
  return 0;
}


unsigned int
MeshSurfaceCell::get_set_size(const std::string setname, const Entity_kind kind,
                              const Parallel_type ptype) const
{
  return get_set_size(
    geometric_model()->FindRegion(setname)->id(), kind, ptype);
}

// Get list of entities of type 'category' in set
void
MeshSurfaceCell::get_set_entities(const Set_ID setid, const Entity_kind kind,
                                  const Parallel_type ptype,
                                  Kokkos::View<Entity_ID*>& entids) const
{
  if (sets_.at(setid)) {
    if (kind == CELL) {
      Kokkos::resize(entids, 1);
    } else {
      Kokkos::resize(entids, nodes_.size());
      for (int i = 0; i != nodes_.size(); ++i) entids(i) = i;
    }
  } else {
    Kokkos::resize(entids, 0);
  }
}

void
MeshSurfaceCell::get_set_entities_and_vofs(const std::string setname,
                                           const Entity_kind kind,
                                           const Parallel_type ptype,
                                           Kokkos::View<Entity_ID*>& entids,
                                           Kokkos::View<double*>* vofs) const
{
  return get_set_entities(
    geometric_model()->FindRegion(setname)->id(), kind, ptype, entids);
}


// Miscellaneous functions
void
MeshSurfaceCell::write_to_exodus_file(const std::string filename) const
{
  Errors::Message mesg("Not implemented");
  Exceptions::amanzi_throw(mesg);
}


// get faces and face dirs of a cell. This can be called by
// cell_get_faces_and_dirs method of the base class and the data
// cached or it can be called directly by the
// cell_get_faces_and_dirs method of this class
void
MeshSurfaceCell::cell_get_faces_and_dirs_internal_(
  const Entity_ID cellid, Kokkos::View<Entity_ID*>& faceids,
  Kokkos::View<int*>& face_dirs) const
{
  AMANZI_ASSERT(cellid == 0);
  Kokkos::resize(faceids, nodes_.size());
  for (int i = 0; i != nodes_.size(); ++i) faceids(i) = i;
  Kokkos::resize(face_dirs, nodes_.size());
}


// Cells connected to a face - this function is implemented in each
// mesh framework. The results are cached in the base class
void
MeshSurfaceCell::face_get_cells_internal_(
  const Entity_ID faceid, const Parallel_type ptype,
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*>& cellids) const
{
  Kokkos::resize(cellids, 1);
  cellids(0) = 0;
}


// edges of a face - this function is implemented in each mesh
// framework. The results are cached in the base class
void
MeshSurfaceCell::face_get_edges_and_dirs_internal_(
  const Entity_ID faceid, Kokkos::View<Entity_ID*>& edgeids,
  Kokkos::View<int*>* edge_dirs, const bool ordered) const
{
  Errors::Message mesg("Not implemented");
  Exceptions::amanzi_throw(mesg);
}

// edges of a cell - this function is implemented in each mesh
// framework. The results are cached in the base class.
void
MeshSurfaceCell::cell_get_edges_internal_(
  const Entity_ID cellid, Kokkos::View<Entity_ID*>& edgeids) const
{
  Errors::Message mesg("Not implemented");
  Exceptions::amanzi_throw(mesg);
}


// edges and directions of a 2D cell - this function is implemented
// in each mesh framework. The results are cached in the base class.
void
MeshSurfaceCell::cell_2D_get_edges_and_dirs_internal_(
  const Entity_ID cellid, Kokkos::View<Entity_ID*>& edgeids,
  Kokkos::View<int*>* edge_dirs) const
{
  Errors::Message mesg("Not implemented");
  Exceptions::amanzi_throw(mesg);
}


} // namespace AmanziMesh
} // namespace Amanzi
