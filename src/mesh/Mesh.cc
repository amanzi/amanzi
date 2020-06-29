/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Rao Garimella
*/

//! Base mesh class for Amanzi

/*

Use the associated mesh factory to create an instance of a
derived class based on a particular mesh framework (like MSTK,
STKmesh etc.)

Design documentation:

This class is designed to be both flexible and performant (and somewhat
succeeds at both).  Many of the low-level methods here are accessed in the
innermost loop of performance-critical physics methods, so must themselves
be as simple as possible.  To meet this need, the low-level geometry ( such
as cell_volume() ) and relational topology ( such as cell_get_faces() ) are
ideally inlined to a single array access.

To accomplish this goal, this class stores a cache of data which is
accessed using the public interface.  This cache is (lazily) built up using
"cache-building" methods (i.e. compute_cell_geometry_() ) which are virtual,
but have a default implementation here.  These default implementations
simply call more virtual methods (i.e. cell_get_faces_internal() ) which
are NOT implemented here and instead use the various mesh frameworks.


There are a few exceptions to this -- Mesh_Simple and classes that inherit
from Mesh_Simple have no internal storage -- they use the cache directly as
their mesh representation.  They are not full-featured, but are useful for
some simple structured tests and mesh manipulations where a mesh is implied
(i.e. MeshLogical, MeshColumn, etc).

Note that not everything is cached, which would be memory overkill.  Only
things which have proven (or more accurately, a few were proven and the
rest were assumed based on those results) to be an issue have been cached.

NOTE: Lazy definition of the cache itself is necessarily "mutable".

*/

#include "AmanziMap.hh"
#include "AmanziVector.hh"
#include "Teuchos_CommHelpers.hpp"

#include "dbc.hh"

#include "Geometry.hh"
#include "RegionLabeledSet.hh"

#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

// Gather and cache cell to face and face to cell connectivity
// info. Both are gathered simultaneously because we assume that if we
// asked for faces we are likely to need both cell->face and
// face->cell info and in the rare case that we don't, we are ready to
// pay the storage price.
//
// Sets up: cell_face_ids_, cell_face_dirs_, face_cell_ids_, face_cell_ptype_

Mesh::Mesh(const Comm_ptr_type& comm,
           const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
           const Teuchos::RCP<const Teuchos::ParameterList>& plist,
           const bool request_faces, const bool request_edges)
  : comm_(comm),
    geometric_model_(gm),
    plist_(plist),
    space_dim_(-1),
    manifold_dim_(-1),
    mesh_type_(GENERAL),
    parent_(Teuchos::null),
    logical_(false),
    mesh_cache_(this,request_faces,request_edges)
{
  if (plist_ == Teuchos::null) {
    plist_ = Teuchos::rcp(new Teuchos::ParameterList("Mesh"));
  }
  vo_ = Teuchos::rcp(new VerboseObject(
      Keys::cleanPListName(plist_->name()), *plist_, comm_));
};

Entity_ID
Mesh::entity_get_parent_type(const Entity_kind kind, const Entity_ID entid) const
{
  Errors::Message mesg(
    "Parent/daughter entities not enabled in this framework.");
  Exceptions::amanzi_throw(mesg);
  return -1;
}

// NOTE: this is a on-processor routine
unsigned int
Mesh::cell_get_max_faces() const
{
  unsigned int n(0);
  int ncells = num_entities(CELL, Parallel_type::OWNED);
  for (int c = 0; c < ncells; ++c) { n = std::max(n, cell_get_num_faces(c)); }
  return n;
}


// NOTE: this is a on-processor routine
unsigned int
Mesh::cell_get_max_nodes() const
{
  unsigned int n(0);
  int ncells = num_entities(CELL, Parallel_type::OWNED);
  for (int c = 0; c < ncells; ++c) {
    Entity_ID_List nodes;
    cell_get_nodes(c, nodes);
    n = std::max(n, (unsigned int)nodes.size());
  }
  return n;
}


// NOTE: this is a on-processor routine
unsigned int
Mesh::cell_get_max_edges() const
{
  unsigned int n(0);
  if (mesh_cache_.edges_requested_) {
    int ncells = num_entities(CELL, Parallel_type::OWNED);
    for (int c = 0; c < ncells; ++c) {
      Kokkos::View<AmanziMesh::Entity_ID*> edges;
      cell_get_edges(c, edges);
      n = std::max(n, (unsigned int)edges.extent(0));
    }
  }
  return n;
}

// Get set ID given the name of the set - return 0 if no match is found
Set_ID
Mesh::set_id_from_name(const std::string setname) const
{
  if (!geometric_model_.get()) return 0;

  unsigned int ngr = geometric_model_->RegionSize();
  for (int i = 0; i < ngr; i++) {
    Teuchos::RCP<const AmanziGeometry::Region> rgn =
      geometric_model_->FindRegion(i);

    if (rgn->name() == setname) return rgn->id();
  }

  return 0;
}


// Get set name given the ID of the set - return 0 if no match is found
std::string
Mesh::set_name_from_id(const int setid) const
{
  std::string nullname("");

  if (!geometric_model_.get()) return nullname;

  unsigned int ngr = geometric_model_->RegionSize();
  for (int i = 0; i < ngr; i++) {
    Teuchos::RCP<const AmanziGeometry::Region> rgn =
      geometric_model_->FindRegion(i);

    if (rgn->id() == setid) return rgn->name();
  }

  return 0;
}


// Is there a set with this id and entity type
bool
Mesh::valid_set_id(Set_ID id, Entity_kind kind) const
{
  if (!geometric_model_.get()) return false;

  unsigned int ngr = geometric_model_->RegionSize();
  for (int i = 0; i < ngr; i++) {
    Teuchos::RCP<const AmanziGeometry::Region> rgn =
      geometric_model_->FindRegion(i);

    unsigned int rdim = rgn->manifold_dimension();

    if (rgn->id() == id) {
      // For regions of type Labeled Set and Color Function, the
      // dimension parameter is not guaranteed to be correct
      if (rgn->type() == AmanziGeometry::LABELEDSET ||
          rgn->type() == AmanziGeometry::COLORFUNCTION)
        return true;

      // If we are looking for a cell set the region has to be
      // of the same topological dimension as the cells
      if (kind == CELL && rdim == manifold_dim_) return true;

      // If we are looking for a side set, the region has to be
      // one topological dimension less than the cells
      if (kind == FACE && rdim == manifold_dim_ - 1) return true;

      // If we are looking for a node set, the region can be of any
      // dimension upto the spatial dimension of the domain
      if (kind == NODE) return true;
    }
  }

  return false;
}


// Is there a set with this name and entity type
bool
Mesh::valid_set_name(std::string name, Entity_kind kind) const
{
  if (!geometric_model_.get()) {
    Errors::Message mesg("Mesh sets not enabled because mesh was created "
                         "without reference to a geometric model");
    Exceptions::amanzi_throw(mesg);
  }

  Teuchos::RCP<const AmanziGeometry::Region> rgn;
  try {
    rgn = geometric_model_->FindRegion(name);
  } catch (...) {
    return false;
  }

  unsigned int rdim = rgn->manifold_dimension();

  // For regions of type Color Function, the dimension
  // parameter is not guaranteed to be correct
  if (rgn->type() == AmanziGeometry::COLORFUNCTION) return true;

  // For regions of type Labeled set, extract some more info and verify
  if (rgn->type() == AmanziGeometry::LABELEDSET) {
    auto lsrgn =
      Teuchos::rcp_dynamic_cast<const AmanziGeometry::RegionLabeledSet>(rgn);
    std::string entity_type = lsrgn->entity_str();

    if (parent() == Teuchos::null) {
      if ((kind == CELL && entity_type == "CELL") ||
          (kind == FACE && entity_type == "FACE") ||
          (kind == EDGE && entity_type == "EDGE") ||
          (kind == NODE && entity_type == "NODE"))
        return true;
    } else {
      if (kind == CELL && entity_type == "FACE") return true;
    }
    return false;
  }

  // If we are looking for a cell set the region has to be
  // of the same topological dimension as the cells or it
  // has to be a point region
  if (kind == CELL && (rdim >= manifold_dim_ || rdim == 1 || rdim == 0))
    return true;

  // If we are looking for a side set, the region has to be
  // one topological dimension less than the cells
  if (kind == FACE && (rdim >= manifold_dim_ - 1 || rdim == 0)) return true;

  // If we are looking for a node set, the region can be of any
  // dimension upto the spatial dimension of the domain
  if (kind == NODE) return true;

  return false;
}


bool
Mesh::point_in_cell(const AmanziGeometry::Point& p,
                    const Entity_ID cellid) const
{
  std::vector<AmanziGeometry::Point> ccoords;

  if (manifold_dim_ == 3) {
    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine
    int nf;
    Kokkos::View<Entity_ID*> faces;
    std::vector<unsigned int> nfnodes;
    Kokkos::View<int*> fdirs;
    std::vector<AmanziGeometry::Point> cfcoords;

    cell_get_faces_and_dirs(cellid, faces, fdirs);

    nf = faces.extent(0);
    nfnodes.resize(nf);

    for (int j = 0; j < nf; j++) {
      std::vector<AmanziGeometry::Point> fcoords;
      face_get_coordinates(faces(j), fcoords);
      nfnodes[j] = fcoords.size();
      if (fdirs(j) == 1) {
        for (int k = 0; k < nfnodes[j]; k++) {
          cfcoords.push_back(fcoords[k]); 
        }
      } else {
        for (int k = nfnodes[j] - 1; k >= 0; k--) {
          cfcoords.push_back(fcoords[k]); 
        }
      }
    }

    cell_get_coordinates(cellid, ccoords);

    return AmanziGeometry::point_in_polyhed(p, ccoords, nf, nfnodes, cfcoords);

  } else if (manifold_dim_ == 2) {
    cell_get_coordinates(cellid, ccoords);
    return AmanziGeometry::point_in_polygon(p, ccoords);
  }

  return false;
}


// synchronize node positions across processors
//
// on host
void
Mesh::update_ghost_node_coordinates()
{
  Exceptions::amanzi_throw("not yet implemented in tpetra");
  // int ndim = space_dimension();

  // Import_type importer(node_map(true), node_map(false));

  // // change last arg to false after debugging
  // MultiVector_type owned_node_coords(node_map(true), ndim, true);


  // AmanziGeometry::Point pnt(ndim);

  // // Fill the owned node coordinates
  // int nnodes_owned = num_entities(NODE, Parallel_type::OWNED);
  // for (int i = 0; i < nnodes_owned; i++) {
  //   node_get_coordinates(i,&pnt);
  //   for (int k = 0; k < ndim; k++)
  //     owned_node_coords[k][i] = pnt[k];
  // }

  // double **data;
  // owned_node_coords.ExtractView(&data);
  // MultiVector_type used_node_coords(View, *node_map(false), data, ndim);

  // used_node_coords.Import(owned_node_coords, importer, Insert);

  // int nnodes_used = num_entities(NODE,Parallel_type::ALL);
  // for (int i = nnodes_owned; i < nnodes_used; i++) {
  //   double xyz[3];
  //   for (int k = 0; k < ndim; k++)
  //     xyz[k] = used_node_coords[k][i];
  //   pnt.set(xyz);
  //   node_set_coordinates(i, pnt);
  // }
}


// Deform the mesh according to a given set of new node positions
// If keep_valid is true, the routine will cut back node displacement
// if the cells connected to a moved node become invalid
//
// CAVEAT: this is not parallel, and so all deformations must be consistent
// across ghost entities, and provided for ghost nodes.  User beware!
int
Mesh::deform(const Kokkos::View<Entity_ID*>& nodeids,
             const Kokkos::View<AmanziGeometry::Point*>& new_positions)
{
  AMANZI_ASSERT(nodeids.extent(0) == new_positions.extent(0));

  int nn = nodeids.extent(0);
  for (int j = 0; j != nn; ++j) {
    node_set_coordinates(nodeids(j), new_positions(j));
  }

  // recompute all geometric quantities
  mesh_cache_.compute_cell_geometric_quantities_();
  if (mesh_cache_.faces_requested_) mesh_cache_.compute_face_geometric_quantities_();
  if (mesh_cache_.edges_requested_) mesh_cache_.compute_edge_geometric_quantities_();

  int nc = num_entities(CELL, Parallel_type::ALL);
  for (int c = 0; c != nc; ++c) {
    if (cell_volume(c, false) < 0.) return 0;
  }
  return 1;
}


// Deform the mesh according to a given set of new node positions
// If keep_valid is true, the routine will cut back node displacement
// if the cells connected to a moved node become invalid
//
// CAVEAT: this is not parallel, and so all deformations must be consistent
// across ghost entities, and provided for ghost nodes.  User beware!
int
Mesh::deform(const Entity_ID_List& nodeids,
             const AmanziGeometry::Point_List& new_positions,
             const bool keep_valid,
             std::vector<AmanziGeometry::Point>& final_positions)
{
  int status = 1;

  AMANZI_ASSERT(nodeids.size() == new_positions.size());

  // Once we start moving nodes around, the precomputed/cached
  // geometric quantities are no longer valid. So any geometric calls
  // must use the "recompute=true" option until the end of this routine
  // where we once again call compute_geometric_quantities
  int nn = nodeids.size();
  final_positions.resize(nn); 

  bool done_outer = false;
  int iter = 0, maxiter = 5;

  while (!done_outer) {
    double totdisp2 = 0.0;

    for (int j = 0; j < nn; j++) {
      Entity_ID node = nodeids[j];

      AmanziGeometry::Point oldcoords, newcoords, dispvec;
      Entity_ID_List cells;

      node_get_coordinates(node, &oldcoords);
      dispvec = new_positions[j] - oldcoords;

      node_get_cells(node, Parallel_type::ALL, cells);
      int nc = cells.size();

      double mult = 1.0;
      bool done = false;
      bool allvalid = true;

      while (!done) {
        newcoords = oldcoords + mult * dispvec;

        node_set_coordinates(node, newcoords);

        if (keep_valid) { // check if the cells remain valid
          allvalid = true;
          for (int k = 0; k < nc; k++)
            if (cell_volume(cells[k], true) < 0.0) {
              allvalid = false;
              break;
            }
        }

        if (allvalid)
          done = true;
        else {
          if (mult < 1.0e-5) done = true;
          mult = mult / 2.0;
        }
      } // while (!done)

      if (!allvalid) { // could not move the node even a bit
        status = 0;    // perhaps the mesh was invalid to start?
        node_set_coordinates(node, oldcoords);
      } else {
        AmanziGeometry::Point actual_dispvec = newcoords - oldcoords;
        totdisp2 += L22(actual_dispvec);
      }
    } // for (j = 0; j < nn; j++)

    if (totdisp2 < 1.0e-12) done_outer = 1;

    if (++iter == maxiter) break;
  } // while (!done_outer)


  for (int j = 0; j < nn; j++) {
    Entity_ID node = nodeids[j];

    AmanziGeometry::Point newcoords;

    node_get_coordinates(node, &newcoords);
    final_positions[j] = newcoords;
  }

  // recompute all geometric quantities
  mesh_cache_.compute_cell_geometric_quantities_();
  if (mesh_cache_.faces_requested_) mesh_cache_.compute_face_geometric_quantities_();
  if (mesh_cache_.edges_requested_) mesh_cache_.compute_edge_geometric_quantities_();

  return status;
}

std::string
Mesh::cell_type_to_name(const Cell_type type)
{
  switch (type) {
  case TRI:
    return "triangle";
  case QUAD:
    return "quad";
  case POLYGON:
    return "polygon";
  case TET:
    return "tetrahedron";
  case PYRAMID:
    return "pyramid";
  case PRISM:
    return "prism";
  case HEX:
    return "hexahedron";
  case POLYHED:
    return "polyhedron";
  default:
    return "unknown";
  }
}


void
Mesh::PrintMeshStatistics() const
{
  if (vo_.get() && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    int ncells =
      num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    int nfaces =
      num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

    int min_out[2], max_out[2], sum_out[2], tmp_in[2] = { ncells, nfaces };
    Teuchos::reduceAll(*get_comm(), Teuchos::REDUCE_MIN, 2, tmp_in, min_out);
    Teuchos::reduceAll(*get_comm(), Teuchos::REDUCE_MAX, 2, tmp_in, max_out);
    Teuchos::reduceAll(*get_comm(), Teuchos::REDUCE_SUM, 2, tmp_in, sum_out);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "cells, tot/min/max: " << sum_out[0] << "/" << min_out[0]
               << "/" << max_out[0] << "\n";
    *vo_->os() << "faces, tot/min/max: " << sum_out[1] << "/" << min_out[1]
               << "/" << max_out[1] << "\n\n";
  }
}

AmanziGeometry::Point
get_coordinate(const Mesh& m, const Entity_kind& kind, Entity_ID id)
{
  switch(kind) {
    case(CELL): return m.cell_centroid(id);
    case(FACE): return m.face_centroid(id);
    case(BOUNDARY_FACE): m.face_centroid(m.face_map(false)->getLocalElement(m.exterior_face_map(false)->getGlobalElement(id)));
    case(EDGE): return AmanziGeometry::Point();
    case(NODE):
      { AmanziGeometry::Point p; m.node_get_coordinates(id, &p); return p; }
    default: return AmanziGeometry::Point();
  }
}


} // namespace AmanziMesh
} // namespace Amanzi
