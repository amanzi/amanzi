//
// Mesh
//
// Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
// Amanzi is released under the three-clause BSD License.
// The terms of use and "as is" disclaimer for this license are
// provided in the top-level COPYRIGHT file.
//
// Base mesh class for Amanzi
//
// Use the associated mesh factory to create an instance of a
// derived class based on a particular mesh framework (like MSTK, MOAB, etc)
//
// Design documentation:
//
// This class is designed to be both flexible and performant (and somewhat
// succeeds at both).  Many of the low-level methods here are accessed in the
// innermost loop of performance-critical physics methods, so must themselves
// be as simple as possible.  To meet this need, the low-level geometry ( such
// as cell_volume() ) and relational topology ( such as cell_get_faces() ) are
// ideally inlined to a single array access.
//
// To accomplish this goal, this class stores a cache of data which is
// accessed using the public interface.  This cache is (lazily) built up using
// "cache-building" methods (i.e. compute_cell_geometry_() ) which are virtual,
// but have a default implementation here.  These default implementations
// simply call more virtual methods (i.e. cell_get_faces_internal() ) which
// are NOT implemented here and instead use the various mesh frameworks.
//
//
// There are a few exceptions to this -- Mesh_Simple and classes that inherit
// from Mesh_Simple have no internal storage -- they use the cache directly as
// their mesh representation.  They are not full-featured, but are useful for
// some simple structured tests and mesh manipulations where a mesh is implied
// (i.e. MeshLogical, MeshColumn, etc).
//
// Note that not everything is cached, which would be memory overkill.  Only
// things which have proven (or more accurately, a few were proven and the
// rest were assumed based on those results) to be an issue have been cached.
//
// NOTE: Lazy definition of the cache itself is necessarily "mutable".
//

#include "Epetra_MultiVector.h"

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
           const bool request_faces,
           const bool request_edges)
    : comm_(comm),
      geometric_model_(gm),
      plist_(plist),
      faces_requested_(request_faces),
      edges_requested_(request_edges),
      space_dim_(-1),
      manifold_dim_(-1),
      mesh_type_(GENERAL),
      cell_geometry_precomputed_(false),
      face_geometry_precomputed_(false),
      edge_geometry_precomputed_(false),
      columns_built_(false),
      cell2face_info_cached_(false),
      face2cell_info_cached_(false),
      cell2edge_info_cached_(false),
      face2edge_info_cached_(false),
      parent_(Teuchos::null),
      logical_(false),
      kdtree_faces_initialized_(false) {
  if (plist_ == Teuchos::null) {
    plist_ = Teuchos::rcp(new Teuchos::ParameterList("Mesh"));
  }
  vo_ = Teuchos::rcp(new VerboseObject(comm_,
          Keys::cleanPListName(plist_->name()), *plist_));
};


void
Mesh::cache_cell_face_info_() const
{
  int ncells = num_entities(CELL, Parallel_type::ALL);
  cell_face_ids_.resize(ncells);
  cell_face_dirs_.resize(ncells);

  int nfaces = num_entities(FACE, Parallel_type::ALL);
  face_cell_ids_.resize(nfaces);
  face_cell_ptype_.resize(nfaces);

  for (int c = 0; c < ncells; c++) {
    cell_get_faces_and_dirs_internal_(c, &(cell_face_ids_[c]),
            &(cell_face_dirs_[c]), false);

    int nf = cell_face_ids_[c].size();
    for (int jf = 0; jf < nf; jf++) {
      Entity_ID f = cell_face_ids_[c][jf];
      int dir = cell_face_dirs_[c][jf];
      face_cell_ids_[f].push_back(dir > 0 ? c : ~c); // store 1s complement of c if dir is -ve
      face_cell_ptype_[f].push_back(entity_get_ptype(CELL, c));
    }
  }

  cell2face_info_cached_ = true;
  face2cell_info_cached_ = true;
  faces_requested_ = true;
}


// Gather and cache face to edge connectivity info.
void
Mesh::cache_face2edge_info_() const
{
  int nfaces = num_entities(FACE,Parallel_type::ALL);
  face_edge_ids_.resize(nfaces);
  face_edge_dirs_.resize(nfaces);

  for (int f = 0; f < nfaces; f++) {
    Entity_ID_List fedgeids;
    std::vector<int> fedgedirs;

    face_get_edges_and_dirs_internal_(f, &(face_edge_ids_[f]),
            &(face_edge_dirs_[f]), true);
  }

  face2edge_info_cached_ = true;
  faces_requested_ = true;
  edges_requested_ = true;
}


// Gather and cache cell to edge connectivity info.
void
Mesh::cache_cell2edge_info_() const
{
  int ncells = num_entities(CELL,Parallel_type::ALL);
  cell_edge_ids_.resize(ncells);

  if (space_dim_ == 2) {
    cell_2D_edge_dirs_.resize(ncells);
    for (int c = 0; c < ncells; c++)
      cell_2D_get_edges_and_dirs_internal_(c, &(cell_edge_ids_[c]),
              &(cell_2D_edge_dirs_[c]));
  }
  else
    for (int c = 0; c < ncells; c++)
      cell_get_edges_internal_(c, &(cell_edge_ids_[c]));

  cell2edge_info_cached_ = true;
}


Entity_ID
Mesh::entity_get_parent(const Entity_kind kind, const Entity_ID entid) const
{
  Errors::Message mesg("Parent/daughter entities not enabled in this framework.");
  Exceptions::amanzi_throw(mesg);
  return -1;
}


unsigned int
Mesh::cell_get_num_faces(const Entity_ID cellid) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2face_info_cached_) cache_cell_face_info_();
  return cell_face_ids_[cellid].size();

#else  // Non-cached version
  Entity_ID_List cfaceids;
  std::vector<int> cfacedirs;

  cell_get_faces_and_dirs_internal_(cellid, &cfaceids, &cfacedirs, false);
  return cfaceids.size();
#endif
}


// NOTE: this is a on-processor routine
unsigned int
Mesh::cell_get_max_faces() const
{
  unsigned int n(0);
  int ncells = num_entities(CELL, Parallel_type::OWNED);
  for (int c = 0; c < ncells; ++c) {
    n = std::max(n, cell_get_num_faces(c));
  }
  return n;
}


// NOTE: this is a on-processor routine
unsigned int
Mesh::cell_get_max_nodes() const
{
  unsigned int n(0);
  int ncells = num_entities(CELL, Parallel_type::OWNED);
  for (int c = 0; c < ncells; ++c) {
    AmanziMesh::Entity_ID_List nodes;
    cell_get_nodes(c, &nodes);
    n = std::max(n, (unsigned int) nodes.size());
  }
  return n;
}


// NOTE: this is a on-processor routine
unsigned int
Mesh::cell_get_max_edges() const
{
  unsigned int n(0);
  if (edges_requested_) {
    int ncells = num_entities(CELL, Parallel_type::OWNED);
    for (int c = 0; c < ncells; ++c) {
      AmanziMesh::Entity_ID_List edges;
      cell_get_edges(c, &edges);
      n = std::max(n, (unsigned int) edges.size());
    }
  }
  return n;
}


void
Mesh::cell_get_faces_and_dirs(const Entity_ID cellid,
                              Entity_ID_List *faceids,
                              std::vector<int> *face_dirs,
                              const bool ordered) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2face_info_cached_) cache_cell_face_info_();

  if (ordered)
    cell_get_faces_and_dirs_internal_(cellid, faceids, face_dirs, ordered);
  else {
    Entity_ID_List &cfaceids = cell_face_ids_[cellid];
    *faceids = cfaceids; // copy operation

    if (face_dirs) {
      std::vector<int> &cfacedirs = cell_face_dirs_[cellid];
      *face_dirs = cfacedirs; // copy operation
    }
  }

#else // Non-cached version
  cell_get_faces_and_dirs_internal_(cellid, faceids, face_dirs, ordered);

#endif
}


// Get the bisectors, i.e. vectors from cell centroid to face centroids.
void Mesh::cell_get_faces_and_bisectors(
    const Entity_ID cellid,
    Entity_ID_List *faceids,
    std::vector<AmanziGeometry::Point> *bisectors,
    const bool ordered) const
{
  cell_get_faces(cellid, faceids, ordered);

  AmanziGeometry::Point cc = cell_centroid(cellid);
  if (bisectors) {
    bisectors->resize(faceids->size());
    for (int i=0; i!=faceids->size(); ++i) {
      (*bisectors)[i] = face_centroid((*faceids)[i]) - cc;
    }
  }
  return;
}


// Cells connected to a face - cache the results the first time it
// is called and then return the cached results subsequently
void Mesh::face_get_cells(const Entity_ID faceid, const Parallel_type ptype,
                          Entity_ID_List *cellids) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!face2cell_info_cached_) cache_cell_face_info_();

  cellids->clear();
  int n = face_cell_ptype_[faceid].size();

  for (int i = 0; i < n; i++) {
    Parallel_type cell_ptype = face_cell_ptype_[faceid][i];
    if (cell_ptype == Parallel_type::PTYPE_UNKNOWN) continue;

    int c = face_cell_ids_[faceid][i];
    if (std::signbit(c)) c = ~c;  // strip dir info by taking 1s complement

    if (ptype == Parallel_type::ALL || ptype == cell_ptype)
      cellids->push_back(c);
  } 

#else  // Non-cached version
  Entity_ID_List fcells;
  face_get_cells_internal_(faceid, Parallel_type::ALL, &fcells);

  cellids->clear();
  int n = face_cell_ptype_[faceid].size();

  switch (ptype) {
    case Parallel_type::ALL:
      for (int i = 0; i < n; i++)
        if (entity_get_ptype(CELL,fcells[i]) != PTYPE_UNKNOWN)
          cellids->push_back(fcells[i]);
      break;
    case Parallel_type::OWNED:
      for (int i = 0; i < n; i++)
        if (entity_get_ptype(CELL,fcells[i]) == Parallel_type::OWNED)
          cellids->push_back(fcells[i]);
      break;
    case Parallel_type::GHOST:
      for (int i = 0; i < n; i++)
        if (entity_get_ptype(CELL,fcells[i]) == Parallel_type::GHOST)
          cellids->push_back(fcells[i]);
      break;
    default:
      break;
  }

#endif
}


void
Mesh::face_get_edges_and_dirs(const Entity_ID faceid,
                              Entity_ID_List *edgeids,
                              std::vector<int> *edge_dirs,
                              const bool ordered) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!face2edge_info_cached_) cache_face2edge_info_();

  *edgeids = face_edge_ids_[faceid]; // copy operation

  if (edge_dirs) {
    std::vector<int> &fedgedirs = face_edge_dirs_[faceid];
    *edge_dirs = fedgedirs; // copy operation
  }

#else  // Non-cached version
  face_get_edges_and_dirs_internal_(faceid, edgeids, edge_dirs, ordered);

#endif
}


// Get the local ID of a face edge in a cell edge list
void
Mesh::face_to_cell_edge_map(const Entity_ID faceid,
                            const Entity_ID cellid,
                            std::vector<int> *map) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!face2edge_info_cached_) cache_face2edge_info_();
  if (!cell2edge_info_cached_) cache_cell2edge_info_();

  map->resize(face_edge_ids_[faceid].size());
  for (int f = 0; f < face_edge_ids_[faceid].size(); ++f) {
    Entity_ID fedge = face_edge_ids_[faceid][f];

    for (int c = 0; c < cell_edge_ids_[cellid].size(); ++c) {
      if (fedge == cell_edge_ids_[cellid][c]) {
        (*map)[f] = c;
        break;
      }
    }
  }

#else // non-cached version

  Entity_ID_List fedgeids, cedgeids;
  std::vector<int> fedgedirs;

  face_get_edges_and_dirs(faceid, &fedgeids, &fedgedirs, true);
  cell_get_edges(cellid, &cedgeids);

  map->resize(fedgeids.size(),-1);
  for (int f = 0; f < fedgeids.size(); ++f) {
    Entity_ID fedge = fedgeids[f];

    for (int c = 0; c < cedgeids.size(); ++c) {
      if (fedge == cedgeids[c]) {
        (*map)[f] = c;
        break;
      }
    }
  }

#endif
}


void
Mesh::cell_get_edges(const Entity_ID cellid,
                     Entity_ID_List *edgeids) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2edge_info_cached_) cache_cell2edge_info_();

  Entity_ID_List &cedgeids = cell_edge_ids_[cellid];
  *edgeids = cell_edge_ids_[cellid]; // copy operation

#else  // Non-cached version
  cell_get_edges_internal_(cellid, edgeids);

#endif
}


void
Mesh::cell_2D_get_edges_and_dirs(const Entity_ID cellid,
                                 Entity_ID_List *edgeids,
                                 std::vector<int> *edgedirs) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2edge_info_cached_) cache_cell2edge_info_();

  *edgeids = cell_edge_ids_[cellid]; // copy operation
  *edgedirs = cell_2D_edge_dirs_[cellid];

#else  // Non-cached version
  cell_2D_get_edges_and_dirs_internal_(cellid, edgeids, edgedirs);

#endif
}


int
Mesh::compute_cell_geometric_quantities_() const
{
  int ncells = num_entities(CELL,Parallel_type::ALL);

  cell_volumes_.resize(ncells);
  cell_centroids_.resize(ncells);
  for (int i = 0; i < ncells; i++) {
    double volume;
    AmanziGeometry::Point centroid(space_dim_);

    compute_cell_geometry_(i,&volume,&centroid);

    cell_volumes_[i] = volume;
    cell_centroids_[i] = centroid;
  }

  cell_geometry_precomputed_ = true;

  return 1;
}


int
Mesh::compute_face_geometric_quantities_() const
{
  if (space_dimension() == 3 && manifold_dimension() == 2) {
    // need cell centroids to compute normals
    if (!cell_geometry_precomputed_)
      compute_cell_geometric_quantities_();
  }

  int nfaces = num_entities(FACE,Parallel_type::ALL);

  face_areas_.resize(nfaces);
  face_centroids_.resize(nfaces);
  face_normals_.resize(nfaces);

  for (int i = 0; i < nfaces; i++) {
    double area;
    AmanziGeometry::Point centroid(space_dim_);
    std::vector<AmanziGeometry::Point> normals;

    // normal0 and normal1 are outward normals of the face with
    // respect to the cell0 and cell1 of the face. The natural normal
    // of the face points out of cell0 and into cell1. If one of these
    // cells do not exist, then the normal is the null vector.
    compute_face_geometry_(i, &area, &centroid, &normals);

    face_areas_[i] = area;
    face_centroids_[i] = centroid;
    face_normals_[i] = normals;
  }

  face_geometry_precomputed_ = true;

  return 1;
}


int
Mesh::compute_edge_geometric_quantities_() const
{
  int nedges = num_entities(EDGE,Parallel_type::ALL);

  edge_vectors_.resize(nedges);
  edge_lengths_.resize(nedges);

  for (int i = 0; i < nedges; i++) {
    double length;
    AmanziGeometry::Point evector(space_dim_);

    compute_edge_geometry_(i,&length,&evector);

    edge_lengths_[i] = length;
    edge_vectors_[i] = evector;
  }

  edge_geometry_precomputed_ = true;

  return 1;
}


int
Mesh::compute_cell_geometry_(const Entity_ID cellid, double *volume,
                             AmanziGeometry::Point *centroid) const
{
  if (manifold_dim_ == 3) {
    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine

    // General polyhedra always need to have an explicit face
    // representation - special elements like hexes can get away
    // without (but we have yet to put in the code for the standard
    // node ordering and computation for these special elements)
    Entity_ID_List faces;
    std::vector<unsigned int> nfnodes;
    std::vector<int> fdirs;
    std::vector<AmanziGeometry::Point> ccoords, cfcoords, fcoords;

    cell_get_faces_and_dirs(cellid,&faces,&fdirs);

    int nf = faces.size();
    nfnodes.resize(nf);

    if (nf == 2) {
      /* special case of column mesh - only top and bottom faces
         are returned */
      AmanziGeometry::Point fcentroid0(space_dim_), fcentroid1(space_dim_);
      AmanziGeometry::Point normal(space_dim_);
      double farea;

      /* compute volume on the assumption that the top and bottom faces form
         a vertical columnar cell or in other words a polygonal prism */

      face_get_coordinates(faces[0],&fcoords);
      AmanziGeometry::polygon_get_area_centroid_normal(fcoords,&farea,
              &fcentroid0,&normal);

      face_get_coordinates(faces[1],&fcoords);
      AmanziGeometry::polygon_get_area_centroid_normal(fcoords,&farea,
              &fcentroid1,&normal);

      *centroid = (fcentroid0+fcentroid1)/2.0;
      double height = norm(fcentroid1-fcentroid0);

      *volume = farea*height;
    }
    else { /* general case */

      for (int j = 0; j < nf; j++) {
        face_get_coordinates(faces[j],&fcoords);
        nfnodes[j] = fcoords.size();

        if (fdirs[j] == 1) {
          for (int k = 0; k < nfnodes[j]; k++)
            cfcoords.push_back(fcoords[k]);
        }
        else {
          for (int k = nfnodes[j]-1; k >=0; k--)
            cfcoords.push_back(fcoords[k]);
        }
      }

      cell_get_coordinates(cellid,&ccoords);

      AmanziGeometry::polyhed_get_vol_centroid(ccoords,nf,nfnodes,
                                               cfcoords,volume,
                                               centroid);
    }
    return 1;
  }
  else if (manifold_dim_ == 2) {
    std::vector<AmanziGeometry::Point> ccoords;
    cell_get_coordinates(cellid,&ccoords);

    AmanziGeometry::Point normal(space_dim_);

    AmanziGeometry::polygon_get_area_centroid_normal(ccoords,volume,centroid,
                                                     &normal);
    return 1;
  }

  return 0;
}


int
Mesh::compute_face_geometry_(const Entity_ID faceid, double *area,
                             AmanziGeometry::Point *centroid,
                             std::vector<AmanziGeometry::Point> *normals) const
{
  AmanziGeometry::Point_List fcoords;
  normals->clear();

  if (manifold_dim_ == 3) {
    // 3D Elements with possibly curved faces

    face_get_coordinates(faceid,&fcoords);

    AmanziGeometry::Point normal(3);
    AmanziGeometry::polygon_get_area_centroid_normal(fcoords,area,centroid,&normal);

    Entity_ID_List cellids;
    face_get_cells(faceid, Parallel_type::ALL, &cellids);
    AMANZI_ASSERT(cellids.size() <= 2);

    normals->resize(cellids.size(), AmanziGeometry::Point(0.0, 0.0, 0.0));
    for (int i = 0; i < cellids.size(); i++) {
      Entity_ID_List cellfaceids;
      std::vector<int> cellfacedirs;
      int dir = 1;

      cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);

      bool found = false;
      for (int j = 0; j < cellfaceids.size(); j++) {
        if (cellfaceids[j] == faceid) {
          found = true;
          dir = cellfacedirs[j];
          break;
        }
      }

      AMANZI_ASSERT(found);

      (*normals)[i] = (dir == 1) ? normal : -normal;
    }

    return 1;
  }
  else if (manifold_dim_ == 2) {
    if (space_dim_ == 2) {   // 2D mesh
      face_get_coordinates(faceid,&fcoords);

      AmanziGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);

      *centroid = 0.5*(fcoords[0]+fcoords[1]);

      AmanziGeometry::Point normal(evec[1],-evec[0]);

      Entity_ID_List cellids;
      face_get_cells(faceid, Parallel_type::ALL, &cellids);
      AMANZI_ASSERT(cellids.size() <= 2);

      normals->resize(cellids.size(), AmanziGeometry::Point(0.0, 0.0));
      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        std::vector<int> cellfacedirs;
        int dir = 1;

        cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);

        bool found = false;
        for (int j = 0; j < cellfaceids.size(); j++) {
          if (cellfaceids[j] == faceid) {
            found = true;
            dir = cellfacedirs[j];
            break;
          }
        }

        AMANZI_ASSERT(found);

        (*normals)[i] = (dir == 1) ? normal : -normal;
      }

      return 1;
    }
    else {  // Surface mesh - cells are 2D, coordinates are 3D

      // Since the edge likely forms a discontinuity in the surface
      // (or may even be the intersection of several surfaces), we
      // have to compute an outward normal to the edge with respect to
      // each face

      face_get_coordinates(faceid,&fcoords);

      AmanziGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);

      *centroid = 0.5*(fcoords[0]+fcoords[1]);

      Entity_ID_List cellids;
      face_get_cells(faceid, Parallel_type::ALL, &cellids);

      normals->resize(cellids.size(), AmanziGeometry::Point(0.0, 0.0, 0.0));
      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        std::vector<int> cellfacedirs;
        int dir = 1;

        cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);

        bool found = false;
        for (int j = 0; j < cellfaceids.size(); j++) {
          if (cellfaceids[j] == faceid) {
            found = true;
            dir = cellfacedirs[j];
            break;
          }
        }

        AMANZI_ASSERT(found);

        AmanziGeometry::Point cvec = fcoords[0]-cell_centroids_[cellids[i]];
        AmanziGeometry::Point trinormal = cvec^evec;

        AmanziGeometry::Point normal = evec^trinormal;

        double len = norm(normal);
        normal /= len;
        normal *= *area;

        (*normals)[i] = normal;  // Always an outward normal as calculated
      }

      return 1;
    }
  }
  return 0;
}


int
Mesh::compute_edge_geometry_(const Entity_ID edgeid, double *edge_length,
                             AmanziGeometry::Point *edge_vector) const
{
  (*edge_vector).set(0.0L);
  *edge_length = 0.0;

  Entity_ID node0, node1;

  edge_get_nodes(edgeid,&node0,&node1);

  AmanziGeometry::Point point0, point1;
  node_get_coordinates(node0,&point0);
  node_get_coordinates(node1,&point1);

  *edge_vector = point1 - point0;
  *edge_length = norm(*edge_vector);

  return 0;
}


// Volume/Area of cell
double
Mesh::cell_volume(const Entity_ID cellid, const bool recompute) const
{
  if (!cell_geometry_precomputed_) {
    compute_cell_geometric_quantities_();
    return cell_volumes_[cellid];
  }
  else {
    if (recompute) {
      double volume;
      AmanziGeometry::Point centroid(space_dim_);
      compute_cell_geometry_(cellid, &volume, &centroid);
      return volume;
    }
    else
      return cell_volumes_[cellid];
  }
}


// Area/length of face
double Mesh::face_area(const Entity_ID faceid, const bool recompute) const
{
  AMANZI_ASSERT(faces_requested_);

  if (!face_geometry_precomputed_) {
    compute_face_geometric_quantities_();
    return face_areas_[faceid];
  }
  else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(space_dim_);
      std::vector<AmanziGeometry::Point> normals;
      compute_face_geometry_(faceid, &area, &centroid, &normals);
      return area;
    }
    else
      return face_areas_[faceid];
  }
}


// Length of an edge
double
Mesh::edge_length(const Entity_ID edgeid, const bool recompute) const
{
  AMANZI_ASSERT(edges_requested_);

  if (!edge_geometry_precomputed_) {
    compute_edge_geometric_quantities_();
    return edge_lengths_[edgeid];
  }
  else {
    if (recompute) {
      double length;
      AmanziGeometry::Point vector(space_dim_);
      compute_edge_geometry_(edgeid, &length, &vector);
      return length;
    }
    else
      return edge_lengths_[edgeid];
  }
}


// Centroid of cell
AmanziGeometry::Point
Mesh::cell_centroid(const Entity_ID cellid,
                    const bool recompute) const
{
  if (!cell_geometry_precomputed_) {
    compute_cell_geometric_quantities_();
    return cell_centroids_[cellid];
  }
  else {
    if (recompute) {
      double volume;
      AmanziGeometry::Point centroid(space_dim_);
      compute_cell_geometry_(cellid, &volume, &centroid);
      return centroid;
    }
    else
      return cell_centroids_[cellid];
  }
}


// Centroid of face
AmanziGeometry::Point
Mesh::face_centroid(const Entity_ID faceid, const bool recompute) const
{
  AMANZI_ASSERT(faces_requested_);

  if (!face_geometry_precomputed_) {
    compute_face_geometric_quantities_();
    return face_centroids_[faceid];
  }
  else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(space_dim_);
      std::vector<AmanziGeometry::Point> normals;
      compute_face_geometry_(faceid, &area, &centroid, &normals);
      return centroid;
    }
    else
      return face_centroids_[faceid];
  }
}


// Centroid of edge
AmanziGeometry::Point
Mesh::edge_centroid(const Entity_ID edgeid) const
{
  Entity_ID p0, p1;
  AmanziGeometry::Point xyz0, xyz1;

  edge_get_nodes(edgeid, &p0, &p1);
  node_get_coordinates(p0, &xyz0);
  node_get_coordinates(p1, &xyz1);
  return (xyz0+xyz1)/2;
}


// Normal to face
// The vector is normalized and then weighted by the area of the face
//
// If recompute is TRUE, then the normal is recalculated using current
// face coordinates but not stored. (If the recomputed normal must be
// stored, then call recompute_geometric_quantities).
//
// If cellid is not specified, the normal is the natural normal of the face
// If cellid is specified, the normal is the outward normal with respect
// to the cell. In planar and solid meshes, the normal with respect to
// the cell on one side of the face is just the negative of the normal
// with respect to the cell on the other side. In general surfaces meshes,
// this will not be true at C1 discontinuities
//
// if cellid is specified, then orientation returns the direction of
// the natural normal of the face with respect to the cell (1 is
// pointing out of the cell and -1 pointing in)
AmanziGeometry::Point
Mesh::face_normal(const Entity_ID faceid,
                  const bool recompute,
                  const Entity_ID cellid,
                  int *orientation) const
{
  AMANZI_ASSERT(faces_requested_);

  std::vector<AmanziGeometry::Point> *fnormals = nullptr;
  std::vector<AmanziGeometry::Point> fnormals_new;

  if (!face_geometry_precomputed_) {
    compute_face_geometric_quantities_();

    fnormals = &(face_normals_[faceid]);
  }
  else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(space_dim_);
      compute_face_geometry_(faceid, &area, &centroid, &fnormals_new);
      fnormals = &fnormals_new;
    }
    else
      fnormals = &(face_normals_[faceid]);
  }

  AMANZI_ASSERT(fnormals->size() > 0);

  if (cellid == -1) {
    // Return the natural normal. This is the normal with respect to
    // the first cell, appropriately adjusted according to whether the
    // face is pointing into the cell (-ve cell id) or out

    int c = face_cell_ids_[faceid][0];
    return std::signbit(c) ? -(*fnormals)[0] : (*fnormals)[0];
  } else {
    // Find the index of 'cellid' in list of cells connected to face

    int dir;
    int irefcell;
    int nfc = face_cell_ids_[faceid].size();
    for (irefcell = 0; irefcell < nfc; irefcell++) {
      int c = face_cell_ids_[faceid][irefcell];
      if (c == cellid || ~c == cellid) {
        dir = std::signbit(c) ? -1 : 1;
        break;
      }
    }
    AMANZI_ASSERT(irefcell < nfc);
    if (orientation) *orientation = dir;  // if orientation was requested
    return (*fnormals)[irefcell];
  }
}


// Direction vector of edge
AmanziGeometry::Point
Mesh::edge_vector(const Entity_ID edgeid,
                  const bool recompute,
                  const Entity_ID pointid,
                  int *orientation) const
{
  AMANZI_ASSERT(edges_requested_);

  AmanziGeometry::Point evector(space_dim_);
  AmanziGeometry::Point& evector_ref = evector; // to avoid extra copying

  if (!edge_geometry_precomputed_)
    compute_edge_geometric_quantities_();

  if (recompute) {
    double length;
    compute_edge_geometry_(edgeid, &length, &evector);
    // evector_ref already points to evector
  }
  else
    evector_ref = edge_vectors_[edgeid];

  if (orientation) *orientation = 1;

  if (pointid == -1)
    return evector_ref;
  else {
    Entity_ID p0, p1;
    edge_get_nodes(edgeid, &p0, &p1);

    if (pointid == p0)
      return evector_ref;
    else {
      if (orientation) *orientation=-1;
      return -evector_ref;
    }
  }
}


// Get set ID given the name of the set - return 0 if no match is found
Set_ID
Mesh::set_id_from_name(const std::string setname) const
{
  if (!geometric_model_.get()) return 0;

  unsigned int ngr = geometric_model_->size();
  for (int i = 0; i < ngr; i++) {
    Teuchos::RCP<const AmanziGeometry::Region> rgn = geometric_model_->FindRegion(i);

    if (rgn->name() == setname)
      return rgn->id();
  }

  return 0;
}


// Get set name given the ID of the set - return 0 if no match is found
std::string
Mesh::set_name_from_id(const int setid) const
{
  std::string nullname("");

  if (!geometric_model_.get()) return nullname;

  unsigned int ngr = geometric_model_->size();
  for (int i = 0; i < ngr; i++) {
    Teuchos::RCP<const AmanziGeometry::Region> rgn = geometric_model_->FindRegion(i);

    if (rgn->id() == setid)
      return rgn->name();
  }

  return 0;
}


// Is there a set with this id and entity type
bool
Mesh::valid_set_id(Set_ID id, Entity_kind kind) const
{
  if (!geometric_model_.get()) return false;

  unsigned int ngr = geometric_model_->size();
  for (int i = 0; i < ngr; i++) {
    Teuchos::RCP<const AmanziGeometry::Region> rgn = geometric_model_->FindRegion(i);

    unsigned int rdim = rgn->manifold_dimension();

    if (rgn->id() == id) {
      // For regions of type Labeled Set and Color Function, the
      // dimension parameter is not guaranteed to be correct
      if (rgn->type() == AmanziGeometry::LABELEDSET ||
          rgn->type() == AmanziGeometry::COLORFUNCTION) return true;

      // If we are looking for a cell set the region has to be
      // of the same topological dimension as the cells
      if (kind == CELL && rdim == manifold_dim_) return true;

      // If we are looking for a side set, the region has to be
      // one topological dimension less than the cells
      if (kind == FACE && rdim == manifold_dim_-1) return true;

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
    Errors::Message mesg("Mesh sets not enabled because mesh was created without reference to a geometric model");
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
    auto lsrgn = Teuchos::rcp_dynamic_cast<const AmanziGeometry::RegionLabeledSet>(rgn);
    std::string entity_type = lsrgn->entity_str();

    if (parent() == Teuchos::null) {
      if ((kind == CELL && entity_type == "CELL") ||
          (kind == FACE && entity_type == "FACE") ||
          (kind == EDGE && entity_type == "EDGE") ||
          (kind == NODE && entity_type == "NODE"))
        return true;
    } else {
      if (kind == CELL && entity_type == "FACE") return true;
      if (kind == FACE && entity_type == "EDGE") return true;
    } 
    return false;
  }

  // If we are looking for a cell set the region has to be
  // of the same topological dimension as the cells or it
  // has to be a point region
  if (kind == CELL && (rdim >= manifold_dim_ || rdim == 1 || rdim == 0)) return true;

  // If we are looking for a side set, the region has to be
  // one topological dimension less than the cells
  if (kind == FACE && (rdim >= manifold_dim_-1 || rdim == 0)) return true;

  // If we are looking for a node set, the region can be of any
  // dimension upto the spatial dimension of the domain
  if (kind == NODE) return true;

  return false;
}


bool
Mesh::point_in_cell(const AmanziGeometry::Point &p, const Entity_ID cellid) const
{
  std::vector<AmanziGeometry::Point> ccoords;

  if (manifold_dim_ == 3) {
    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine
    int nf;
    Entity_ID_List faces;
    std::vector<unsigned int> nfnodes;
    std::vector<int> fdirs;
    std::vector<AmanziGeometry::Point> cfcoords;

    cell_get_faces_and_dirs(cellid,&faces,&fdirs);

    nf = faces.size();
    nfnodes.resize(nf);

    for (int j = 0; j < nf; j++) {
      std::vector<AmanziGeometry::Point> fcoords;

      face_get_coordinates(faces[j],&fcoords);
      nfnodes[j] = fcoords.size();

      if (fdirs[j] == 1) {
        for (int k = 0; k < nfnodes[j]; k++)
          cfcoords.push_back(fcoords[k]);
      }
      else {
        for (int k = nfnodes[j]-1; k >=0; k--)
          cfcoords.push_back(fcoords[k]);
      }
    }

    cell_get_coordinates(cellid,&ccoords);

    return AmanziGeometry::point_in_polyhed(p,ccoords,nf,nfnodes,cfcoords);

  }
  else if (manifold_dim_ == 2) {
    cell_get_coordinates(cellid,&ccoords);
    return AmanziGeometry::point_in_polygon(p,ccoords);
  }

  return false;
}


// synchronize node positions across processors
void
Mesh::update_ghost_node_coordinates()
{
  int ndim = space_dimension();

  Epetra_Map owned_node_map = node_map(false);
  Epetra_Map used_node_map = node_map(true);
  Epetra_Import importer(used_node_map, owned_node_map);

  // change last arg to false after debugging
  Epetra_MultiVector owned_node_coords(node_map(true),ndim,true);

  AmanziGeometry::Point pnt(ndim);

  // Fill the owned node coordinates
  int nnodes_owned = num_entities(NODE,Parallel_type::OWNED);
  for (int i = 0; i < nnodes_owned; i++) {
    node_get_coordinates(i,&pnt);
    for (int k = 0; k < ndim; k++)
      owned_node_coords[k][i] = pnt[k];
  }

  double **data;
  owned_node_coords.ExtractView(&data);
  Epetra_MultiVector used_node_coords(View, owned_node_map, data, ndim);

  used_node_coords.Import(owned_node_coords, importer, Insert);

  int nnodes_used = num_entities(NODE,Parallel_type::ALL);
  for (int i = nnodes_owned; i < nnodes_used; i++) {
    double xyz[3];
    for (int k = 0; k < ndim; k++)
      xyz[k] = used_node_coords[k][i];
    pnt.set(xyz);
    node_set_coordinates(i, pnt);
  }
}


// Deform the mesh according to a given set of new node positions
// If keep_valid is true, the routine will cut back node displacement
// if the cells connected to a moved node become invalid
//
// CAVEAT: this is not parallel, and so all deformations must be consistent
// across ghost entities, and provided for ghost nodes.  User beware!
int
Mesh::deform(const Entity_ID_List& nodeids,
             const AmanziGeometry::Point_List& new_positions) {
  AMANZI_ASSERT(nodeids.size() == new_positions.size());

  int nn = nodeids.size();
  for (int j=0; j!=nn; ++j) {
    node_set_coordinates(nodeids[j], new_positions[j]);
  }

  // recompute all geometric quantities
  compute_cell_geometric_quantities_();
  if (faces_requested_) compute_face_geometric_quantities_();
  if (edges_requested_) compute_edge_geometric_quantities_();

  int nc = num_entities(CELL,Parallel_type::ALL);
  for (int c=0; c!=nc; ++c) {
    if (cell_volume(c) < 0.) return 0;
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
             AmanziGeometry::Point_List *final_positions)
{
  int status = 1;

  AMANZI_ASSERT(nodeids.size() == new_positions.size());
  AMANZI_ASSERT(final_positions != NULL);

  // Once we start moving nodes around, the precomputed/cached
  // geometric quantities are no longer valid. So any geometric calls
  // must use the "recompute=true" option until the end of this routine
  // where we once again call compute_geometric_quantities
  int nn = nodeids.size();
  final_positions->resize(nn);

  bool done_outer = false;
  int iter = 0, maxiter = 5;

  while (!done_outer) {
    double totdisp2 = 0.0;

    for (int j = 0; j < nn; j++) {
      Entity_ID node = nodeids[j];

      AmanziGeometry::Point oldcoords, newcoords, dispvec;
      Entity_ID_List cells;

      node_get_coordinates(node,&oldcoords);
      dispvec = new_positions[j]-oldcoords;

      node_get_cells(node,Parallel_type::ALL,&cells);
      int nc = cells.size();

      double mult = 1.0;
      bool done = false;
      bool allvalid = true;

      while (!done) {
        newcoords = oldcoords + mult*dispvec;

        node_set_coordinates(node,newcoords);

        if (keep_valid) { // check if the cells remain valid
          allvalid = true;
          for (int k = 0; k < nc; k++)
            if (cell_volume(cells[k],true) < 0.0) {
              allvalid = false;
              break;
            }
        }

        if (allvalid)
          done = true;
        else {
          if (mult < 1.0e-5)
            done = true;
          mult = mult/2.0;
        }
      } // while (!done)

      if (!allvalid) { // could not move the node even a bit
        status = 0;    // perhaps the mesh was invalid to start?
        node_set_coordinates(node,oldcoords);
      }
      else {
        AmanziGeometry::Point actual_dispvec = newcoords - oldcoords;
        totdisp2 += L22(actual_dispvec);
      }
    } // for (j = 0; j < nn; j++)

    if (totdisp2 < 1.0e-12)
      done_outer = 1;

    if (++iter == maxiter)
      break;
  } // while (!done_outer)


  for (int j = 0; j < nn; j++) {
    Entity_ID node = nodeids[j];

    AmanziGeometry::Point newcoords;

    node_get_coordinates(node,&newcoords);
    (*final_positions)[j] = newcoords;
  }

  // recompute all geometric quantities
  compute_cell_geometric_quantities_();
  if (faces_requested_) compute_face_geometric_quantities_();
  if (edges_requested_) compute_edge_geometric_quantities_();

  return status;
}


// Figure out columns of cells in a semi-structured mesh and cache the
// information for later.
//
// The columns are defined by identifying all boundary faces in a provided
// set, then collecting cells and faces upward.
//
// NOTE: Currently ghost columns are built too because we didn't know that
// they weren't necessary. --etc
//
// NOTE: this could be const thanks to only changing mutable things, but we
// choose not to be since it is directly a user asking for provided column
// structure.
int
Mesh::build_columns(const std::string& setname) const
{

  if (columns_built_) return 0;

  // Allocate space and initialize.
  int nn = num_entities(NODE,Parallel_type::ALL);
  int nf = num_entities(FACE,Parallel_type::ALL);
  int nf_owned = num_entities(FACE,Parallel_type::OWNED);
  int nc = num_entities(CELL,Parallel_type::ALL);
  int nc_owned = num_entities(CELL,Parallel_type::OWNED);

  columnID_.resize(nc);
  cell_cellbelow_.resize(nc);
  cell_cellbelow_.assign(nc,-1);
  cell_cellabove_.resize(nc);
  cell_cellabove_.assign(nc,-1);
  node_nodeabove_.resize(nn);
  node_nodeabove_.assign(nn,-1);

  Entity_ID_List top_faces;
  get_set_entities(setname, FACE, Parallel_type::ALL, &top_faces);

  int ncolumns = top_faces.size();
  num_owned_cols_ = get_set_size(setname, FACE, Parallel_type::OWNED);

  int success = 1;
  for (int i = 0; i < ncolumns; i++) {
    Entity_ID f = top_faces[i];
    Entity_ID_List fcells;
    face_get_cells(f,Parallel_type::ALL,&fcells);

    // not a boundary face?
    if (fcells.size() != 1) {
      std::cerr << "Mesh: Provided set for build_columns() includes faces that are not exterior faces.\n";
      success = 0;
      break;
    }

    // check that the normal points upward
    AmanziGeometry::Point normal = face_normal(f,false,fcells[0]);
    normal /= norm(normal);
    if (normal[2] < 1.e-10) {
      std::cerr << "Mesh: Provided set for build_columns() includes faces that don't point upward.\n";
      success = 0;
      break;
    }

    success &= build_single_column_(i, f);
  }

  int min_success;
  get_comm()->MinAll(&success, &min_success, 1);
  columns_built_ = (min_success == 1);
  return columns_built_ ? 1 : 0;
}


// Figure out columns of cells in a semi-structured mesh and cache the
// information for later.
//
// The columns are defined by identifying all boundary faces which
// have a negative-z-direction normal, then collecting cells and faces
// while travelling downward through the columns.  As a
// semi-structured mesh requires that all lateral faces have 0
// z-normal, and all top surface faces have positive z-normal, the set
// of bottom faces is defined with no user input.
//
// NOTE: Currently ghost columns are built too because we didn't know that
// they weren't necessary. --etc
int
Mesh::build_columns() const
{
  if (columns_built_) return 1;

  // Allocate space and initialize.
  int nn = num_entities(NODE,Parallel_type::ALL);
  int nf = num_entities(FACE,Parallel_type::ALL);
  int nf_owned = num_entities(FACE,Parallel_type::OWNED);
  int nc = num_entities(CELL,Parallel_type::ALL);
  int nc_owned = num_entities(CELL, Parallel_type::OWNED);

  columnID_.resize(nc);
  cell_cellbelow_.resize(nc);
  cell_cellbelow_.assign(nc,-1);
  cell_cellabove_.resize(nc);
  cell_cellabove_.assign(nc,-1);
  node_nodeabove_.resize(nn);
  node_nodeabove_.assign(nn,-1);

  // Find the faces at the top of the domain. We assume that these are all
  // the boundary faces whose normal points in the positive z-direction
  //
  int success = 1;
  int ncolumns = 0;
  num_owned_cols_ = 0;
  for (int i = 0; i < nf; i++) {
    Entity_ID_List fcells;
    face_get_cells(i,Parallel_type::ALL,&fcells);

    // Is it a boundary face?
    if (fcells.size() != 1) continue;

    // Is it pointing up?
    AmanziGeometry::Point normal = face_normal(i,false,fcells[0]);
    normal /= norm(normal);

    AmanziGeometry::Point zvec(space_dim_);
    if (space_dim_ == 2)
      zvec.set(0.0,1.0);
    else if (space_dim_ == 3)
      zvec.set(0.0,0.0,1.0);

    double dp = zvec*normal;

    // Check the normal:
    //  1) n dot z = 0 --> lateral face
    //  2) n dot z < 0 --> downard pointing face
    if (dp < 1.e-10) continue;

    success &= build_single_column_(ncolumns, i);
    ncolumns++;
    if (i < nf_owned) num_owned_cols_++;
  }

  int min_success;
  get_comm()->MinAll(&success, &min_success, 1);
  columns_built_ = (min_success == 1);
  return columns_built_ ? 1 : 0;
}


int
Mesh::build_single_column_(int colnum, Entity_ID top_face) const
{
  Entity_ID_List fcells;
  face_get_cells(top_face,Parallel_type::ALL,&fcells);

  // Walk through the cells until we get to the bottom of the domain
  Entity_ID cur_cell = fcells[0];
  bool is_ghost_column = (entity_get_ptype(CELL, cur_cell) == Parallel_type::GHOST);
  Entity_ID bot_face = -1;
  Entity_ID_List fcells2, cfaces, colcells, colfaces;
  std::vector<int> cfdirs;

  AmanziGeometry::Point negzvec(space_dim_);
  if (space_dim_ == 2)
    negzvec.set(0.0,-1.0);
  else if (space_dim_ == 3)
    negzvec.set(0.0,0.0,-1.0);

  int success = 1;
  bool done = false;
  while (!done) {
    bool is_ghost_cell = (entity_get_ptype(CELL, cur_cell) == Parallel_type::GHOST);
    if (is_ghost_column != is_ghost_cell) {
      //      Errors::Message mesg("A column contains cells from different mesh partitions!");
      //      Exceptions::amanzi_throw(mesg);
      std::cerr << "A column contains cells from different mesh partitions" << std::endl;
      success = 0;
      break;

    }
    columnID_[cur_cell] = colnum;
    colcells.push_back(cur_cell);
    colfaces.push_back(top_face);

    // Faces of current cell
    cell_get_faces_and_dirs(cur_cell,&cfaces,&cfdirs);

    // Find the bottom face of the cell as the face whose outward
    // normal from the current cell is most aligned with the -Z
    // direction
    double mindp = -999.0;
    bot_face = -1;
    for (int j = 0; j < cfaces.size(); j++) {
      AmanziGeometry::Point normal = face_normal(cfaces[j]);
      if (cfdirs[j] == -1) normal *= -1;
      normal /= norm(normal);

      double dp = normal*negzvec;
      if (dp > mindp) {
        mindp = dp;
        bot_face = cfaces[j];
      }
    }

    if (bot_face == top_face) {
      std::cerr << "Build Columns broke:" << std::endl
                << "  on column = " << colnum << std::endl
                << "  cell / face = " << cur_cell << "," << bot_face << std::endl
                << "  candidates = " << cfaces[cfaces.size()-2] << "," << cfaces[cfaces.size()-1] << std::endl;
      success = 0;
      break;
    }
    AMANZI_ASSERT(bot_face != top_face);
    AMANZI_ASSERT(bot_face != -1);

    // record the cell above and cell below
    face_get_cells(bot_face,Parallel_type::ALL,&fcells2);
    if (fcells2.size() == 2) {
      if (cell_cellbelow_[cur_cell] != -1) {  // intersecting column of cells
        std::cerr << "Intersecting column of cells\n";
        success = 0;
        break;
      }

      if (fcells2[0] == cur_cell) {
        cell_cellbelow_[cur_cell] = fcells2[1];
        cell_cellabove_[fcells2[1]] = cur_cell;
        cur_cell = fcells2[1];
      }
      else if (fcells2[1] == cur_cell) {
        cell_cellbelow_[cur_cell] = fcells2[0];
        cell_cellabove_[fcells2[0]] = cur_cell;
        cur_cell = fcells2[0];
      }
      else {
        std::cerr << "Unlikely problem in face to cell connectivity\n";
        success = 0;
        break;
      }
    }
    else {
      done = true;
    }

    // record the node below for each of the top face nodes
    // start node of bottom face
    Entity_ID_List botnodes, topnodes, sidenodes;

    face_get_nodes(top_face,&topnodes);
    Entity_ID topnode0 = topnodes[0];

    // nodes of the top face
    face_get_nodes(bot_face,&botnodes);

    if (topnodes.size() != botnodes.size()) {
      std::cerr << "Top and bottom face of columnar cell have different number of nodes.\n";
      success = 0;
      break;
    }

    // match a node below to a node above
    bool found = false;
    int ind = -1;
    int nfvtop = topnodes.size();

    AmanziGeometry::Point topnode0c;
    node_get_coordinates(topnode0, &topnode0c);

    for (int k = 0; k < nfvtop; k++) {
      AmanziGeometry::Point kc;
      node_get_coordinates(botnodes[k], &kc);

      double horiz_dist = 0.;
      for (int m=0; m!=space_dim_-1; ++m) {
        horiz_dist += std::abs(topnode0c[m]-kc[m]);
      }

      if (horiz_dist < 1.e-6) {
        found = true;
        ind = k;
        break;
      }
    }

    if (!found) {
      std::cerr << "Could not find the right structure in mesh\n";
      success = 0;
      break;
    }

    // We have a matching topnode and botnode - now match up the rest
    // even or odd handedness?
    double even_odd_product = face_normal(bot_face)[space_dim_-1]
        * face_normal(top_face)[space_dim_-1];
    AMANZI_ASSERT(std::abs(even_odd_product) > 0);
    int even_odd = even_odd_product >= 0. ? 1 : -1;

    for (int k = 0; k < nfvtop; k++) {
      Entity_ID topnode = topnodes[k];
      int bot_i = (ind+even_odd*k)%nfvtop;
      if (bot_i < 0) bot_i += nfvtop;
      Entity_ID botnode = botnodes[bot_i];
      node_nodeabove_[botnode] = topnode;

      // AMANZI_ASSERT used in debugging
      // AmanziGeometry::Point bc;
      // AmanziGeometry::Point tc;
      // node_get_coordinates(topnode, &bc);
      // node_get_coordinates(botnode, &tc);
      // double horiz_dist = 0.;
      // for (int m=0; m!=space_dim_-1; ++m) {
      //   horiz_dist += std::abs(bc[m]-tc[m]);
      // }
      // AMANZI_ASSERT(horz_dist < 1.e-10);
    }

    top_face = bot_face;
  } // while (!done)

  if (success) {
    colfaces.push_back(bot_face);
    column_cells_.push_back(colcells);
    column_faces_.push_back(colfaces);
  }

  return success;
}


std::string
Mesh::cell_type_to_name(const Cell_type type)
{
  switch (type)
  {
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


void Mesh::PrintMeshStatistics() const
{
  if (vo_.get() && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    int ncells = num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    int nfaces = num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

    int min_out[2], max_out[2], sum_out[2], tmp_in[2] = {ncells, nfaces};
    get_comm()->MinAll(tmp_in, min_out, 2);
    get_comm()->MaxAll(tmp_in, max_out, 2);
    get_comm()->SumAll(tmp_in, sum_out, 2);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "cells, tot/min/max: " << sum_out[0] << "/" << min_out[0] << "/" << max_out[0] << "\n";
    *vo_->os() << "faces, tot/min/max: " << sum_out[1] << "/" << min_out[1] << "/" << max_out[1] << "\n\n";
  }
}

}  // namespace AmanziMesh
}  // namespace Amanzi
