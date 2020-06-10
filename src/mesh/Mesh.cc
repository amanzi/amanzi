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
    cell_get_faces_and_bisectors_precomputed_(false), 
    parents_precomputed_(false),
    parent_(Teuchos::null),
    logical_(false),
    kdtree_faces_initialized_(false)
{
  if (plist_ == Teuchos::null) {
    plist_ = Teuchos::rcp(new Teuchos::ParameterList("Mesh"));
  }
  vo_ = Teuchos::rcp(new VerboseObject(
      Keys::cleanPListName(plist_->name()), *plist_, comm_));
};


void
Mesh::cache_cell_face_info_() const
{
  int ncells = num_entities(CELL, Parallel_type::ALL);
  Kokkos::resize(cell_face_ids_.row_map, ncells + 1);
  Kokkos::resize(cell_face_dirs_.row_map, ncells + 1);
  cell_face_ids_.row_map(0) = 0;
  cell_face_dirs_.row_map(0) = 0;

  int nfaces = num_entities(FACE, Parallel_type::ALL);
  Kokkos::resize(face_cell_ids_.row_map, nfaces + 1);
  Kokkos::resize(face_cell_ptype_.row_map, nfaces + 1);
  face_cell_ids_.row_map(0) = 0;
  face_cell_ptype_.row_map(0) = 0;

  int cell_face_ids_size = 0; 
  int cell_face_dirs_size = 0; 
  for (int c = 0; c < ncells; c++) {
    Entity_ID_List cell_face_ids_view;
    std::vector<int> cell_face_dirs_view;
    cell_get_faces_and_dirs_internal_(
      c, cell_face_ids_view, cell_face_dirs_view);
    int nf = cell_face_ids_view.size();

    for (int jf = 0; jf < nf; jf++) {
      Entity_ID f = cell_face_ids_view[jf];
      int dir = cell_face_dirs_view[jf];
      face_cell_ids_.row_map(f)++;
      face_cell_ptype_.row_map(f)++;
    }

    cell_face_ids_size += cell_face_ids_view.size();
    cell_face_dirs_size += cell_face_dirs_view.size();
  }

  Kokkos::resize(cell_face_ids_.entries,cell_face_ids_size);
  Kokkos::resize(cell_face_dirs_.entries,cell_face_dirs_size); 

  for (int c = 0; c < ncells; c++) {
    Entity_ID_List cell_face_ids_view;
    std::vector<int> cell_face_dirs_view;
    cell_get_faces_and_dirs_internal_(
      c, cell_face_ids_view, cell_face_dirs_view);
    int nf = cell_face_ids_view.size();

    // Save to the CRS
    cell_face_ids_.row_map(c + 1) =
      cell_face_ids_.row_map(c) + cell_face_ids_view.size();
    for (int fi = 0; fi < cell_face_ids_view.size(); ++fi) {
      cell_face_ids_.entries(cell_face_ids_.row_map(c) + fi) =
        cell_face_ids_view[fi];
    }
    cell_face_dirs_.row_map(c + 1) =
      cell_face_dirs_.row_map(c) + cell_face_dirs_view.size();
    for (int fd = 0; fd < cell_face_dirs_view.size(); ++fd) {
      cell_face_dirs_.entries(cell_face_dirs_.row_map(c) + fd) =
        cell_face_dirs_view[fd];
    }
  }

  // Prefix sum
  Entity_ID tmp1 = face_cell_ids_.row_map(0);
  face_cell_ids_.row_map(0) = 0;
  for (int f = 0; f < nfaces; ++f) {
    Entity_ID tmp2 = face_cell_ids_.row_map(f + 1);
    face_cell_ids_.row_map(f + 1) = face_cell_ids_.row_map(f) + tmp1;
    tmp1 = tmp2;
  }
  Kokkos::deep_copy(face_cell_ptype_.row_map, face_cell_ids_.row_map);
  Kokkos::resize(face_cell_ids_.entries,
                 face_cell_ids_.row_map(face_cell_ids_.row_map.extent(0) - 1));
  Kokkos::resize(
    face_cell_ptype_.entries,
    face_cell_ptype_.row_map(face_cell_ptype_.row_map.extent(0) - 1));


  Kokkos::View<int*> offset("", nfaces);

  for (int c = 0; c < ncells; c++) {
    int nf = cell_face_ids_.row_map(c + 1) - cell_face_ids_.row_map(c);
    for (int jf = 0; jf < nf; jf++) {
      Entity_ID f = cell_face_ids_.entries(jf + cell_face_ids_.row_map(c));
      int dir = cell_face_dirs_.entries(jf + cell_face_dirs_.row_map(c));
      face_cell_ids_.entries(face_cell_ids_.row_map(f) + offset(f)) =
        dir > 0 ? c : ~c;
      face_cell_ptype_.entries(face_cell_ptype_.row_map(f) + offset(f)) =
        entity_get_ptype(CELL, c);
      offset(f)++;
    }
  }

  // Sort the cells based on parallel type
  // Basic select sort for now \TODO improve the sort
  for (int f = 0; f < nfaces; f++) {
    int begin = face_cell_ids_.row_map(f);
    int end = face_cell_ids_.row_map(f + 1);
    // Sort face_cell_ids_ regarding face_cell_ptype_
    // Search for the parallel types one by one
    int pos = begin;
    for (int pt = 0; pt < static_cast<int>(Parallel_type::PARALLEL_TYPE_SIZE);
         ++pt) {
      for (int c = pos; c < end; ++c) {
        if (static_cast<int>(face_cell_ptype_.entries(c)) == pt) {
          std::swap(face_cell_ptype_.entries(c), face_cell_ptype_.entries(pos));
          std::swap(face_cell_ids_.entries(c), face_cell_ids_.entries(pos));
          ++pos;
        }
      }
    }
    assert(pos <= end);
  }

  // Copy to the dirs
  Kokkos::resize(face_cell_ids_dirs_, face_cell_ids_.entries.extent(0));
  for (int i = 0; i < face_cell_ids_.entries.extent(0); ++i) {
    if (face_cell_ids_.entries(i) < 0)
      face_cell_ids_dirs_(i) = ~face_cell_ids_.entries(i);
    else
      face_cell_ids_dirs_(i) = face_cell_ids_.entries(i);
  }

  cell2face_info_cached_ = true;
  face2cell_info_cached_ = true;
  faces_requested_ = true;
}

void
Mesh::display_cache()
{
  std::cout << std::endl;
  if (!cell2face_info_cached_) cache_cell_face_info_();
  if (!cell_geometry_precomputed_) compute_cell_geometric_quantities_();
  if (!face_geometry_precomputed_) compute_face_geometric_quantities_();
  // Display all cache
  std::cout << "face_centroid: " << std::endl;
  for (int i = 0; i < face_centroids_.extent(0); ++i)
    std::cout << face_centroids_(i)[0] << " " << face_centroids_(i)[1] << " "
              << face_centroids_(i)[2] << std::endl;
  std::cout << std::endl;

  std::cout << "face_normal: " << std::endl;
  for (int i = 0; i < face_normals_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = face_normals_.row_map(i); j < face_normals_.row_map(i + 1);
         ++j) {
      std::cout << face_normals_.entries(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "face_cell_ptype_: " << std::endl;
  for (int i = 0; i < face_cell_ptype_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = face_cell_ptype_.row_map(i);
         j < face_cell_ptype_.row_map(i + 1);
         ++j) {
      std::cout << static_cast<int>(face_cell_ptype_.entries(j)) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;


  std::cout << "face_cell_ids_: " << std::endl;
  for (int i = 0; i < face_cell_ids_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = face_cell_ids_.row_map(i); j < face_cell_ids_.row_map(i + 1);
         ++j) {
      std::cout << face_cell_ids_.entries(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "cell_face_ids_: " << std::endl;
  for (int i = 0; i < cell_face_ids_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = cell_face_ids_.row_map(i); j < cell_face_ids_.row_map(i + 1);
         ++j) {
      std::cout << cell_face_ids_.entries(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "cell_face_dirs_: " << std::endl;
  for (int i = 0; i < cell_face_dirs_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = cell_face_dirs_.row_map(i); j < cell_face_dirs_.row_map(i + 1);
         ++j) {
      std::cout << cell_face_dirs_.entries(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

// Gather and cache face to edge connectivity info.
void
Mesh::cache_face2edge_info_() const
{
  int nfaces = num_entities(FACE, Parallel_type::ALL);
  Kokkos::resize(face_edge_ids_.row_map, nfaces + 1);
  Kokkos::resize(face_edge_dirs_.row_map, nfaces + 1);
  face_edge_ids_.row_map(0) = 0;
  face_edge_dirs_.row_map(0) = 0;

  int face_edge_ids_size = 0; 
  for (int f = 0; f < nfaces; f++) {
    Entity_ID_List fedgeids;
    std::vector<int> fedgedirs;
    face_get_edges_and_dirs_internal_(f, fedgeids, &fedgedirs, true);
    face_edge_ids_.row_map(f + 1) =
      face_edge_ids_.row_map(f) + fedgeids.size();
    face_edge_dirs_.row_map(f + 1) =
      face_edge_dirs_.row_map(f) + fedgedirs.size();
    face_edge_ids_size += fedgeids.size();
  }

  Kokkos::resize(face_edge_ids_.entries,face_edge_ids_size); 
  Kokkos::resize(face_edge_dirs_.entries,face_edge_ids_size); 

  for (int f = 0; f < nfaces; f++) {
    Entity_ID_List fedgeids;
    std::vector<int> fedgedirs;
    face_get_edges_and_dirs_internal_(f, fedgeids, &fedgedirs, true);
    face_edge_ids_.row_map(f + 1) =
      face_edge_ids_.row_map(f) + fedgeids.size();
    face_edge_dirs_.row_map(f + 1) =
      face_edge_dirs_.row_map(f) + fedgedirs.size();
    for (int i = 0; i < fedgeids.size(); ++i) {
      face_edge_ids_.entries(face_edge_ids_.row_map(f) + i) = fedgeids[i];
      face_edge_dirs_.entries(face_edge_dirs_.row_map(f) + i) = fedgedirs[i];
    }
  }

  face2edge_info_cached_ = true;
  faces_requested_ = true;
  edges_requested_ = true;
}


// Gather and cache cell to edge connectivity info.
void
Mesh::cache_cell2edge_info_() const
{
  int ncells = num_entities(CELL, Parallel_type::ALL);
  Kokkos::resize(cell_edge_ids_.row_map, ncells + 1);

  if (space_dim_ == 2) {
    Kokkos::resize(cell_2D_edge_dirs_.row_map, ncells + 1);
    int cell_edge_ids_entries_size = 0; 
    int cell_2D_edge_dirs_entries_size = 0; 
    for (int c = 0; c < ncells; c++) {
      Entity_ID_List cell_edge_ids_tmp;
      std::vector<int> cell_2D_edge_dirs_tmp;
      cell_2D_get_edges_and_dirs_internal_(
        c, cell_edge_ids_tmp, &cell_2D_edge_dirs_tmp);
      cell_edge_ids_.row_map(c + 1) =
        cell_edge_ids_.row_map(c) + cell_edge_ids_tmp.size();
      cell_2D_edge_dirs_.row_map(c + 1) =
        cell_2D_edge_dirs_.row_map(c) + cell_2D_edge_dirs_tmp.size();
      cell_edge_ids_entries_size += cell_edge_ids_tmp.size(); 
      cell_2D_edge_dirs_entries_size += cell_2D_edge_dirs_tmp.size();
    }
    Kokkos::resize(cell_edge_ids_.entries,cell_edge_ids_entries_size); 
    Kokkos::resize(cell_2D_edge_dirs_.entries,cell_2D_edge_dirs_entries_size);
    for (int c = 0; c < ncells; c++) {
      Entity_ID_List cell_edge_ids_tmp;
      std::vector<int> cell_2D_edge_dirs_tmp;
      cell_2D_get_edges_and_dirs_internal_(
        c, cell_edge_ids_tmp, &cell_2D_edge_dirs_tmp);
      cell_edge_ids_.row_map(c + 1) =
        cell_edge_ids_.row_map(c) + cell_edge_ids_tmp.size();
      cell_2D_edge_dirs_.row_map(c + 1) =
        cell_2D_edge_dirs_.row_map(c) + cell_2D_edge_dirs_tmp.size();
      for (int i = 0; i < cell_edge_ids_tmp.size(); ++i) {
        cell_edge_ids_.entries(cell_edge_ids_.row_map(c) + i) =
          cell_edge_ids_tmp[i];
        cell_2D_edge_dirs_.entries(cell_2D_edge_dirs_.row_map(c) + i) =
          cell_2D_edge_dirs_tmp[i];
      }
    }
  } else  {
    int cell_edge_ids_size = 0; 
    for (int c = 0; c < ncells; c++) {
      Entity_ID_List cell_edge_ids_tmp;
      cell_get_edges_internal_(c, cell_edge_ids_tmp);
      cell_edge_ids_.row_map(c + 1) =
        cell_edge_ids_.row_map(c) + cell_edge_ids_tmp.size();
      cell_edge_ids_size += cell_edge_ids_tmp.size();
    }
    Kokkos::resize(cell_edge_ids_.entries,cell_edge_ids_size); 
    for (int c = 0; c < ncells; c++) {
      Entity_ID_List cell_edge_ids_tmp;
      cell_get_edges_internal_(c, cell_edge_ids_tmp);
      cell_edge_ids_.row_map(c + 1) =
        cell_edge_ids_.row_map(c) + cell_edge_ids_tmp.size();
      for (int i = 0; i < cell_edge_ids_tmp.size(); ++i) {
        cell_edge_ids_.entries(cell_edge_ids_.row_map(c) + i) =
          cell_edge_ids_tmp[i];
      }
    }
  }
  cell2edge_info_cached_ = true;
}

Entity_ID
Mesh::entity_get_parent_type(const Entity_kind kind, const Entity_ID entid) const
{
  Errors::Message mesg(
    "Parent/daughter entities not enabled in this framework.");
  Exceptions::amanzi_throw(mesg);
  return -1;
}

void 
Mesh::cache_parents_info_() const 
{
  // CELLS
  int ncells = num_entities(CELL, Parallel_type::ALL);
  Kokkos::resize(cells_parent_,ncells); 
  for(int i = 0 ; i < ncells; ++i){
    cells_parent_[i] = entity_get_parent_type(CELL,i); 
  }
  // FACES 
  int nfaces = num_entities(FACE, Parallel_type::ALL);
  Kokkos::resize(faces_parent_,nfaces);
  for(int i = 0 ; i < nfaces; ++i){
    faces_parent_[i] = entity_get_parent_type(FACE,i); 
  }
  // NODES 
  int nnodes = num_entities(NODE, Parallel_type::ALL);
  Kokkos::resize(nodes_parent_,nnodes);
  for(int i = 0 ; i < nnodes; ++i){
    nodes_parent_[i] = entity_get_parent_type(NODE,i); 
  }
  // EDGES
  if(edges_requested_){
    int nedges = num_entities(EDGE, Parallel_type::ALL);
    Kokkos::resize(edges_parent_,nedges);
    for(int i = 0 ; i < nedges; ++i){
      edges_parent_[i] = entity_get_parent_type(EDGE,i); 
    }
  }
  parents_precomputed_ = true; 
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
  if (edges_requested_) {
    int ncells = num_entities(CELL, Parallel_type::OWNED);
    for (int c = 0; c < ncells; ++c) {
      Kokkos::View<AmanziMesh::Entity_ID*> edges;
      cell_get_edges(c, edges);
      n = std::max(n, (unsigned int)edges.extent(0));
    }
  }
  return n;
}

void
Mesh::init_cache()
{
  assert(!cell2face_info_cached_);
  cache_cell_face_info_();
  cell2face_info_cached_ = true;
  assert(!cell_geometry_precomputed_);
  compute_cell_geometric_quantities_();
  cell_geometry_precomputed_ = true;
  assert(!face_geometry_precomputed_);
  compute_face_geometric_quantities_();
  face_geometry_precomputed_ = true;
  if(parent().get()){
    assert(!parents_precomputed_); 
    cache_parents_info_(); 
    parents_precomputed_ = true; 
  }
  assert(!cell_get_faces_and_bisectors_precomputed_); 
  cache_cell_get_faces_and_bisectors_(); 
  cell_get_faces_and_bisectors_precomputed_ = true;
}

void 
Mesh::cache_cell_get_faces_and_bisectors_() const {
  int ncells = num_entities(CELL, Parallel_type::ALL);
  Kokkos::resize(cell_faces_bisectors_.row_map,ncells+1); 
  cell_faces_bisectors_.row_map(0) = 0; 
  
  int entries_size = 0; 
  for(int i = 0 ; i < ncells; ++i){
    Kokkos::View<Entity_ID*> faceids;
    cell_get_faces(i, faceids);
    entries_size += faceids.size(); 
  }
  Kokkos::resize(cell_faces_bisectors_.entries,entries_size); 
  for(int i = 0 ; i < ncells; ++i){
    Kokkos::View<Entity_ID*> faceids;
    cell_get_faces(i, faceids);
    AmanziGeometry::Point cc = cell_centroid(i);
    cell_faces_bisectors_.row_map(i + 1) = faceids.extent(0) +
       cell_faces_bisectors_.row_map(i);
    for (int j = 0; j < faceids.extent(0); ++j) {
      cell_faces_bisectors_.entries(cell_faces_bisectors_.row_map(i)+j) = 
        face_centroid(faceids(j)) - cc;
    }
  }
}

void
Mesh::face_get_edges_and_dirs(const Entity_ID faceid,
                              Kokkos::View<Entity_ID*>& edgeids,
                              Kokkos::View<int*>* edge_dirs,
                              const bool ordered) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!face2edge_info_cached_) cache_face2edge_info_();

  edgeids = Kokkos::subview(face_edge_ids_.entries,
                            std::make_pair(face_edge_ids_.row_map(faceid),
                                           face_edge_ids_.row_map(faceid + 1)));
  //*edgeids = face_edge_ids_[faceid]; // copy operation

  if (edge_dirs) {
    *edge_dirs =
      Kokkos::subview(face_edge_dirs_.entries,
                      std::make_pair(face_edge_dirs_.row_map(faceid),
                                     face_edge_dirs_.row_map(faceid + 1)));
    // std::vector<int> &fedgedirs = face_edge_dirs_[faceid];
    //*edge_dirs = fedgedirs; // copy operation
  }

#else // Non-cached version
  face_get_edges_and_dirs_internal_(faceid, edgeids, edge_dirs, ordered);

#endif
}


// Get the local ID of a face edge in a cell edge list
void
Mesh::face_to_cell_edge_map(const Entity_ID faceid, const Entity_ID cellid,
                            std::vector<int>* map) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!face2edge_info_cached_) cache_face2edge_info_();
  if (!cell2edge_info_cached_) cache_cell2edge_info_();

  size_t faceid_size =
    face_edge_ids_.row_map(faceid + 1) - face_edge_ids_.row_map(faceid);
  map->resize(faceid_size);
  for (int f = 0; f < faceid_size; ++f) {
    Entity_ID fedge =
      face_edge_ids_.entries(face_edge_ids_.row_map(faceid) + f);

    size_t cellid_size =
      cell_edge_ids_.row_map(cellid + 1) - cell_edge_ids_.row_map(cellid);
    for (int c = 0; c < cellid_size; ++c) {
      if (fedge == cell_edge_ids_.entries(cell_edge_ids_.row_map(cellid) + c)) {
        (*map)[f] = c;
        break;
      }
    }
  }

#else // non-cached version

  Kokkos::View<Entity_ID*> fedgeids, cedgeids;
  Kokkos::View<int*> fedgedirs;

  face_get_edges_and_dirs(faceid, fedgeids, &fedgedirs, true);
  cell_get_edges(cellid, cedgeids);

  map->resize(fedgeids.extent(0), -1);
  for (int f = 0; f < fedgeids.extent(0); ++f) {
    Entity_ID fedge = fedgeids(f);

    for (int c = 0; c < cedgeids.extent(0); ++c) {
      if (fedge == cedgeids(c)) {
        (*map)[f] = c;
        break;
      }
    }
  }

#endif
}


void
Mesh::cell_get_edges(const Entity_ID cellid,
                     Kokkos::View<Entity_ID*>& edgeids) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2edge_info_cached_) cache_cell2edge_info_();

  edgeids = Kokkos::subview(cell_edge_ids_.entries,
                            std::make_pair(cell_edge_ids_.row_map(cellid),
                                           cell_edge_ids_.row_map(cellid + 1)));
  // Entity_ID_List &cedgeids = cell_edge_ids_[cellid];
  //*edgeids = cell_edge_ids_[cellid]; // copy operation

#else // Non-cached version
  cell_get_edges_internal_(cellid, edgeids);

#endif
}


void
Mesh::cell_2D_get_edges_and_dirs(const Entity_ID cellid,
                                 Kokkos::View<Entity_ID*>& edgeids,
                                 Kokkos::View<int*>* edgedirs) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2edge_info_cached_) cache_cell2edge_info_();

  edgeids = Kokkos::subview(cell_edge_ids_.entries,
                            std::make_pair(cell_edge_ids_.row_map(cellid),
                                           cell_edge_ids_.row_map(cellid + 1)));
  *edgedirs =
    Kokkos::subview(cell_2D_edge_dirs_.entries,
                    std::make_pair(cell_2D_edge_dirs_.row_map(cellid),
                                   cell_2D_edge_dirs_.row_map(cellid + 1)));

#else // Non-cached version
  cell_2D_get_edges_and_dirs_internal_(cellid, edgeids, edgedirs);

#endif
}


int
Mesh::compute_cell_geometric_quantities_() const
{
  int ncells = num_entities(CELL, Parallel_type::ALL);

  Kokkos::resize(cell_volumes_, ncells);
  Kokkos::resize(cell_centroids_, ncells);
  for (int i = 0; i < ncells; i++) {
    double volume;
    AmanziGeometry::Point centroid(space_dim_);

    compute_cell_geometry_(i, &volume, &centroid);

    cell_volumes_(i) = volume;
    cell_centroids_(i) = centroid;
  }

  cell_geometry_precomputed_ = true;

  return 1;
}


int
Mesh::compute_face_geometric_quantities_() const
{
  if (space_dimension() == 3 && manifold_dimension() == 2) {
    // need cell centroids to compute normals
    if (!cell_geometry_precomputed_) compute_cell_geometric_quantities_();
  }

  int nfaces = num_entities(FACE, Parallel_type::ALL);

  Kokkos::resize(face_areas_, nfaces);
  Kokkos::resize(face_centroids_, nfaces);

  // Temporary views
  Kokkos::resize(face_normals_.row_map, nfaces + 1);
  face_normals_.row_map(0) = 0;

  // Find size 
  int entries_size = 0; 
  for (int i = 0; i < nfaces; i++) {
    double area;
    AmanziGeometry::Point centroid(space_dim_);
    std::vector<AmanziGeometry::Point> normals;
    // normal0 and normal1 are outward normals of the face with
    // respect to the cell0 and cell1 of the face. The natural normal
    // of the face points out of cell0 and into cell1. If one of these
    // cells do not exist, then the normal is the null vector.
    compute_face_geometry_(i, &area, &centroid, normals);
    face_areas_(i) = area;
    face_centroids_(i) = centroid;
    entries_size += normals.size(); 
  }

  Kokkos::resize(face_normals_.entries,entries_size); 

  for (int i = 0; i < nfaces; i++) {
    double area;
    AmanziGeometry::Point centroid(space_dim_);
    std::vector<AmanziGeometry::Point> normals;
    compute_face_geometry_(i, &area, &centroid, normals);
    for (int j = 0; j < normals.size(); ++j) {
      face_normals_.entries(face_normals_.row_map(i) + j) = normals[j];
    }
    face_normals_.row_map(i + 1) = normals.size() + face_normals_.row_map(i);
  }

  face_geometry_precomputed_ = true;

  return 1;
}


int
Mesh::compute_edge_geometric_quantities_() const
{
  int nedges = num_entities(EDGE, Parallel_type::ALL);

  Kokkos::resize(edge_vectors_, nedges);
  Kokkos::resize(edge_lengths_, nedges);

  for (int i = 0; i < nedges; i++) {
    double length;
    AmanziGeometry::Point evector(space_dim_);

    compute_edge_geometry_(i, &length, &evector);

    edge_lengths_(i) = length;
    edge_vectors_(i) = evector;
  }

  edge_geometry_precomputed_ = true;

  return 1;
}


int
Mesh::compute_cell_geometry_(const Entity_ID cellid, double* volume,
                             AmanziGeometry::Point* centroid) const
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
    Kokkos::View<Entity_ID*> faces;
    std::vector<unsigned int> nfnodes;
    Kokkos::View<int*> fdirs;
    std::vector<AmanziGeometry::Point> fcoords, ccoords, cfcoords; 

    cell_get_faces_and_dirs(cellid, faces, fdirs);

    int nf = faces.extent(0);
    nfnodes.resize(nf);

    if (nf == 2) {
      /* special case of column mesh - only top and bottom faces
         are returned */
      AmanziGeometry::Point fcentroid0(space_dim_), fcentroid1(space_dim_);
      AmanziGeometry::Point normal(space_dim_);
      double farea;

      /* compute volume on the assumption that the top and bottom faces form
         a vertical columnar cell or in other words a polygonal prism */

      face_get_coordinates(faces(0), fcoords);
      AmanziGeometry::polygon_get_area_centroid_normal(
        fcoords, &farea, &fcentroid0, &normal);

      face_get_coordinates(faces(1), fcoords);
      AmanziGeometry::polygon_get_area_centroid_normal(
        fcoords, &farea, &fcentroid1, &normal);

      *centroid = (fcentroid0 + fcentroid1) / 2.0;
      double height = norm(fcentroid1 - fcentroid0);

      *volume = farea * height;
    } else { /* general case */

      size_t cfcoords_size = 0;
      for (int j = 0; j < nf; j++) {
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

      AmanziGeometry::polyhed_get_vol_centroid(
        ccoords, nf, nfnodes, cfcoords, volume, centroid);
    }
    return 1;
  } else if (manifold_dim_ == 2) {
    std::vector<AmanziGeometry::Point> ccoords;
    cell_get_coordinates(cellid, ccoords);

    AmanziGeometry::Point normal(space_dim_);

    AmanziGeometry::polygon_get_area_centroid_normal(
      ccoords, volume, centroid, &normal);
    return 1;
  }

  return 0;
}


int
Mesh::compute_face_geometry_(
  const Entity_ID faceid, double* area, AmanziGeometry::Point* centroid,
  std::vector<AmanziGeometry::Point>& normals) const
{
  std::vector<AmanziGeometry::Point> fcoords;
  // normals->clear();

  if (manifold_dim_ == 3) {
    // 3D Elements with possibly curved faces

    face_get_coordinates(faceid, fcoords);

    AmanziGeometry::Point normal(3);
    AmanziGeometry::polygon_get_area_centroid_normal(
      fcoords, area, centroid, &normal);

    Kokkos::View<Entity_ID*> cellids;
    face_get_cells(faceid, Parallel_type::ALL, cellids);
    AMANZI_ASSERT(cellids.extent(0) <= 2);

    normals.resize(cellids.size()); 
    for (int i = 0; i < cellids.extent(0); i++) {
      Kokkos::View<Entity_ID*> cellfaceids;
      Kokkos::View<int*> cellfacedirs;
      int dir = 1;

      cell_get_faces_and_dirs(cellids(i), cellfaceids, cellfacedirs);

      bool found = false;
      for (int j = 0; j < cellfaceids.extent(0); j++) {
        if (cellfaceids(j) == faceid) {
          found = true;
          dir = cellfacedirs(j);
          break;
        }
      }

      AMANZI_ASSERT(found);

      normals[i] = (dir == 1) ? normal : -normal;
    }

    return 1;
  } else if (manifold_dim_ == 2) {
    if (space_dim_ == 2) { // 2D mesh

      face_get_coordinates(faceid, fcoords);

      AmanziGeometry::Point evec = fcoords[1] - fcoords[0];
      *area = sqrt(evec * evec);

      *centroid = 0.5 * (fcoords[0] + fcoords[1]);

      AmanziGeometry::Point normal(evec[1], -evec[0]);

      Kokkos::View<Entity_ID*> cellids;
      face_get_cells(faceid, Parallel_type::ALL, cellids);
      AMANZI_ASSERT(cellids.extent(0) <= 2);

      normals.resize(cellids.size()); 
      // normals->resize(cellids.size(), AmanziGeometry::Point(0.0, 0.0));
      for (int i = 0; i < cellids.extent(0); i++) {
        Kokkos::View<Entity_ID*> cellfaceids;
        Kokkos::View<int*> cellfacedirs;
        int dir = 1;

        cell_get_faces_and_dirs(cellids(i), cellfaceids, cellfacedirs);

        bool found = false;
        for (int j = 0; j < cellfaceids.extent(0); j++) {
          if (cellfaceids(j) == faceid) {
            found = true;
            dir = cellfacedirs(j);
            break;
          }
        }

        AMANZI_ASSERT(found);

        normals[i] = (dir == 1) ? normal : -normal;
      }

      return 1;
    } else { // Surface mesh - cells are 2D, coordinates are 3D

      // Since the edge likely forms a discontinuity in the surface
      // (or may even be the intersection of several surfaces), we
      // have to compute an outward normal to the edge with respect to
      // each face

      face_get_coordinates(faceid, fcoords);

      AmanziGeometry::Point evec = fcoords[1] - fcoords[0];
      *area = sqrt(evec * evec);

      *centroid = 0.5 * (fcoords[0] + fcoords[1]);

      Kokkos::View<Entity_ID*> cellids;
      face_get_cells(faceid, Parallel_type::ALL, cellids);

      normals.resize(cellids.extent(0),  AmanziGeometry::Point(0.0, 0.0, 0.0)); 
      for (int i = 0; i < cellids.extent(0); i++) {
        Kokkos::View<Entity_ID*> cellfaceids;
        Kokkos::View<int*> cellfacedirs;
        // int dir = 1; un-used

        cell_get_faces_and_dirs(cellids(i), cellfaceids, cellfacedirs);

        bool found = false;
        for (int j = 0; j < cellfaceids.extent(0); j++) {
          if (cellfaceids(j) == faceid) {
            found = true;
            cellfacedirs(j);
            break;
          }
        }

        AMANZI_ASSERT(found);

        AmanziGeometry::Point cvec = fcoords[0] - cell_centroids_(cellids(i));
        AmanziGeometry::Point trinormal = cvec ^ evec;

        AmanziGeometry::Point normal = evec ^ trinormal;

        double len = norm(normal);
        normal /= len;
        normal *= *area;

        normals[i] = normal; // Always an outward normal as calculated
      }

      return 1;
    }
  }
  return 0;
}


int
Mesh::compute_edge_geometry_(const Entity_ID edgeid, double* edge_length,
                             AmanziGeometry::Point* edge_vector) const
{
  (*edge_vector).set(0.0L);
  *edge_length = 0.0;

  Entity_ID node0, node1;

  edge_get_nodes(edgeid, &node0, &node1);

  AmanziGeometry::Point point0, point1;
  node_get_coordinates(node0, &point0);
  node_get_coordinates(node1, &point1);

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
    return cell_volumes_(cellid);
  } else {
    if (recompute) {
      double volume;
      AmanziGeometry::Point centroid(space_dim_);
      compute_cell_geometry_(cellid, &volume, &centroid);
      return volume;
    } else
      return cell_volumes_(cellid);
  }
}

// Length of an edge
double
Mesh::edge_length(const Entity_ID edgeid, const bool recompute) const
{
  AMANZI_ASSERT(edges_requested_);

  if (!edge_geometry_precomputed_) {
    compute_edge_geometric_quantities_();
    return edge_lengths_(edgeid);
  } else {
    if (recompute) {
      double length;
      AmanziGeometry::Point vector(space_dim_);
      compute_edge_geometry_(edgeid, &length, &vector);
      return length;
    } else
      return edge_lengths_(edgeid);
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
  return (xyz0 + xyz1) / 2;
}


// Direction vector of edge
AmanziGeometry::Point
Mesh::edge_vector(const Entity_ID edgeid, const bool recompute,
                  const Entity_ID pointid, int* orientation) const
{
  AMANZI_ASSERT(edges_requested_);

  AmanziGeometry::Point evector(space_dim_);
  AmanziGeometry::Point& evector_ref = evector; // to avoid extra copying

  if (!edge_geometry_precomputed_) compute_edge_geometric_quantities_();

  if (recompute) {
    double length;
    compute_edge_geometry_(edgeid, &length, &evector);
    // evector_ref already points to evector
  } else
    evector_ref = edge_vectors_(edgeid);

  if (orientation) *orientation = 1;

  if (pointid == -1)
    return evector_ref;
  else {
    Entity_ID p0, p1;
    edge_get_nodes(edgeid, &p0, &p1);

    if (pointid == p0)
      return evector_ref;
    else {
      if (orientation) *orientation = -1;
      return -evector_ref;
    }
  }
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
  compute_cell_geometric_quantities_();
  if (faces_requested_) compute_face_geometric_quantities_();
  if (edges_requested_) compute_edge_geometric_quantities_();

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
  int nn = num_entities(NODE, Parallel_type::ALL);
  int nf = num_entities(FACE, Parallel_type::ALL);
  int nf_owned = num_entities(FACE, Parallel_type::OWNED);
  int nc = num_entities(CELL, Parallel_type::ALL);
  int nc_owned = num_entities(CELL, Parallel_type::OWNED);

  Kokkos::resize(columnID_, nc);
  Kokkos::resize(cell_cellbelow_, nc);
  Kokkos::resize(cell_cellabove_, nc);
  Kokkos::resize(node_nodeabove_, nn);

  Kokkos::parallel_for(
    "Mesh::build_columns loop 1",
    nc, KOKKOS_LAMBDA(const int& i) {
    cell_cellbelow_(i) = -1;
    cell_cellabove_(i) = -1;
  });
  Kokkos::parallel_for(
    "Mesh::build_columns loop 2",
    nn, KOKKOS_LAMBDA(const int& i) { node_nodeabove_(i) = -1; });

  Entity_ID_List top_faces;
  get_set_entities(setname, FACE, Parallel_type::ALL, top_faces);

  int ncolumns = top_faces.size();
  num_owned_cols_ = get_set_size(setname, FACE, Parallel_type::OWNED);

  int success = 1;
  for (int i = 0; i < ncolumns; i++) {
    Entity_ID f = top_faces[i];
    Kokkos::View<Entity_ID*> fcells;
    face_get_cells(f, Parallel_type::ALL, fcells);

    // not a boundary face?
    if (fcells.extent(0) != 1) {
      std::cerr << "Mesh: Provided set for build_columns() includes faces that "
                   "are not exterior faces.\n";
      success = 0;
      break;
    }

    // check that the normal points upward
    AmanziGeometry::Point normal = face_normal(f, false, fcells(0));
    normal /= norm(normal);
    if (normal[2] < 1.e-10) {
      std::cerr << "Mesh: Provided set for build_columns() includes faces that "
                   "don't point upward.\n";
      success = 0;
      break;
    }

    success &= build_single_column_(i, f);
  }

  int min_success;
  Teuchos::reduceAll(
    *get_comm(), Teuchos::REDUCE_MIN, 1, &success, &min_success);
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
  int nn = num_entities(NODE, Parallel_type::ALL);
  int nf = num_entities(FACE, Parallel_type::ALL);
  int nf_owned = num_entities(FACE, Parallel_type::OWNED);
  int nc = num_entities(CELL, Parallel_type::ALL);
  int nc_owned = num_entities(CELL, Parallel_type::OWNED);

  Kokkos::resize(columnID_, nc);
  Kokkos::resize(cell_cellbelow_, nc);
  Kokkos::resize(cell_cellabove_, nc);
  Kokkos::resize(node_nodeabove_, nn);

  Kokkos::parallel_for(
    "Mesh::build_columns loop 1",
    nc, KOKKOS_LAMBDA(const int& i) {
    cell_cellbelow_(i) = -1;
    cell_cellabove_(i) = -1;
  });
  Kokkos::parallel_for(
    "Mesh::build_columns loop 2",
    nn, KOKKOS_LAMBDA(const int& i) { node_nodeabove_(i) = -1; });

  // Find the faces at the top of the domain. We assume that these are all
  // the boundary faces whose normal points in the positive z-direction
  //
  int success = 1;
  int ncolumns = 0;
  num_owned_cols_ = 0;
  for (int i = 0; i < nf; i++) {
    Kokkos::View<Entity_ID*> fcells;
    face_get_cells(i, Parallel_type::ALL, fcells);

    // Is it a boundary face?
    if (fcells.extent(0) != 1) continue;

    // Is it pointing up?
    AmanziGeometry::Point normal = face_normal(i, false, fcells(0));
    normal /= norm(normal);

    AmanziGeometry::Point zvec(space_dim_);
    if (space_dim_ == 2)
      zvec.set(0.0, 1.0);
    else if (space_dim_ == 3)
      zvec.set(0.0, 0.0, 1.0);

    double dp = zvec * normal;

    // Check the normal:
    //  1) n dot z = 0 --> lateral face
    //  2) n dot z < 0 --> downard pointing face
    if (dp < 1.e-10) continue;

    success &= build_single_column_(ncolumns, i);
    ncolumns++;
    if (i < nf_owned) num_owned_cols_++;
  }

  int min_success;
  Teuchos::reduceAll(
    *get_comm(), Teuchos::REDUCE_MIN, 1, &success, &min_success);
  columns_built_ = (min_success == 1);
  return columns_built_ ? 1 : 0;
}


int
Mesh::build_single_column_(int colnum, Entity_ID top_face) const
{
  Kokkos::View<Entity_ID*> fcells;
  face_get_cells(top_face, Parallel_type::ALL, fcells);

  // Walk through the cells until we get to the bottom of the domain
  Entity_ID cur_cell = fcells(0);
  bool is_ghost_column =
    (entity_get_ptype(CELL, cur_cell) == Parallel_type::GHOST);
  Entity_ID bot_face = -1;
  Entity_ID_List colcells, colfaces;
  Kokkos::View<Entity_ID*> cfaces, fcells2;
  Kokkos::View<int*> cfdirs;

  AmanziGeometry::Point negzvec(space_dim_);
  if (space_dim_ == 2)
    negzvec.set(0.0, -1.0);
  else if (space_dim_ == 3)
    negzvec.set(0.0, 0.0, -1.0);

  int success = 1;
  bool done = false;
  while (!done) {
    bool is_ghost_cell =
      (entity_get_ptype(CELL, cur_cell) == Parallel_type::GHOST);
    if (is_ghost_column != is_ghost_cell) {
      //      Errors::Message mesg("A column contains cells from different mesh
      //      partitions!"); Exceptions::amanzi_throw(mesg);
      std::cerr << "A column contains cells from different mesh partitions"
                << std::endl;
      success = 0;
      break;
    }
    columnID_(cur_cell) = colnum;
    colcells.push_back(cur_cell);
    colfaces.push_back(top_face);

    // Faces of current cell
    cell_get_faces_and_dirs(cur_cell, cfaces, cfdirs);

    // Find the bottom face of the cell as the face whose outward
    // normal from the current cell is most aligned with the -Z
    // direction
    double mindp = -999.0;
    bot_face = -1;
    for (int j = 0; j < cfaces.extent(0); j++) {
      AmanziGeometry::Point normal = face_normal(cfaces(j));
      if (cfdirs(j) == -1) normal *= -1;
      normal /= norm(normal);

      double dp = normal * negzvec;
      if (dp > mindp) {
        mindp = dp;
        bot_face = cfaces(j);
      }
    }

    if (bot_face == top_face) {
      std::cerr << "Build Columns broke:" << std::endl
                << "  on column = " << colnum << std::endl
                << "  cell / face = " << cur_cell << "," << bot_face
                << std::endl
                << "  candidates = " << cfaces(cfaces.extent(0) - 2) << ","
                << cfaces(cfaces.extent(0) - 1) << std::endl;
      success = 0;
      break;
    }
    AMANZI_ASSERT(bot_face != top_face);
    AMANZI_ASSERT(bot_face != -1);

    // record the cell above and cell below
    face_get_cells(bot_face, Parallel_type::ALL, fcells2);
    if (fcells2.extent(0) == 2) {
      if (cell_cellbelow_(cur_cell) != -1) { // intersecting column of cells
        std::cerr << "Intersecting column of cells\n";
        success = 0;
        break;
      }

      if (fcells2(0) == cur_cell) {
        cell_cellbelow_(cur_cell) = fcells2(1);
        cell_cellabove_(fcells2(1)) = cur_cell;
        cur_cell = fcells2(1);
      } else if (fcells2(1) == cur_cell) {
        cell_cellbelow_(cur_cell) = fcells2(0);
        cell_cellabove_(fcells2(0)) = cur_cell;
        cur_cell = fcells2(0);
      } else {
        std::cerr << "Unlikely problem in face to cell connectivity\n";
        success = 0;
        break;
      }
    } else {
      done = true;
    }

    // record the node below for each of the top face nodes
    // start node of bottom face
    Entity_ID_List sidenodes;
    Entity_ID_List topnodes, botnodes;

    face_get_nodes(top_face, topnodes);
    Entity_ID topnode0 = topnodes[0];

    // nodes of the top face
    face_get_nodes(bot_face, botnodes);

    if (topnodes.size() != botnodes.size()) {
      std::cerr << "Top and bottom face of columnar cell have different number "
                   "of nodes.\n";
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
      for (int m = 0; m != space_dim_ - 1; ++m) {
        horiz_dist += std::abs(topnode0c[m] - kc[m]);
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
    double even_odd_product = face_normal(bot_face)[space_dim_ - 1] *
                              face_normal(top_face)[space_dim_ - 1];
    AMANZI_ASSERT(std::abs(even_odd_product) > 0);
    int even_odd = even_odd_product >= 0. ? 1 : -1;

    for (int k = 0; k < nfvtop; k++) {
      Entity_ID topnode = topnodes[k];
      int bot_i = (ind + even_odd * k) % nfvtop;
      if (bot_i < 0) bot_i += nfvtop;
      Entity_ID botnode = botnodes[bot_i];
      node_nodeabove_(botnode) = topnode;
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

} // namespace AmanziMesh
} // namespace Amanzi
