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

MeshCache::MeshCache(Mesh* mesh, const bool request_faces, const bool request_edges)
  : faces_requested_(request_faces),
    edges_requested_(request_edges),
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
    kdtree_faces_initialized_(false),
    mesh_(mesh)
{};


void
MeshCache::cache_cell_face_info_() const
{
  int ncells = mesh_->num_entities(CELL, Parallel_type::ALL);
  cell_face_ids_.row_map.resize(ncells + 1);
  cell_face_dirs_.row_map.resize(ncells + 1);
  cell_face_ids_.row_map.view_host()(0) = 0;
  cell_face_dirs_.row_map.view_host()(0) = 0;

  int nfaces = mesh_->num_entities(FACE, Parallel_type::ALL);
  face_cell_ids_.row_map.resize(nfaces + 1);
  face_cell_ptype_.row_map.resize(nfaces + 1);
  face_cell_ids_.row_map.view_host()(0) = 0;
  face_cell_ptype_.row_map.view_host()(0) = 0;

  int cell_face_ids_size = 0; 
  int cell_face_dirs_size = 0; 
  for (int c = 0; c < ncells; c++) {
    Entity_ID_List cell_face_ids_view;
    std::vector<int> cell_face_dirs_view;
    mesh_->cell_get_faces_and_dirs_internal_(
      c, cell_face_ids_view, cell_face_dirs_view);
    int nf = cell_face_ids_view.size();

    for (int jf = 0; jf < nf; jf++) {
      Entity_ID f = cell_face_ids_view[jf];
      int dir = cell_face_dirs_view[jf];
      face_cell_ids_.row_map.view_host()(f)++;
      face_cell_ptype_.row_map.view_host()(f)++;
    }

    cell_face_ids_size += cell_face_ids_view.size();
    cell_face_dirs_size += cell_face_dirs_view.size();
  }

  cell_face_ids_.entries.resize(cell_face_ids_size);
  cell_face_dirs_.entries.resize(cell_face_dirs_size); 

  for (int c = 0; c < ncells; c++) {
    Entity_ID_List cell_face_ids_view;
    std::vector<int> cell_face_dirs_view;
    mesh_->cell_get_faces_and_dirs_internal_(
      c, cell_face_ids_view, cell_face_dirs_view);
    int nf = cell_face_ids_view.size();

    // Save to the CRS
    cell_face_ids_.row_map.view_host()(c + 1) =
      cell_face_ids_.row_map.view_host()(c) + cell_face_ids_view.size();
    for (int fi = 0; fi < cell_face_ids_view.size(); ++fi) {
      cell_face_ids_.entries.view_host()(cell_face_ids_.row_map.view_host()(c) + fi) =
        cell_face_ids_view[fi];
    }
    cell_face_dirs_.row_map.view_host()(c + 1) =
      cell_face_dirs_.row_map.view_host()(c) + cell_face_dirs_view.size();
    for (int fd = 0; fd < cell_face_dirs_view.size(); ++fd) {
      cell_face_dirs_.entries.view_host()(cell_face_dirs_.row_map.view_host()(c) + fd) =
        cell_face_dirs_view[fd];
    }
  }

  // Prefix sum
  Entity_ID tmp1 = face_cell_ids_.row_map.view_host()(0);
  face_cell_ids_.row_map.view_host()(0) = 0;
  for (int f = 0; f < nfaces; ++f) {
    Entity_ID tmp2 = face_cell_ids_.row_map.view_host()(f + 1);
    face_cell_ids_.row_map.view_host()(f + 1) = face_cell_ids_.row_map.view_host()(f) + tmp1;
    tmp1 = tmp2;
  }
  Kokkos::deep_copy(face_cell_ptype_.row_map.view_host(), face_cell_ids_.row_map.view_host());
  face_cell_ids_.entries.resize(
                 face_cell_ids_.row_map.view_host()(face_cell_ids_.row_map.extent(0) - 1));
  face_cell_ptype_.entries.resize(
    face_cell_ptype_.row_map.view_host()(face_cell_ptype_.row_map.view_host().extent(0) - 1));


  Kokkos::View<int*, Kokkos::HostSpace> offset("", nfaces);

  for (int c = 0; c < ncells; c++) {
    int nf = cell_face_ids_.row_map.view_host()(c + 1) - cell_face_ids_.row_map.view_host()(c);
    for (int jf = 0; jf < nf; jf++) {
      Entity_ID f = cell_face_ids_.entries.view_host()(jf + cell_face_ids_.row_map.view_host()(c));
      int dir = cell_face_dirs_.entries.view_host()(jf + cell_face_dirs_.row_map.view_host()(c));
      face_cell_ids_.entries.view_host()(face_cell_ids_.row_map.view_host()(f) + offset(f)) =
        dir > 0 ? c : ~c;
      face_cell_ptype_.entries.view_host()(face_cell_ptype_.row_map.view_host()(f) + offset(f)) =
        mesh_->entity_get_ptype(CELL, c);
      offset(f)++;
    }
  }

  // Sort the cells based on parallel type
  // Basic select sort for now \TODO improve the sort
  for (int f = 0; f < nfaces; f++) {
    int begin = face_cell_ids_.row_map.view_host()(f);
    int end = face_cell_ids_.row_map.view_host()(f + 1);
    // Sort face_cell_ids_ regarding face_cell_ptype_
    // Search for the parallel types one by one
    int pos = begin;
    for (int pt = 0; pt < static_cast<int>(Parallel_type::PARALLEL_TYPE_SIZE);
         ++pt) {
      for (int c = pos; c < end; ++c) {
        if (static_cast<int>(face_cell_ptype_.entries.view_host()(c)) == pt) {
          std::swap(face_cell_ptype_.entries.view_host()(c), face_cell_ptype_.entries.view_host()(pos));
          std::swap(face_cell_ids_.entries.view_host()(c), face_cell_ids_.entries.view_host()(pos));
          ++pos;
        }
      }
    }
    assert(pos <= end);
  }

  // Copy to the dirs
  Kokkos::resize(face_cell_ids_dirs_, face_cell_ids_.entries.extent(0));
  for (int i = 0; i < face_cell_ids_.entries.extent(0); ++i) {
    if (face_cell_ids_.entries.view_host()(i) < 0)
      face_cell_ids_dirs_.view_host()(i) = ~face_cell_ids_.entries.view_host()(i);
    else
      face_cell_ids_dirs_.view_host()(i) = face_cell_ids_.entries.view_host()(i);
  }

  cell2face_info_cached_ = true;
  face2cell_info_cached_ = true;
  faces_requested_ = true;
}

void
MeshCache::display_cache()
{
  std::cout << std::endl;
  if (!cell2face_info_cached_) cache_cell_face_info_();
  if (!cell_geometry_precomputed_) compute_cell_geometric_quantities_();
  if (!face_geometry_precomputed_) compute_face_geometric_quantities_();
  // Display all cache
  std::cout << "face_centroid: " << std::endl;
  for (int i = 0; i < face_centroids_.extent(0); ++i)
    std::cout << face_centroids_.view_host()(i)[0] << " " << face_centroids_.view_host()(i)[1] << " "
              << face_centroids_.view_host()(i)[2] << std::endl;
  std::cout << std::endl;

  std::cout << "face_normal: " << std::endl;
  for (int i = 0; i < face_normals_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = face_normals_.row_map.view_host()(i); j < face_normals_.row_map.view_host()(i + 1);
         ++j) {
      std::cout << face_normals_.entries.view_host()(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "face_cell_ptype_: " << std::endl;
  for (int i = 0; i < face_cell_ptype_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = face_cell_ptype_.row_map.view_host()(i);
         j < face_cell_ptype_.row_map.view_host()(i + 1);
         ++j) {
      std::cout << static_cast<int>(face_cell_ptype_.entries.view_host()(j)) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;


  std::cout << "face_cell_ids_: " << std::endl;
  for (int i = 0; i < face_cell_ids_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = face_cell_ids_.row_map.view_host()(i); j < face_cell_ids_.row_map.view_host()(i + 1);
         ++j) {
      std::cout << face_cell_ids_.entries.view_host()(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "cell_face_ids_: " << std::endl;
  for (int i = 0; i < cell_face_ids_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = cell_face_ids_.row_map.view_host()(i); j < cell_face_ids_.row_map.view_host()(i + 1);
         ++j) {
      std::cout << cell_face_ids_.entries.view_host()(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "cell_face_dirs_: " << std::endl;
  for (int i = 0; i < cell_face_dirs_.row_map.extent(0) - 1; ++i) {
    std::cout << i << ": ";
    for (int j = cell_face_dirs_.row_map.view_host()(i); j < cell_face_dirs_.row_map.view_host()(i + 1);
         ++j) {
      std::cout << cell_face_dirs_.entries.view_host()(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

// Gather and cache face to edge connectivity info.
void
MeshCache::cache_face2edge_info_() const
{
  int nfaces = mesh_->num_entities(FACE, Parallel_type::ALL);
  face_edge_ids_.row_map.resize(nfaces + 1);
  face_edge_dirs_.row_map.resize(nfaces + 1);
  face_edge_ids_.row_map.view_host()(0) = 0;
  face_edge_dirs_.row_map.view_host()(0) = 0;

  int face_edge_ids_size = 0; 
  for (int f = 0; f < nfaces; f++) {
    Entity_ID_List fedgeids;
    std::vector<int> fedgedirs;
    mesh_->face_get_edges_and_dirs_internal_(f, fedgeids, &fedgedirs, true);
    face_edge_ids_.row_map.view_host()(f + 1) =
      face_edge_ids_.row_map.view_host()(f) + fedgeids.size();
    face_edge_dirs_.row_map.view_host()(f + 1) =
      face_edge_dirs_.row_map.view_host()(f) + fedgedirs.size();
    face_edge_ids_size += fedgeids.size();
  }

  face_edge_ids_.entries.resize(face_edge_ids_size); 
  face_edge_dirs_.entries.resize(face_edge_ids_size); 

  for (int f = 0; f < nfaces; f++) {
    Entity_ID_List fedgeids;
    std::vector<int> fedgedirs;
    mesh_->face_get_edges_and_dirs_internal_(f, fedgeids, &fedgedirs, true);
    face_edge_ids_.row_map.view_host()(f + 1) =
      face_edge_ids_.row_map.view_host()(f) + fedgeids.size();
    face_edge_dirs_.row_map.view_host()(f + 1) =
      face_edge_dirs_.row_map.view_host()(f) + fedgedirs.size();
    for (int i = 0; i < fedgeids.size(); ++i) {
      face_edge_ids_.entries.view_host()(face_edge_ids_.row_map.view_host()(f) + i) = fedgeids[i];
      face_edge_dirs_.entries.view_host()(face_edge_dirs_.row_map.view_host()(f) + i) = fedgedirs[i];
    }
  }

  face2edge_info_cached_ = true;
  faces_requested_ = true;
  edges_requested_ = true;
}


// Gather and cache cell to edge connectivity info.
void
MeshCache::cache_cell2edge_info_() const
{
  int ncells = mesh_->num_entities(CELL, Parallel_type::ALL);
  Kokkos::resize(cell_edge_ids_.row_map, ncells + 1);

  if (mesh_->space_dimension() == 2) {
    Kokkos::resize(cell_2D_edge_dirs_.row_map, ncells + 1);
    int cell_edge_ids_entries_size = 0; 
    int cell_2D_edge_dirs_entries_size = 0; 
    for (int c = 0; c < ncells; c++) {
      Entity_ID_List cell_edge_ids_tmp;
      std::vector<int> cell_2D_edge_dirs_tmp;
      mesh_->cell_2D_get_edges_and_dirs_internal_(
        c, cell_edge_ids_tmp, &cell_2D_edge_dirs_tmp);
      cell_edge_ids_.row_map.view_host()(c + 1) =
        cell_edge_ids_.row_map.view_host()(c) + cell_edge_ids_tmp.size();
      cell_2D_edge_dirs_.row_map.view_host()(c + 1) =
        cell_2D_edge_dirs_.row_map.view_host()(c) + cell_2D_edge_dirs_tmp.size();
      cell_edge_ids_entries_size += cell_edge_ids_tmp.size(); 
      cell_2D_edge_dirs_entries_size += cell_2D_edge_dirs_tmp.size();
    }
    Kokkos::resize(cell_edge_ids_.entries,cell_edge_ids_entries_size); 
    Kokkos::resize(cell_2D_edge_dirs_.entries,cell_2D_edge_dirs_entries_size);
    for (int c = 0; c < ncells; c++) {
      Entity_ID_List cell_edge_ids_tmp;
      std::vector<int> cell_2D_edge_dirs_tmp;
      mesh_->cell_2D_get_edges_and_dirs_internal_(
        c, cell_edge_ids_tmp, &cell_2D_edge_dirs_tmp);
      cell_edge_ids_.row_map.view_host()(c + 1) =
        cell_edge_ids_.row_map.view_host()(c) + cell_edge_ids_tmp.size();
      cell_2D_edge_dirs_.row_map.view_host()(c + 1) =
        cell_2D_edge_dirs_.row_map.view_host()(c) + cell_2D_edge_dirs_tmp.size();
      for (int i = 0; i < cell_edge_ids_tmp.size(); ++i) {
        cell_edge_ids_.entries.view_host()(cell_edge_ids_.row_map.view_host()(c) + i) =
          cell_edge_ids_tmp[i];
        cell_2D_edge_dirs_.entries.view_host()(cell_2D_edge_dirs_.row_map.view_host()(c) + i) =
          cell_2D_edge_dirs_tmp[i];
      }
    }
  } else  {
    int cell_edge_ids_size = 0; 
    for (int c = 0; c < ncells; c++) {
      Entity_ID_List cell_edge_ids_tmp;
      mesh_->cell_get_edges_internal_(c, cell_edge_ids_tmp);
      cell_edge_ids_.row_map.view_host()(c + 1) =
        cell_edge_ids_.row_map.view_host()(c) + cell_edge_ids_tmp.size();
      cell_edge_ids_size += cell_edge_ids_tmp.size();
    }
    Kokkos::resize(cell_edge_ids_.entries,cell_edge_ids_size); 
    for (int c = 0; c < ncells; c++) {
      Entity_ID_List cell_edge_ids_tmp;
      mesh_->cell_get_edges_internal_(c, cell_edge_ids_tmp);
      cell_edge_ids_.row_map.view_host()(c + 1) =
        cell_edge_ids_.row_map.view_host()(c) + cell_edge_ids_tmp.size();
      for (int i = 0; i < cell_edge_ids_tmp.size(); ++i) {
        cell_edge_ids_.entries.view_host()(cell_edge_ids_.row_map.view_host()(c) + i) =
          cell_edge_ids_tmp[i];
      }
    }
  }
  cell2edge_info_cached_ = true;
}

void 
MeshCache::cache_parents_info_() const 
{
  // CELLS
  int ncells = mesh_->num_entities(CELL, Parallel_type::ALL);
  Kokkos::resize(cells_parent_,ncells); 
  for(int i = 0 ; i < ncells; ++i){
    cells_parent_.view_host()[i] = mesh_->entity_get_parent_type(CELL,i); 
  }
  // FACES 
  int nfaces = mesh_->num_entities(FACE, Parallel_type::ALL);
  Kokkos::resize(faces_parent_,nfaces);
  for(int i = 0 ; i < nfaces; ++i){
    faces_parent_.view_host()[i] = mesh_->entity_get_parent_type(FACE,i); 
  }
  // NODES 
  int nnodes = mesh_->num_entities(NODE, Parallel_type::ALL);
  Kokkos::resize(nodes_parent_,nnodes);
  for(int i = 0 ; i < nnodes; ++i){
    nodes_parent_.view_host()[i] = mesh_->entity_get_parent_type(NODE,i); 
  }
  // EDGES
  if(edges_requested_){
    int nedges = mesh_->num_entities(EDGE, Parallel_type::ALL);
    Kokkos::resize(edges_parent_,nedges);
    for(int i = 0 ; i < nedges; ++i){
      edges_parent_.view_host()[i] = mesh_->entity_get_parent_type(EDGE,i); 
    }
  }
  parents_precomputed_ = true; 
}

void
MeshCache::init()
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
  if(mesh_->parent().get()){
    assert(!parents_precomputed_); 
    cache_parents_info_(); 
    parents_precomputed_ = true; 
  }
  assert(!cell_get_faces_and_bisectors_precomputed_); 
  cache_cell_get_faces_and_bisectors_(); 
  cell_get_faces_and_bisectors_precomputed_ = true;

  // Sync GPU 
  // Dual views 
  cell_volumes_.sync<Kokkos::DefaultExecutionSpace>();
  face_areas_.sync<Kokkos::DefaultExecutionSpace>();
  edge_lengths_.sync<Kokkos::DefaultExecutionSpace>();
  cell_centroids_.sync<Kokkos::DefaultExecutionSpace>();
  face_centroids_.sync<Kokkos::DefaultExecutionSpace>();
  edge_vectors_.sync<Kokkos::DefaultExecutionSpace>();
  cells_parent_.sync<Kokkos::DefaultExecutionSpace>();
  faces_parent_.sync<Kokkos::DefaultExecutionSpace>();
  edges_parent_.sync<Kokkos::DefaultExecutionSpace>();
  nodes_parent_.sync<Kokkos::DefaultExecutionSpace>();
  cell_cellabove_.sync<Kokkos::DefaultExecutionSpace>();
  cell_cellbelow_.sync<Kokkos::DefaultExecutionSpace>();
  node_nodeabove_.sync<Kokkos::DefaultExecutionSpace>();
  columnID_.sync<Kokkos::DefaultExecutionSpace>();
  face_cell_ids_dirs_.sync<Kokkos::DefaultExecutionSpace>();

  
  // Dual CRS 
  face_normals_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  face_normals_.entries.sync<Kokkos::DefaultExecutionSpace>();
  cell_faces_bisectors_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  cell_faces_bisectors_.entries.sync<Kokkos::DefaultExecutionSpace>();
  cell_face_ids_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  cell_face_ids_.entries.sync<Kokkos::DefaultExecutionSpace>();
  cell_face_dirs_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  cell_face_dirs_.entries.sync<Kokkos::DefaultExecutionSpace>();
  face_cell_ids_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  face_cell_ids_.entries.sync<Kokkos::DefaultExecutionSpace>();
  face_cell_ids_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  face_cell_ids_.entries.sync<Kokkos::DefaultExecutionSpace>();
  face_cell_ptype_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  face_cell_ptype_.entries.sync<Kokkos::DefaultExecutionSpace>();
  cell_edge_ids_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  cell_edge_ids_.entries.sync<Kokkos::DefaultExecutionSpace>();
  cell_2D_edge_dirs_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  cell_2D_edge_dirs_.entries.sync<Kokkos::DefaultExecutionSpace>();
  face_edge_ids_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  face_edge_ids_.entries.sync<Kokkos::DefaultExecutionSpace>();
  face_edge_dirs_.row_map.sync<Kokkos::DefaultExecutionSpace>();
  face_edge_dirs_.entries.sync<Kokkos::DefaultExecutionSpace>();
}

void 
MeshCache::cache_cell_get_faces_and_bisectors_() const {

  int ncells = mesh_->num_entities(CELL, Parallel_type::ALL);
  Kokkos::resize(cell_faces_bisectors_.row_map,ncells+1); 
  cell_faces_bisectors_.row_map.view_host()(0) = 0; 

  int entries_size = 0; 
  for(int i = 0 ; i < ncells; ++i){
    Kokkos::View<Entity_ID*,Kokkos::HostSpace> faceids;
    cell_get_faces_host(i, faceids);
    entries_size += faceids.size(); 
  }

  cell_faces_bisectors_.entries.resize(entries_size); 
  for(int i = 0 ; i < ncells; ++i){
    Kokkos::View<Entity_ID*,Kokkos::HostSpace> faceids;
    cell_get_faces_host(i, faceids);
    AmanziGeometry::Point cc = cell_centroid_host(i);
    cell_faces_bisectors_.row_map.view_host()(i + 1) = faceids.extent(0) +
       cell_faces_bisectors_.row_map.view_host()(i);
    for (int j = 0; j < faceids.extent(0); ++j) {
      cell_faces_bisectors_.entries.view_host()(cell_faces_bisectors_.row_map.view_host()(i)+j) = 
        face_centroid_host(faceids(j)) - cc;
    }
  }
  cell_faces_bisectors_.entries.sync_device();
  cell_faces_bisectors_.row_map.sync_device();

}

int
MeshCache::compute_cell_geometric_quantities_() const
{
  int ncells = mesh_->num_entities(CELL, Parallel_type::ALL);

  cell_volumes_.resize(ncells);
  cell_centroids_.resize(ncells);
  for (int i = 0; i < ncells; i++) {
    double volume;
    AmanziGeometry::Point centroid(mesh_->space_dimension());

    compute_cell_geometry_(i, &volume, &centroid);

    cell_volumes_.view_host()(i) = volume;
    cell_centroids_.view_host()(i) = centroid;
  }

  cell_geometry_precomputed_ = true;

  cell_volumes_.sync_device(); 
  cell_centroids_.sync_device(); 

  return 1;
}



int
MeshCache::compute_face_geometric_quantities_() const
{
  if (mesh_->space_dimension() == 3 && mesh_->manifold_dimension() == 2) {
    // need cell centroids to compute normals
    if (!cell_geometry_precomputed_) compute_cell_geometric_quantities_();
  }

  int nfaces = mesh_->num_entities(FACE, Parallel_type::ALL);

  face_areas_.resize(nfaces);
  face_centroids_.resize(nfaces);

  // Temporary views
  face_normals_.row_map.resize(nfaces + 1);
  face_normals_.row_map.view_host()(0) = 0;

  // Find size 
  int entries_size = 0; 
  for (int i = 0; i < nfaces; i++) {
    double area;
    AmanziGeometry::Point centroid(mesh_->space_dimension());
    std::vector<AmanziGeometry::Point> normals;
    // normal0 and normal1 are outward normals of the face with
    // respect to the cell0 and cell1 of the face. The natural normal
    // of the face points out of cell0 and into cell1. If one of these
    // cells do not exist, then the normal is the null vector.
    compute_face_geometry_(i, &area, &centroid, normals);
    face_areas_.view_host()(i) = area;
    face_centroids_.view_host()(i) = centroid;
    entries_size += normals.size(); 
  }

  face_normals_.entries.resize(entries_size); 

  for (int i = 0; i < nfaces; i++) {
    double area;
    AmanziGeometry::Point centroid(mesh_->space_dimension());
    std::vector<AmanziGeometry::Point> normals;
    compute_face_geometry_(i, &area, &centroid, normals);
    for (int j = 0; j < normals.size(); ++j) {
      face_normals_.entries.view_host()(face_normals_.row_map.view_host()(i) + j) = normals[j];
    }
    face_normals_.row_map.view_host()(i + 1) = normals.size() + face_normals_.row_map.view_host()(i);
  }

  face_geometry_precomputed_ = true;

  face_normals_.row_map.sync_device(); 
  face_areas_.sync_device(); 
  face_centroids_.sync_device(); 
  face_normals_.entries.sync_device(); 

  return 1;
}


int
MeshCache::compute_edge_geometric_quantities_() const
{
  int nedges = mesh_->num_entities(EDGE, Parallel_type::ALL);

  edge_vectors_.resize(nedges);
  edge_lengths_.resize(nedges);

  for (int i = 0; i < nedges; i++) {
    double length;
    AmanziGeometry::Point evector(mesh_->space_dimension());

    compute_edge_geometry_(i, &length, &evector);

    edge_lengths_.view_host()(i) = length;
    edge_vectors_.view_host()(i) = evector;
  }

  edge_geometry_precomputed_ = true;
  
  edge_lengths_.sync_device(); 
  edge_vectors_.sync_device(); 
  
  return 1;
}


int
MeshCache::compute_cell_geometry_(const Entity_ID cellid, double* volume,
                             AmanziGeometry::Point* centroid) const
{
  if (mesh_->manifold_dimension() == 3) {
    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine

    // General polyhedra always need to have an explicit face
    // representation - special elements like hexes can get away
    // without (but we have yet to put in the code for the standard
    // node ordering and computation for these special elements)
    Kokkos::View<Entity_ID*,Kokkos::HostSpace> faces;
    std::vector<unsigned int> nfnodes;
    Kokkos::View<int*,Kokkos::HostSpace> fdirs;
    std::vector<AmanziGeometry::Point> fcoords, ccoords, cfcoords; 

    cell_get_faces_and_dirs_host(cellid, faces, fdirs);

    int nf = faces.extent(0);
    nfnodes.resize(nf);

    if (nf == 2) {
      /* special case of column mesh - only top and bottom faces
         are returned */
      AmanziGeometry::Point fcentroid0(mesh_->space_dimension()), fcentroid1(mesh_->space_dimension());
      AmanziGeometry::Point normal(mesh_->space_dimension());
      double farea;

      /* compute volume on the assumption that the top and bottom faces form
         a vertical columnar cell or in other words a polygonal prism */

      mesh_->face_get_coordinates(faces(0), fcoords);
      AmanziGeometry::polygon_get_area_centroid_normal(
        fcoords, &farea, &fcentroid0, &normal);

      mesh_->face_get_coordinates(faces(1), fcoords);
      AmanziGeometry::polygon_get_area_centroid_normal(
        fcoords, &farea, &fcentroid1, &normal);

      *centroid = (fcentroid0 + fcentroid1) / 2.0;
      double height = norm(fcentroid1 - fcentroid0);

      *volume = farea * height;
    } else { /* general case */

      size_t cfcoords_size = 0;
      for (int j = 0; j < nf; j++) {
        mesh_->face_get_coordinates(faces(j), fcoords);
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

      mesh_->cell_get_coordinates(cellid, ccoords);

      AmanziGeometry::polyhed_get_vol_centroid(
        ccoords, nf, nfnodes, cfcoords, volume, centroid);
    }
    return 1;
  } else if (mesh_->manifold_dimension() == 2) {
    std::vector<AmanziGeometry::Point> ccoords;
    mesh_->cell_get_coordinates(cellid, ccoords);

    AmanziGeometry::Point normal(mesh_->space_dimension());

    AmanziGeometry::polygon_get_area_centroid_normal(
      ccoords, volume, centroid, &normal);
    return 1;
  }

  return 0;
}


int
MeshCache::compute_face_geometry_(
  const Entity_ID faceid, double* area, AmanziGeometry::Point* centroid,
  std::vector<AmanziGeometry::Point>& normals) const
{
  std::vector<AmanziGeometry::Point> fcoords;
  // normals->clear();

  if (mesh_->manifold_dimension() == 3) {
    // 3D Elements with possibly curved faces

    mesh_->face_get_coordinates(faceid, fcoords);

    AmanziGeometry::Point normal(3);
    AmanziGeometry::polygon_get_area_centroid_normal(
      fcoords, area, centroid, &normal);

    Kokkos::View<Entity_ID*,Kokkos::HostSpace> cellids;
    face_get_cells_host(faceid, Parallel_type::ALL, cellids);
    AMANZI_ASSERT(cellids.extent(0) <= 2);

    normals.resize(cellids.size()); 
    for (int i = 0; i < cellids.extent(0); i++) {
      Kokkos::View<Entity_ID*,Kokkos::HostSpace> cellfaceids;
      Kokkos::View<int*,Kokkos::HostSpace> cellfacedirs;
      int dir = 1;

      cell_get_faces_and_dirs_host(cellids(i), cellfaceids, cellfacedirs);

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
  } else if (mesh_->manifold_dimension() == 2) {
    if (mesh_->space_dimension() == 2) { // 2D mesh

      mesh_->face_get_coordinates(faceid, fcoords);

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

      mesh_->face_get_coordinates(faceid, fcoords);

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

        AmanziGeometry::Point cvec = fcoords[0] - cell_centroids_.view_host()(cellids(i));
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
MeshCache::compute_edge_geometry_(const Entity_ID edgeid, double* edge_length,
                             AmanziGeometry::Point* edge_vector) const
{
  (*edge_vector).set(0.0L);
  *edge_length = 0.0;

  Entity_ID node0, node1;

  mesh_->edge_get_nodes(edgeid, &node0, &node1);

  AmanziGeometry::Point point0, point1;
  mesh_->node_get_coordinates(node0, &point0);
  mesh_->node_get_coordinates(node1, &point1);

  *edge_vector = point1 - point0;
  *edge_length = norm(*edge_vector);

  return 0;
}

// Volume/Area of cell
double
MeshCache::cell_volume(const Entity_ID cellid, const bool recompute) const
{
  if (!cell_geometry_precomputed_) {
    compute_cell_geometric_quantities_();
    return cell_volumes_.view_host()(cellid);
  } else {
    if (recompute) {
      double volume;
      AmanziGeometry::Point centroid(mesh_->space_dimension());
      compute_cell_geometry_(cellid, &volume, &centroid);
      return volume;
    } else
      return cell_volumes_.view_host()(cellid);
  }
}

// Length of an edge
double
MeshCache::edge_length(const Entity_ID edgeid, const bool recompute) const
{
  AMANZI_ASSERT(edges_requested_);

  if (!edge_geometry_precomputed_) {
    compute_edge_geometric_quantities_();
    return edge_lengths_.view_host()(edgeid);
  } else {
    if (recompute) {
      double length;
      AmanziGeometry::Point vector(mesh_->space_dimension());
      compute_edge_geometry_(edgeid, &length, &vector);
      return length;
    } else
      return edge_lengths_.view_host()(edgeid);
  }
}

// Centroid of edge
AmanziGeometry::Point
MeshCache::edge_centroid(const Entity_ID edgeid) const
{
  Entity_ID p0, p1;
  AmanziGeometry::Point xyz0, xyz1;

  mesh_->edge_get_nodes(edgeid, &p0, &p1);
  mesh_->node_get_coordinates(p0, &xyz0);
  mesh_->node_get_coordinates(p1, &xyz1);
  return (xyz0 + xyz1) / 2;
}


// Direction vector of edge
AmanziGeometry::Point
MeshCache::edge_vector(const Entity_ID edgeid, const bool recompute,
                  const Entity_ID pointid, int* orientation) const
{
  AMANZI_ASSERT(edges_requested_);

  AmanziGeometry::Point evector(mesh_->space_dimension());
  AmanziGeometry::Point& evector_ref = evector; // to avoid extra copying

  if (!edge_geometry_precomputed_) compute_edge_geometric_quantities_();

  if (recompute) {
    double length;
    compute_edge_geometry_(edgeid, &length, &evector);
    // evector_ref already points to evector
  } else
    evector_ref = edge_vectors_.view_host()(edgeid);

  if (orientation) *orientation = 1;

  if (pointid == -1)
    return evector_ref;
  else {
    Entity_ID p0, p1;
    mesh_->edge_get_nodes(edgeid, &p0, &p1);

    if (pointid == p0)
      return evector_ref;
    else {
      if (orientation) *orientation = -1;
      return -evector_ref;
    }
  }
}


bool
MeshCache::point_in_cell(const AmanziGeometry::Point& p,
                    const Entity_ID cellid) const
{
  std::vector<AmanziGeometry::Point> ccoords;

  if (mesh_->manifold_dimension() == 3) {
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
      mesh_->face_get_coordinates(faces(j), fcoords);
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

    mesh_->cell_get_coordinates(cellid, ccoords);

    return AmanziGeometry::point_in_polyhed(p, ccoords, nf, nfnodes, cfcoords);

  } else if (mesh_->manifold_dimension() == 2) {
    mesh_->cell_get_coordinates(cellid, ccoords);
    return AmanziGeometry::point_in_polygon(p, ccoords);
  }

  return false;
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
MeshCache::build_columns(const std::string& setname) const
{
  if (columns_built_) return 0;

  // Allocate space and initialize.
  int nn = mesh_->num_entities(NODE, Parallel_type::ALL);
  int nf = mesh_->num_entities(FACE, Parallel_type::ALL);
  int nf_owned = mesh_->num_entities(FACE, Parallel_type::OWNED);
  int nc = mesh_->num_entities(CELL, Parallel_type::ALL);
  int nc_owned = mesh_->num_entities(CELL, Parallel_type::OWNED);

  Kokkos::resize(columnID_, nc);
  Kokkos::resize(cell_cellbelow_, nc);
  Kokkos::resize(cell_cellabove_, nc);
  Kokkos::resize(node_nodeabove_, nn);

  Kokkos::parallel_for(
    "Mesh::build_columns loop 1",
    nc, KOKKOS_LAMBDA(const int& i) {
    cell_cellbelow_.view_host()(i) = -1;
    cell_cellabove_.view_host()(i) = -1;
  });
  Kokkos::parallel_for(
    "Mesh::build_columns loop 2",
    nn, KOKKOS_LAMBDA(const int& i) { node_nodeabove_.view_host()(i) = -1; });

  Entity_ID_List top_faces;
  mesh_->get_set_entities(setname, FACE, Parallel_type::ALL, top_faces);

  int ncolumns = top_faces.size();
  num_owned_cols_ = mesh_->get_set_size(setname, FACE, Parallel_type::OWNED);

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
    *(mesh_->get_comm()), Teuchos::REDUCE_MIN, 1, &success, &min_success);
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
MeshCache::build_columns() const
{
  if (columns_built_) return 1;

  // Allocate space and initialize.
  int nn = mesh_->num_entities(NODE, Parallel_type::ALL);
  int nf = mesh_->num_entities(FACE, Parallel_type::ALL);
  int nf_owned = mesh_->num_entities(FACE, Parallel_type::OWNED);
  int nc = mesh_->num_entities(CELL, Parallel_type::ALL);
  int nc_owned = mesh_->num_entities(CELL, Parallel_type::OWNED);

  Kokkos::resize(columnID_, nc);
  Kokkos::resize(cell_cellbelow_, nc);
  Kokkos::resize(cell_cellabove_, nc);
  Kokkos::resize(node_nodeabove_, nn);

  Kokkos::parallel_for(
    "Mesh::build_columns loop 1",
    nc, KOKKOS_LAMBDA(const int& i) {
    cell_cellbelow_.view_host()(i) = -1;
    cell_cellabove_.view_host()(i) = -1;
  });
  Kokkos::parallel_for(
    "Mesh::build_columns loop 2",
    nn, KOKKOS_LAMBDA(const int& i) { node_nodeabove_.view_host()(i) = -1; });

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

    AmanziGeometry::Point zvec(mesh_->space_dimension());
    if (mesh_->space_dimension() == 2)
      zvec.set(0.0, 1.0);
    else if (mesh_->space_dimension() == 3)
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
    *(mesh_->get_comm()), Teuchos::REDUCE_MIN, 1, &success, &min_success);
  columns_built_ = (min_success == 1);
  return columns_built_ ? 1 : 0;
}


int
MeshCache::build_single_column_(int colnum, Entity_ID top_face) const
{
  Kokkos::View<Entity_ID*> fcells;
  face_get_cells(top_face, Parallel_type::ALL, fcells);

  // Walk through the cells until we get to the bottom of the domain
  Entity_ID cur_cell = fcells(0);
  bool is_ghost_column =
    (mesh_->entity_get_ptype(CELL, cur_cell) == Parallel_type::GHOST);
  Entity_ID bot_face = -1;
  Entity_ID_List colcells, colfaces;
  Kokkos::View<Entity_ID*> cfaces, fcells2;
  Kokkos::View<int*> cfdirs;

  AmanziGeometry::Point negzvec(mesh_->space_dimension());
  if (mesh_->space_dimension() == 2)
    negzvec.set(0.0, -1.0);
  else if (mesh_->space_dimension() == 3)
    negzvec.set(0.0, 0.0, -1.0);

  int success = 1;
  bool done = false;
  while (!done) {
    bool is_ghost_cell =
      (mesh_->entity_get_ptype(CELL, cur_cell) == Parallel_type::GHOST);
    if (is_ghost_column != is_ghost_cell) {
      //      Errors::Message mesg("A column contains cells from different mesh
      //      partitions!"); Exceptions::amanzi_throw(mesg);
      std::cerr << "A column contains cells from different mesh partitions"
                << std::endl;
      success = 0;
      break;
    }
    columnID_.view_host()(cur_cell) = colnum;
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
      if (cell_cellbelow_.view_host()(cur_cell) != -1) { // intersecting column of cells
        std::cerr << "Intersecting column of cells\n";
        success = 0;
        break;
      }

      if (fcells2(0) == cur_cell) {
        cell_cellbelow_.view_host()(cur_cell) = fcells2(1);
        cell_cellabove_.view_host()(fcells2(1)) = cur_cell;
        cur_cell = fcells2(1);
      } else if (fcells2(1) == cur_cell) {
        cell_cellbelow_.view_host()(cur_cell) = fcells2(0);
        cell_cellabove_.view_host()(fcells2(0)) = cur_cell;
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

    mesh_->face_get_nodes(top_face, topnodes);
    Entity_ID topnode0 = topnodes[0];

    // nodes of the top face
    mesh_->face_get_nodes(bot_face, botnodes);

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
    mesh_->node_get_coordinates(topnode0, &topnode0c);

    for (int k = 0; k < nfvtop; k++) {
      AmanziGeometry::Point kc;
      mesh_->node_get_coordinates(botnodes[k], &kc);

      double horiz_dist = 0.;
      for (int m = 0; m != mesh_->space_dimension() - 1; ++m) {
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
    double even_odd_product = face_normal(bot_face)[mesh_->space_dimension() - 1] *
                              face_normal(top_face)[mesh_->space_dimension() - 1];
    AMANZI_ASSERT(std::abs(even_odd_product) > 0);
    int even_odd = even_odd_product >= 0. ? 1 : -1;

    for (int k = 0; k < nfvtop; k++) {
      Entity_ID topnode = topnodes[k];
      int bot_i = (ind + even_odd * k) % nfvtop;
      if (bot_i < 0) bot_i += nfvtop;
      Entity_ID botnode = botnodes[bot_i];
      node_nodeabove_.view_host()(botnode) = topnode;
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

} // namespace AmanziMesh
} // namespace Amanzi
