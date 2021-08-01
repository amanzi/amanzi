#ifndef AMANZI_MESH_CACHE_HH_
#define AMANZI_MESH_CACHE_HH_

#include <map>
#include <vector>
#include <string>

#include "nanoflann.hpp"

#include "errors.hh"
#include "AmanziTypes.hh"
#include "VerboseObject.hh"
#include "Key.hh"

#include "GeometricModel.hh"
#include "Point.hh"
#include "Region.hh"

#include "CellTopology.hh"
#include "KDTree.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

class Mesh; 

class MeshCache {

 public:
  MeshCache(Mesh* mesh, const bool request_faces, const bool request_edges);

  ~MeshCache() = default;

  //
  // Initialization of the mesh cache
  // --------------------------------
  // Calls all the function to initialize the mesh cache
  // this is called in the meshfactory once the mesh is created
  void init();


  // Downward Adjacencies
  //---------------------

  // Get number of faces of a cell.
  //
  // On a distributed mesh, this will count all the faces of the
  // cell, OWNED or GHOST.
  KOKKOS_INLINE_FUNCTION unsigned int cell_get_num_faces(const Entity_ID cellid) const{
      assert(cell2face_info_cached_);
      return cell_face_ids_.row_map(cellid + 1) - cell_face_ids_.row_map(cellid);
  }

  // Get faces of a cell.
  //
  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order
  //
  // EXTENSIONS: MSTK FRAMEWORK: by the way the parallel partitioning,
  // send-receive protocols and mesh query operators are designed, a side
  // effect of this is that master and ghost entities will have the same
  // hierarchical topology.
  KOKKOS_INLINE_FUNCTION void
  cell_get_faces(const Entity_ID cellid,
                 Kokkos::View<Entity_ID*>& faceids) const
  {
    assert(cell2face_info_cached_);
    faceids =
      Kokkos::subview(cell_face_ids_.entries,
                      Kokkos::make_pair(cell_face_ids_.row_map(cellid),
                                        cell_face_ids_.row_map(cellid + 1)));
  }

  // Get faces of a cell and directions in which the cell uses the face
  //
  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order
  //
  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise
  KOKKOS_INLINE_FUNCTION void
  cell_get_faces_and_dirs(const Entity_ID cellid,
                          Kokkos::View<Entity_ID*>& faceids,
                          Kokkos::View<int*>& face_dirs) const
  {
    assert(cell2face_info_cached_);
    faceids =
      Kokkos::subview(cell_face_ids_.entries,
                      Kokkos::make_pair(cell_face_ids_.row_map(cellid),
                                        cell_face_ids_.row_map(cellid + 1)));
    face_dirs =
      Kokkos::subview(cell_face_dirs_.entries,
                      Kokkos::make_pair(cell_face_dirs_.row_map(cellid),
                                        cell_face_dirs_.row_map(cellid + 1)));
  }
  
  // Get the bisectors, i.e. vectors from cell centroid to face centroids.
  KOKKOS_INLINE_FUNCTION void 
  cell_get_faces_and_bisectors(const Entity_ID cellid, 
                      Kokkos::View<Entity_ID*>& faceids,
                      Kokkos::View<AmanziGeometry::Point*>& bisectors) const
  {
    assert(cell_get_faces_and_bisectors_precomputed_); 
    faceids =
      Kokkos::subview(cell_face_ids_.entries,
                      Kokkos::make_pair(cell_face_ids_.row_map(cellid),
                                        cell_face_ids_.row_map(cellid + 1)));
    bisectors = 
      Kokkos::subview(cell_faces_bisectors_.entries,
                    Kokkos::make_pair(cell_faces_bisectors_.row_map(cellid),
                                      cell_faces_bisectors_.row_map(cellid + 1)));
  }

  // Cells connected to a face
  //
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors
  KOKKOS_INLINE_FUNCTION void
  face_get_cells(const Entity_ID faceid, const Parallel_type ptype,
                 Kokkos::View<Entity_ID*>& cellids) const
  {
    int ncellids = 0;
    int csr_start = face_cell_ptype_.row_map(faceid);
    int csr_end = face_cell_ptype_.row_map(faceid + 1);

    // Create a subview
    int start = -1;
    int stop = -1;
    // Return all cells except the UNKNOWN
    if (ptype == Parallel_type::ALL) {
      bool done = false;
      int i = csr_start;
      stop = csr_end;
      while (!done && i < csr_end) {
        Parallel_type cell_ptype = face_cell_ptype_.entries(i);
        if (cell_ptype == Parallel_type::PTYPE_UNKNOWN) {
          ++i;
          continue;
        } else {
          start = i;
          done = true;
        }
      }
    } else if (ptype == Parallel_type::PTYPE_UNKNOWN) {
      assert(false);
      return;
    } else {
      bool done = false;
      int i = csr_start;
      while (!done && i < csr_end) {
        Parallel_type cell_ptype = face_cell_ptype_.entries(i);
        if (cell_ptype == Parallel_type::PTYPE_UNKNOWN) {
          ++i;
          continue;
        }
        if (cell_ptype == ptype && start == -1) { start = i; }
        if (cell_ptype != ptype && start != -1) { stop = i; }
        ++i;
      }
      if (stop == -1) { stop = i; }
    }
    // Generate the subview
    cellids =
      Kokkos::subview(face_cell_ids_dirs_, Kokkos::make_pair(start, stop));
  }


  //
  // Column information
  //-------------------

  // Special adjacency information for geological domains with a
  // semi-structured mesh. Will return -1 if there is no suitable cell
  // The code currently makes the assumption that the "bottom" of the
  // is a flat surface in the XY plane. It then builds up information
  // about the cell above and cell below for each cell based on the
  // orientation of the face normals w.r.t the z direction. If the
  // mesh is highly warped, this could lead to ambiguities. Also,
  // intersecting columns in an unstructured mesh will lead to an
  // exception being thrown. These data structures are never populated
  // if these operators are never called. The above and below cells
  // are computed for all cells the first time one of these routines
  // is called and then cached

  //
  int build_columns(const std::string& Msetname) const;

  // Build columns over the entire mesh. The columns are defined by
  // starting from boundary faces which have a negative-z-direction
  // normal, then collecting cells and faces while traveling downward
  // through the columns.

  int build_columns() const;

  // Parent entity in the source mesh if mesh was derived from another mesh
  KOKKOS_INLINE_FUNCTION Entity_ID
  entity_get_parent(const Entity_kind kind, const Entity_ID entid) const
  {
    assert(parents_precomputed_); 
    switch(kind){
      case CELL:
        if(entid > cells_parent_.extent(0))
          return -1; 
        return cells_parent_(entid);  
        break;
      case FACE: 
        if(entid > faces_parent_.extent(0))
          return -1; 
        return faces_parent_(entid);
        break; 
      case EDGE: 
        assert(edges_requested_); 
        if(entid > edges_parent_.extent(0))
          return -1; 
        return edges_parent_(entid);
        break; 
      case NODE: 
        if(entid > nodes_parent_.extent(0))
          return -1; 
        return nodes_parent_(entid);
      default: {} 
    }
    return -1;
  }

  //
  // Mesh entity geometry
  // --------------------

  // Volume/Area of cell
  double cell_volume(const Entity_ID cellid, const bool recompute) const;
  KOKKOS_INLINE_FUNCTION double cell_volume(const Entity_ID cellid) const
  {
    return cell_volumes_(cellid);
  }

  // Area/length of face
  KOKKOS_INLINE_FUNCTION double
  face_area(const Entity_ID faceid, const bool recompute = false) const
  {
    return face_areas_(faceid);
  }

  // Length of edge
  double
  edge_length(const Entity_ID edgeid, const bool recompute = false) const;

  // Centroid of cell (center of gravity not just average of node coordinates)
  //
  // The cell centroid is computed as the volume weighted average of the
  // centroids of tetrahedra from a symmetric tetrahedral
  // decomposition of the cell. The tetrahedral decomposition is
  // formed by connecting the cell center (average of cell nodes), a
  // face center (average of face nodes) and the two nodes of an edge
  // of the face
  void build_cell_centroid() const;

  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
  cell_centroid(const Entity_ID cellid) const
  {
    assert(cell_geometry_precomputed_);
    return cell_centroids_(cellid);
  }

  // Centroid of face (center of gravity not just the average of node
  // coordinates)
  //
  // The face centroid is computed as the area weighted average of the
  // centroids of the triangles from a symmetric triangular
  // decomposition of the face. Each triangular facet is formed by the
  // connecting the face center (average of face nodes) to the two
  // nodes of an edge of the face
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
  face_centroid(const Entity_ID faceid, const bool recompute = false) const
  {
    return face_centroids_(faceid);
  }

  void display_cache();


  // Centroid of edge
  AmanziGeometry::Point edge_centroid(const Entity_ID edgeid) const;

  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
  face_normal(const Entity_ID faceid, const bool recompute = false,
              const Entity_ID cellid = -1, int* orientation = NULL) const
  {
    assert(face_geometry_precomputed_);
    assert(recompute == false);
    assert(faces_requested_);

    Kokkos::View<AmanziGeometry::Point*> fnormals, fnormals_new;
    fnormals =
      Kokkos::subview(face_normals_.entries,
                      Kokkos::make_pair(face_normals_.row_map(faceid),
                                        face_normals_.row_map(faceid + 1)));

    assert(fnormals.extent(0) > 0);

    if (cellid == -1) {
      // Return the natural normal. This is the normal with respect to
      // the first cell, appropriately adjusted according to whether the
      // face is pointing into the cell (-ve cell id) or out

      int c = face_cell_ids_.entries(face_cell_ids_.row_map(faceid));
      return c < 0 ? -fnormals(0) : fnormals(0);
    } else {
      // Find the index of 'cellid' in list of cells connected to face

      int dir;
      int irefcell;
      int nfc =
        face_cell_ids_.row_map(faceid + 1) - face_cell_ids_.row_map(faceid);
      for (irefcell = 0; irefcell < nfc; irefcell++) {
        int c =
          face_cell_ids_.entries(face_cell_ids_.row_map(faceid) + irefcell);
        if (c == cellid || ~c == cellid) {
          dir = c < 0 ? -1 : 1;
          break;
        }
      }
      assert(irefcell < nfc);
      if (orientation) *orientation = dir; // if orientation was requested
      return fnormals(irefcell);
    }
  }


  // Edge vector
  //
  // Vector length equals to edge length.
  //
  // If recompute is TRUE, then the vector is recalculated using current
  // edge coordinates but not stored. (If the recomputed vector must be
  // stored, then call recompute_geometric_quantities).
  //
  // If pointid is not specified, the vector is the natural direction of
  // the edge (from point0 to point1).  On the other hand, if pointid
  // is specified (has to be a point of the face), the vector is from
  // specified point to opposite point of edge.
  //
  // if pointid is specified, then orientation returns the direction of
  // the natural direction of the edge with respect to the point (1 is
  // away from the point and -1 is towards)
  AmanziGeometry::Point
  edge_vector(const Entity_ID edgeid, const bool recompute = false,
              const Entity_ID pointid = -1, int* orientation = NULL) const;

  // Is a point in a given cell?
  bool
  point_in_cell(const AmanziGeometry::Point& p, const Entity_ID cellid) const;
  //
  // Miscellaneous functions
  //------------------------

  // Helper function to build columns
  int build_single_column_(int colnum, Entity_ID top_face) const;

  //
  // Cache management
  //-----------------

  // Get edges of a face and directions in which the face uses the edges.
  //
  // In 3D, edge direction is 1 when it is oriented counter clockwise
  // with respect to the face natural normal.
  //
  // On a distributed mesh, this will return all the edges of the
  // face, OWNED or GHOST. If ordered = true, the edges will be
  // returned in a ccw order around the face as it is naturally defined.
  //
  // IMPORTANT NOTE IN 2D CELLS: In meshes where the cells are two
  // dimensional, faces and edges are identical. For such cells, this
  // operator will return a single edge and a direction of 1. However,
  // this direction cannot be relied upon to compute, say, a contour
  // integral around the 2D cell.
  void
  face_get_edges_and_dirs(const Entity_ID faceid,
                                Kokkos::View<Entity_ID*>& edgeids,
                                Kokkos::View<int*>* edge_dirs,
                                const bool ordered = false) const
  {
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
  }


  // Get the local index of a face edge in a cell edge list
  // Example:
  //
  // face_get_edges(face=5) --> {20, 21, 35, 9, 10}
  // cell_get_edges(cell=18) --> {1, 2, 3, 5, 8, 9, 10, 13, 21, 35, 20, 37, 40}
  // face_to_cell_edge_map(face=5,cell=18) --> {10, 8, 9, 5, 6}
  void
  face_to_cell_edge_map(const Entity_ID faceid, const Entity_ID cellid,
                              std::vector<int>* map) const
  {
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
  }

  // Get edges of a cell
  void cell_get_edges(const Entity_ID cellid,
                      Kokkos::View<Entity_ID*>& edgeids) const
  {
    if (!cell2edge_info_cached_) cache_cell2edge_info_();

    edgeids = Kokkos::subview(cell_edge_ids_.entries,
                              std::make_pair(cell_edge_ids_.row_map(cellid),
                                            cell_edge_ids_.row_map(cellid + 1)));
    // Entity_ID_List &cedgeids = cell_edge_ids_[cellid];
    //*edgeids = cell_edge_ids_[cellid]; // copy operation
  }

  // edges and directions of a 2D cell
  void
  cell_2D_get_edges_and_dirs(const Entity_ID cellid,
                                  Kokkos::View<Entity_ID*>& edgeids,
                                  Kokkos::View<int*>* edgedirs) const
  {
    if (!cell2edge_info_cached_) cache_cell2edge_info_();

    edgeids = Kokkos::subview(cell_edge_ids_.entries,
                              std::make_pair(cell_edge_ids_.row_map(cellid),
                                            cell_edge_ids_.row_map(cellid + 1)));
    *edgedirs =
      Kokkos::subview(cell_2D_edge_dirs_.entries,
                      std::make_pair(cell_2D_edge_dirs_.row_map(cellid),
                                    cell_2D_edge_dirs_.row_map(cellid + 1)));
  }

  Kokkos::View<double*> cell_volumes() const { return cell_volumes_; }

  // Virtual methods to fill the cache with geometric quantities.
  //
  // Default implementations use _internal() methods below.
  int compute_cell_geometric_quantities_() const;
  int compute_face_geometric_quantities_() const;
  int compute_edge_geometric_quantities_() const;

  // Virtual methods to fill the cache with topological quantities.
  //
  // Default implementations use _internal() methods below.
  void cache_cell_face_info_() const;
  void cache_cell2edge_info_() const;
  void cache_face2edge_info_() const;

  // Virtual methods to fill the cache with geometric quantities.
  //
  // Convenience methods that wrap multiple calls.
  //
  // These are declared const since they do not modify the
  // mesh but just modify cached variables declared as mutable
  int compute_cell_geometry_(const Entity_ID cellid, double* volume,
                                     AmanziGeometry::Point* centroid) const;
  int
  compute_face_geometry_(const Entity_ID faceid, double* area,
                         AmanziGeometry::Point* centroid,
                         std::vector<AmanziGeometry::Point>& normals) const;

  int compute_edge_geometry_(const Entity_ID edgeid, double* length,
                                     AmanziGeometry::Point* edge_vector) const;

  void cache_parents_info_() const; 

  void cache_cell_get_faces_and_bisectors_() const; 

  // the cache
  // -- geometry
  mutable Kokkos::View<double*> cell_volumes_, face_areas_, edge_lengths_;
  mutable Kokkos::View<AmanziGeometry::Point*> cell_centroids_, face_centroids_;

  // -- Have to account for the fact that a "face" for a non-manifold
  // surface mesh can have more than one cell connected to
  // it. Therefore, we have to store as many normals for a face as
  // there are cells connected to it. For a given face, its normal to
  // face_get_cells()[i] is face_normals_[i]
  mutable Kokkos::Crs<AmanziGeometry::Point, Kokkos::DefaultExecutionSpace>
    face_normals_;

  mutable Kokkos::View<AmanziGeometry::Point*> edge_vectors_;

  mutable Kokkos::View<Entity_ID*> cells_parent_; 
  mutable Kokkos::View<Entity_ID*> faces_parent_; 
  mutable Kokkos::View<Entity_ID*> edges_parent_; 
  mutable Kokkos::View<Entity_ID*> nodes_parent_; 

  mutable Kokkos::Crs<AmanziGeometry::Point, Kokkos::DefaultExecutionSpace> cell_faces_bisectors_;


  // -- column information, only created if columns are requested
  mutable Kokkos::View<Entity_ID*> cell_cellabove_, cell_cellbelow_,
    node_nodeabove_;
  mutable std::vector<Entity_ID_List> column_cells_; 
  mutable std::vector<Entity_ID_List> column_faces_; 
  //mutable Kokkos::Crs<Entity_ID, Kokkos::DefaultExecutionSpace> column_cells_;
  //mutable Kokkos::Crs<Entity_ID, Kokkos::DefaultExecutionSpace> column_faces_;
  mutable Kokkos::View<Entity_ID*> columnID_;
  mutable int num_owned_cols_;
  mutable bool columns_built_;

  // -- topology
  mutable Kokkos::Crs<Entity_ID, Kokkos::DefaultExecutionSpace> cell_face_ids_;
  mutable Kokkos::Crs<int, Kokkos::DefaultExecutionSpace>
    cell_face_dirs_; // 1 or -1

  // 1s complement if face is pointing out of cell; cannot use 0 as
  // cellid can be 0
  mutable Kokkos::Crs<Entity_ID, Kokkos::DefaultExecutionSpace> face_cell_ids_;
  mutable Kokkos::View<Entity_ID*> face_cell_ids_dirs_;

  mutable Kokkos::Crs<Parallel_type, Kokkos::DefaultExecutionSpace>
    face_cell_ptype_;
  mutable Kokkos::Crs<Entity_ID, Kokkos::DefaultExecutionSpace> cell_edge_ids_;
  mutable Kokkos::Crs<int, Kokkos::DefaultExecutionSpace> cell_2D_edge_dirs_;
  mutable Kokkos::Crs<Entity_ID, Kokkos::DefaultExecutionSpace> face_edge_ids_;
  mutable Kokkos::Crs<int, Kokkos::DefaultExecutionSpace> face_edge_dirs_;

  // -- flags to indicate what part of cache is up-to-date
  mutable bool cell2face_info_cached_, face2cell_info_cached_;
  mutable bool cell2edge_info_cached_, face2edge_info_cached_;
  mutable bool cell_geometry_precomputed_, face_geometry_precomputed_,
    edge_geometry_precomputed_;
  mutable bool parents_precomputed_, cell_get_faces_and_bisectors_precomputed_; 

  // -- region data
  mutable std::map<std::string, std::vector<int>> region_ids;
  mutable std::map<std::string, std::vector<double>> region_vofs;

  // probably should not be mutable?  these should be set by constructor and not
  // changed! --etc
  mutable bool faces_requested_, edges_requested_;

  // friend classes change the cache?  why is this necessary? --etc
  friend class MeshEmbeddedLogical;

  // fast search tools
  mutable bool kdtree_faces_initialized_;
  mutable KDTree kdtree_faces_;

  Mesh* mesh_; 
};

} // namespace AmanziMesh
} // namespace Amanzi

#endif
