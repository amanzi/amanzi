#pragma once

#include <string>

#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "GeometricModel.hh"

#include "MeshDefs.hh"
#include "MeshSets.hh"
#include "MeshColumns.hh"


namespace Amanzi {
namespace AmanziMesh {

class MeshFramework;
struct MeshAlgorithms;
struct MeshCacheDevice; 


struct MeshCacheData {
  // flags
  bool cell_geometry_cached = false;
  bool cell_faces_cached = false;
  bool cell_edges_cached = false;
  bool cell_nodes_cached = false;
  bool cell_coordinates_cached = false;

  bool face_geometry_cached = false;
  bool face_cells_cached = false;
  bool face_edges_cached = false;
  bool face_nodes_cached = false;
  bool face_coordinates_cached = false;

  bool edge_geometry_cached = false;
  bool edge_cells_cached = false;
  bool edge_faces_cached = false;
  bool edge_nodes_cached = false;
  bool edge_coordinates_cached = false;
  bool edge_lengths_cached = false;

  bool node_cells_cached = false;
  bool node_faces_cached = false;
  bool node_edges_cached = false;
  bool node_coordinates_cached = false;

  bool parent_entities_cached = false;
  bool cell_cellbelow_cached = false;

  bool cell_global_indices_cached = false; 

  // geometry
  Point_DualView node_coordinates;
  Point_DualView cell_centroids;
  Point_DualView face_centroids;
  Point_DualView edge_centroids;
  RaggedArray_DualView<AmanziGeometry::Point> cell_coordinates;
  RaggedArray_DualView<AmanziGeometry::Point> face_coordinates;
  RaggedArray_DualView<AmanziGeometry::Point> edge_coordinates;

  Double_DualView cell_volumes;
  Double_DualView face_areas;
  Double_DualView edge_lengths;

  RaggedArray_DualView<AmanziGeometry::Point> face_normals;
  RaggedArray_DualView<Direction_type> face_normal_orientations;
  Point_DualView edge_vectors;

  // downward adjacencies
  RaggedArray_DualView<Entity_ID> cell_faces;
  RaggedArray_DualView<Direction_type> cell_face_directions;
  RaggedArray_DualView<AmanziGeometry::Point> cell_face_bisectors;

  RaggedArray_DualView<Entity_ID> cell_edges;
  RaggedArray_DualView<Entity_ID> cell_nodes;
  RaggedArray_DualView<Entity_ID> face_edges;
  RaggedArray_DualView<Direction_type> face_edge_directions;
  RaggedArray_DualView<Entity_ID> face_nodes;
  RaggedArray_DualView<Entity_ID> edge_nodes;
  Entity_ID_DualView cell_cellbelow;

  // upward adjacencies
  RaggedArray_DualView<Entity_ID> face_cells;
  RaggedArray_DualView<Entity_ID> edge_cells;
  RaggedArray_DualView<Entity_ID> edge_faces;
  RaggedArray_DualView<Entity_ID> node_cells;
  RaggedArray_DualView<Entity_ID> node_faces;
  RaggedArray_DualView<Entity_ID> node_edges;

  // parent entities
  Entity_ID_DualView parent_nodes;
  Entity_ID_DualView parent_edges;
  Entity_ID_DualView parent_faces;
  Entity_ID_DualView parent_cells;

  // Map
  Entity_GID_DualView cell_global_indices;

  MeshCacheData(): node_coordinates("node_coordinates", 0),
                   cell_centroids("cell_centroids", 0),
                   face_centroids("face_centroids", 0), 
                   edge_centroids("edge_centroids", 0),
                   cell_volumes("cell_volumes", 0),
                   face_areas("face_areas", 0),
                   edge_lengths("edge_lengths", 0),
                   cell_cellbelow("cell_cellbelow", 0){}


};


struct MeshCacheBase {
 protected:
  MeshCacheBase()
    : ncells_owned(-1),
      ncells_all(-1),
      nfaces_owned(-1),
      nfaces_all(-1),
      nedges_owned(-1),
      nedges_all(-1),
      nnodes_owned(-1),
      nnodes_all(-1),
      nboundary_faces_owned(-1),
      nboundary_faces_all(-1),
      nboundary_nodes_owned(-1),
      nboundary_nodes_all(-1),
      space_dim_(-1),
      manifold_dim_(-1),
      is_ordered_(false),
      is_logical_(false),
      has_edges_(false),
      has_nodes_(true),
      has_node_faces_(true),
      is_sfm_(false)
  {}


  MeshCacheBase(const Teuchos::RCP<Teuchos::ParameterList>& plist) : MeshCacheBase()
  {
    if (plist == Teuchos::null) {
      plist_ = Teuchos::rcp(new Teuchos::ParameterList());
    } else {
      plist_ = plist;
    }
  }

  MeshCacheBase(const Teuchos::RCP<MeshFramework>& framework_mesh,
                const Teuchos::RCP<MeshAlgorithms>& framework_algorithms,
                const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : MeshCacheBase(plist)
  {
    if (framework_algorithms == Teuchos::null) {
      Errors::Message msg("MeshAlgorithms is Teuchos::null.");
      Exceptions::amanzi_throw(msg);
    }
    algorithms_ = framework_algorithms;
  }

  
  MeshCacheBase(const MeshCacheBase& other) = default;

  ~MeshCacheBase() = default; 
  
 public:
  // sizes
  Entity_ID ncells_owned, ncells_all;
  Entity_ID nfaces_owned, nfaces_all;
  Entity_ID nedges_owned, nedges_all;
  Entity_ID nnodes_owned, nnodes_all;
  Entity_ID nboundary_faces_owned, nboundary_faces_all;
  Entity_ID nboundary_nodes_owned, nboundary_nodes_all;

  //
  // Baseline mesh functionality
  // =============================================
  // mesh properties
  KOKKOS_INLINE_FUNCTION bool isOrdered() const { return is_ordered_; }
  KOKKOS_INLINE_FUNCTION bool isLogical() const { return is_logical_; }
  KOKKOS_INLINE_FUNCTION bool isSFM() const { return is_sfm_; } // single face mesh -- special case

  KOKKOS_INLINE_FUNCTION bool hasNodes() const { return has_nodes_; }
  KOKKOS_INLINE_FUNCTION bool hasEdges() const { return has_edges_; }
  KOKKOS_INLINE_FUNCTION bool hasNodeFaces() const { return has_node_faces_; }

  void hasEdgesOrThrow() const
  {
    if (!hasEdges()) {
      Errors::Message msg("MeshFramework does not include edges.");
      Exceptions::amanzi_throw(msg);
    }
  }

  // ----------------------
  // Accessors and Mutators
  // ----------------------
  Teuchos::RCP<const MeshFramework> getMeshFramework() const { return framework_mesh_; }
  Teuchos::RCP<MeshFramework> getMeshFramework() { return framework_mesh_; }
  void destroyFramework() { framework_mesh_ = Teuchos::null; }
  Teuchos::RCP<Teuchos::ParameterList> getParameterList() const { return plist_; }

  // algorithms class, a partner to the framework, that defines how to compute things
  Teuchos::RCP<const MeshAlgorithms> getAlgorithms() const { return algorithms_; }

  Comm_ptr_type getComm() const { return comm_; }
  void setComm(const Comm_ptr_type& comm) { comm_ = comm; }

  Teuchos::RCP<const AmanziGeometry::GeometricModel> getGeometricModel() const { return gm_; }
  void setGeometricModel(const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm) { gm_ = gm; }

  // space dimension describes the dimension of coordinates in space
  
  KOKKOS_INLINE_FUNCTION std::size_t getSpaceDimension() const { return space_dim_; }
  KOKKOS_INLINE_FUNCTION void setSpaceDimension(unsigned int dim) { space_dim_ = dim; }

  // manifold dimension describes the dimensionality of the corresponding R^n
  // manifold onto which this mesh can be projected.
  
  KOKKOS_INLINE_FUNCTION std::size_t getManifoldDimension() const { return manifold_dim_; }
  void setManifoldDimension(const unsigned int dim) { manifold_dim_ = dim; }

  // column structure
  MeshColumns columns;

 protected:
  // standard things
  Comm_ptr_type comm_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_;
  int space_dim_;
  int manifold_dim_;
  bool is_ordered_;
  bool is_logical_;
  bool has_edges_, has_nodes_;
  bool has_node_faces_;
  bool is_sfm_;

  // related meshes
  Teuchos::RCP<MeshFramework> framework_mesh_;
  Teuchos::RCP<const MeshAlgorithms> algorithms_;

  // helper classes
  MeshCacheData data_;
  MeshMaps maps_;
  mutable MeshSets sets_;
  mutable MeshSetVolumeFractions set_vol_fracs_;
};


// -----------------------------------------------------------------------------
// Memory space transfer
// -----------------------------------------------------------------------------
const MeshCacheDevice
onMemDevice(const MeshCacheHost& mc_in);
  
Teuchos::RCP<const MeshCacheDevice>
onMemDevice(const Teuchos::RCP<const MeshCacheHost>& mc_in);

Teuchos::RCP<MeshCacheDevice>
onMemDevice(const Teuchos::RCP<MeshCacheHost>& mc_in);

const MeshCacheHost
onMemHost(const MeshCacheDevice& mc_in);
  
Teuchos::RCP<const MeshCacheHost>
onMemHost(const Teuchos::RCP<const MeshCacheDevice>& mc_in);

Teuchos::RCP<MeshCacheHost>
onMemHost(const Teuchos::RCP<MeshCacheDevice>& mc_in);
  
}
}
