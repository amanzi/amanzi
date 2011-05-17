// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: Mesh_STK.hh
// -------------------------------------------------------------
/**
 * @file   Mesh_STK.hh
 * @author William A. Perkins
 * @date Tue May 17 11:42:40 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// Created May  2, 2011 by William A. Perkins
// Last Change: Tue May 17 11:42:40 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _Mesh_STK_hh_
#define _Mesh_STK_hh_

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MpiComm.h>
#include <memory>

#include "Mesh.hh"
#include "Data_structures.hh"
#include "Entity_map.hh"

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------
//  class Mesh_STK
// -------------------------------------------------------------
class Mesh_STK : public Amanzi::AmanziMesh::Mesh {
 public:

  /// Construct a mesh from a
  explicit Mesh_STK(STK::Mesh_STK_Impl_p mesh);

  /// Construct hexahedral mesh of the given size and spacing
  Mesh_STK(const Epetra_MpiComm& comm, 
           const unsigned int& ni, const unsigned int& nj, const unsigned int& nk,
           const double& xorigin = 0.0, 
           const double& yorigin = 0.0, 
           const double& zorigin = 0.0, 
           const double& xdelta = 1.0, 
           const double& ydelta = 1.0, 
           const double& zdelta = 1.0);

  /// Construct hexahedral mesh (Mesh_maps_simple alternative)
  Mesh_STK(double x0, double y0, double z0,
           double x1, double y1, double z1,
           int nx, int ny, int nz, 
           Epetra_MpiComm *communicator);

  /// Construct a mesh from a Exodus II file or file set
  Mesh_STK(const Epetra_MpiComm& comm, 
           const std::string& fname);

  /// Construct a mesh from a Exodus II file or file set
  Mesh_STK(const char *filename, MPI_Comm comm);

  /// Destructor
  ~Mesh_STK(void);

  /// Get parallel type of eneity
  Parallel_type entity_get_ptype(const Entity_kind kind, 
                                 const Entity_ID entid) const;

  // Get cell type
  Cell_type cell_get_type(const Entity_ID cellid) const;

  //
  // General mesh information
  // -------------------------
  //
  
  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, USED)
  unsigned int num_entities (const Entity_kind kind,
                             const Parallel_type ptype) const;
    
  // Global ID of any entity
  unsigned int GID(const Entity_ID lid, const Entity_kind kind) const;
    
  /// Get entity local id (0-based) given a global id (0-based)
  Entity_ID LID(const Entity_ID& gid, const Entity_kind& kind) const;

  //
  // Mesh Entity Adjacencies 
  //-------------------------

  // Downward Adjacencies
  //---------------------

  // Get faces of a cell.
  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. The faces will be returned in a standard
  // order according to Exodus II convention.
  void cell_get_faces (const Entity_ID cellid, 
                       Entity_ID_List *faceids);
    
    
  // Get directions in which a cell uses face
  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise
  void cell_get_face_dirs (const Entity_ID cellid, 
                           std::vector<int> *face_dirs);

  // Get nodes of cell 
  // On a distributed mesh, all nodes (OWNED or GHOST) of the cell 
  // are returned
  // Nodes are returned in a standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES in 3D
  // For a general polyhedron this will return the nodes in
  // arbitrary order
  // In 2D, the nodes of the polygon will be returned in ccw order 
  // consistent with the face normal
  void cell_get_nodes (const Entity_ID cellid, 
                       Entity_ID_List *nodeids);

  // Get nodes of face 
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face 
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  void face_get_nodes (const Entity_ID faceid, 
                       Entity_ID_List *nodeids);

  // Upward adjacencies
  //-------------------
  
  // Cells of type 'ptype' connected to a node
  void node_get_cells (const Entity_ID nodeid, 
                       const Parallel_type ptype,
                       Entity_ID_List *cellids);
    
  // Faces of type 'ptype' connected to a node
  void node_get_faces (const Entity_ID nodeid, 
                       const Parallel_type ptype,
                       Entity_ID_List *faceids);
    
  // Get faces of ptype of a particular cell that are connected to the
  // given node
  void node_get_cell_faces (const Entity_ID nodeid, 
                            const Entity_ID cellid,
                            const Parallel_type ptype,
                            Entity_ID_List *faceids);    
    
  // Cells connected to a face
  void face_get_cells (const Entity_ID faceid, 
                       const Parallel_type ptype,
                       Entity_ID_List *cellids);
    


  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = USED, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces
  void cell_get_face_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List *fadj_cellids);

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order
  void cell_get_node_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List *nadj_cellids);

  //
  // Mesh Topology for viz  
  //----------------------
  //
  // We need a special function because certain types of degenerate
  // hexes will not be recognized as any standard element type (hex,
  // pyramid, prism or tet). The original topology of this element 
  // without any collapsed nodes will be returned by this call.
    
  // Original cell type 
  Cell_type cell_get_type_4viz(const Entity_ID cellid) const;
    
    
  // See cell_get_nodes for details on node ordering
  void cell_get_nodes_4viz (const Entity_ID cellid, 
                            Entity_ID_List *nodeids);
    
    
  //
  // Mesh entity geometry
  //--------------
  //
    
  // Node coordinates - 3 in 3D and 2 in 2D
  void node_get_coordinates (const Entity_ID nodeid, 
			     AmanziGeometry::Point *ncoord);
    
    
  // Face coordinates - conventions same as face_to_nodes call 
  // Number of nodes is the vector size divided by number of spatial dimensions
  void face_get_coordinates (const Entity_ID faceid, 
			     std::vector<AmanziGeometry::Point> *fcoords); 
    
  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions
  void cell_get_coordinates (const Entity_ID cellid, 
			     std::vector<AmanziGeometry::Point> *ccoords);


  //
  // Epetra maps
  //------------
  const Epetra_Map& cell_epetra_map (const bool include_ghost) const;
  const Epetra_Map& face_epetra_map (const bool include_ghost) const; 
  const Epetra_Map& node_epetra_map (const bool include_ghost) const;
    
  //
  // Boundary Conditions or Sets
  //----------------------------
  //
    
  // Number of sets containing entities of type 'kind' in mesh
  unsigned int num_sets(const Entity_kind kind) const;
    
    
  // Ids of sets containing entities of 'kind'
  void get_set_ids (const Entity_kind kind, 
                    Set_ID_List *setids);


  // Is this is a valid ID of a set containing entities of 'kind'
  bool valid_set_id (const Set_ID setid, 
                     const Entity_kind kind) const;


  // Get number of entities of type 'category' in set
  unsigned int get_set_size (const Set_ID setid, 
                             const Entity_kind kind,
                             const Parallel_type ptype);


  // Get list of entities of type 'category' in set
  void get_set_entities (const Set_ID setid, 
                         const Entity_kind kind, 
                         const Parallel_type ptype, 
                         Entity_ID_List *entids); 

  /// Make a cell-to-cell graph from the mesh
  Teuchos::RCP<Epetra_CrsGraph> cellgraph(void) const;

  /// Redistribute the mesh according to the specified cell map
  void redistribute(const Epetra_Map& cellmap);

  /// Repartition and redistribute the mesh according to the specified parameters
  void redistribute(const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"));

 private:

  /// A list of supported entity kinds
  static const unsigned int num_kinds_;
  static const Entity_kind kinds_[];

  /// The parallel environment
  std::auto_ptr<Epetra_MpiComm> communicator_;

  /// The mesh implementation
  STK::Mesh_STK_Impl_p mesh_;

  /// A thing to relate stk::mesh ranks to Entity_kind
  STK::Entity_map entity_map_;

  // Maps, Accessors and setters.
  // ----------------------------
  
  /// A thing to relate Epetra_Maps to mesh entities
  typedef std::map< Entity_kind, Teuchos::RCP<Epetra_Map> > MapSet;

  MapSet map_owned_;          /**< The Epetra_Map's for owned entities */
  MapSet map_used_;           /**< The Epetra_Map's for used (owned + ghost) entities */
    
  /// Build a mesh from a Exodus II file or file set
  void read_exodus_(const std::string& fname);
  
  /// Generate a hexahedral mesh 
  void generate_(const unsigned int& ni, const unsigned int& nj, const unsigned int& nk,
                 const double& xorigin, 
                 const double& yorigin, 
                 const double& zorigin, 
                 const double& xdelta, 
                 const double& ydelta, 
                 const double& zdelta);

  /// Build and store the required Epetra_Map instances
  void build_maps_ ();

  /// Get the appropriate map for the specified @c kind
  const Epetra_Map& get_map_(const Entity_kind& kind, 
                             const bool& include_ghost) const;

};

} // namespace AmanziMesh
} // namespace Amanzi

#endif
