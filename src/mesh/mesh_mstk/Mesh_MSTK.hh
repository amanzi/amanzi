#ifndef _AMANZI_MESH_MSTK_H_
#define _AMANZI_MESH_MSTK_H_

#include <memory>
#include <vector>
#include <sstream>
#include <Epetra_MpiComm.h>

#include <MSTK.h>

#include <Mesh.hh>
#include <Point.hh>
#include <GeometricModel.hh>
#include <LabeledSetRegion.hh>
#include <PointRegion.hh>
#include <GenerationSpec.hh>
#include <dbc.hh>
#include <errors.hh>

namespace Amanzi 
{
namespace AmanziMesh 
{

class Mesh_MSTK : public Mesh
{
      
private:

  MPI_Comm mpicomm;
  int myprocid, numprocs;

  Mesh_ptr mesh;
  //    MSTKComm *mstkcomm; Not defined - will define if we need it

  int serial_run;
  

  // Local handles to entity lists (Vertices, "Faces", "Cells")
  
  // For a surface mesh, "Faces" refers to mesh edges and "Cells"
  // refers to mesh faces
  //
  // For a solid mesh, "Faces" refers to mesh faces and "Cells"
  // refers to mesh regions
  
  
  // These are MSTK's definitions of types of parallel mesh entities
  // These definitions are slightly different from what Amanzi has defined
  //
  // There are 2 types of entities relevant to this code - Owned and Ghost
  //
  // 1. OWNED - owned by this processor
  //
  // 2. GHOST - not owned by this processor
  
  // ALL = OWNED + GHOST
  

  MSet_ptr AllVerts, OwnedVerts, NotOwnedVerts;
  
  MSet_ptr AllFaces, OwnedFaces, NotOwnedFaces;
  
  MSet_ptr AllCells, OwnedCells, GhostCells;


  // Attribute handle to store cell_type and original cell_type

  MAttrib_ptr celltype_att, orig_celltype_att, orig_celltopo_att;
      

  // Local ID to MSTK handle map

  std::vector<MEntity_ptr> vtx_id_to_handle;
  std::vector<MEntity_ptr> face_id_to_handle;
  std::vector<MEntity_ptr> cell_id_to_handle;


  // Maps

  Epetra_Map *cell_map_wo_ghosts_, *face_map_wo_ghosts_, *node_map_wo_ghosts_;
  Epetra_Map *cell_map_w_ghosts_, *face_map_w_ghosts_, *node_map_w_ghosts_;
  
  
  // flag whether to flip a face dir or not when returning nodes of a face
  
  bool *faceflip;
  
  // Private methods
  // ----------------------------
  
  void clear_internals_();

  void pre_create_steps_(const int space_dimension, const Epetra_MpiComm *incomm, 
                         const AmanziGeometry::GeometricModelPtr& gm);
  void post_create_steps_();

  void collapse_degen_edges();
  void label_celltype();

  void init_pvert_lists();
  void init_pface_lists();
  void init_pcell_lists();
  void init_pface_dirs();
  
  void init_id_handle_maps();
  void init_global_ids();
  
  void init_cell_map();
  void init_face_map();
  void init_node_map();
  
  void init_set_info();

  int  generate_regular_mesh(Mesh_ptr mesh, double x0, double y0, double z0,
			     double x1, double y1, double z1, int nx,
			     int ny, int nz);
  int  generate_regular_mesh(Mesh_ptr mesh, double x0, double y0,
			     double x1, double y1, int nx, int ny);
public:

  // Constructors that read the mesh from a file

  Mesh_MSTK (const char *filename, const Epetra_MpiComm *incomm,
	     const AmanziGeometry::GeometricModelPtr& gm = 
	     (AmanziGeometry::GeometricModelPtr) NULL);

  Mesh_MSTK (const char *filename, const Epetra_MpiComm *incomm, int space_dimension,
	     const AmanziGeometry::GeometricModelPtr& gm = 
	     (AmanziGeometry::GeometricModelPtr) NULL);

  // Constructors that generate a mesh internally (regular hexahedral mesh only)

  // 3D
  Mesh_MSTK(const double x0, const double y0, const double z0,
	    const double x1, const double y1, const double z1,
	    const unsigned int nx, const unsigned int ny, const unsigned int nz, 
           const Epetra_MpiComm *incomm,
           const AmanziGeometry::GeometricModelPtr& gm = 
           (AmanziGeometry::GeometricModelPtr) NULL);


  // 2D
  Mesh_MSTK(const double x0, const double y0,
	    const double x1, const double y1,
	    const int nx, const int ny, 
	    const Epetra_MpiComm *comm,
	    const AmanziGeometry::GeometricModelPtr& gm = 
	    (AmanziGeometry::GeometricModelPtr) NULL);

  // Construct a hexahedral mesh from specs 
  Mesh_MSTK(const GenerationSpec& gspec,
	    const Epetra_MpiComm *comm,
	    const AmanziGeometry::GeometricModelPtr& gm = 
	    (AmanziGeometry::GeometricModelPtr) NULL);


  ~Mesh_MSTK ();


  // Get parallel type of entity
    
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
		       Entity_ID_List *faceids) const;
    
    
  // Get directions in which a cell uses face
  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise
    
  void cell_get_face_dirs (const Entity_ID cellid, 
			   std::vector<int> *face_dirs) const;
    
    
    
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
		       Entity_ID_List *nodeids) const;
    
    
  // Get nodes of face 
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face 
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
    
  void face_get_nodes (const Entity_ID faceid, 
		       Entity_ID_List *nodeids) const;
    


  // Upward adjacencies
  //-------------------
    
  // Cells of type 'ptype' connected to a node
    
  void node_get_cells (const Entity_ID nodeid, 
		       const Parallel_type ptype,
		       Entity_ID_List *cellids) const;
    
  // Faces of type 'ptype' connected to a node
    
  void node_get_faces (const Entity_ID nodeid, 
		       const Parallel_type ptype,
		       Entity_ID_List *faceids) const;
    
  // Get faces of ptype of a particular cell that are connected to the
  // given node
    
  void node_get_cell_faces (const Entity_ID nodeid, 
			    const Entity_ID cellid,
			    const Parallel_type ptype,
			    Entity_ID_List *faceids) const;    
    
  // Cells connected to a face
    
  void face_get_cells (const Entity_ID faceid, 
		       const Parallel_type ptype,
		       Entity_ID_List *cellids) const;
    


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
			       Entity_ID_List *fadj_cellids) const;

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order

  void cell_get_node_adj_cells(const Entity_ID cellid,
			       const Parallel_type ptype,
			       Entity_ID_List *nadj_cellids) const;


    
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
			    Entity_ID_List *nodeids) const;
    
    
    
  //
  // Mesh entity geometry
  //--------------
  //
    
  // Node coordinates - 3 in 3D and 2 in 2D
    
  void node_get_coordinates (const Entity_ID nodeid, 
			     AmanziGeometry::Point *ncoord) const;
    
    
  // Face coordinates - conventions same as face_to_nodes call 
  // Number of nodes is the vector size divided by number of spatial dimensions
    
  void face_get_coordinates (const Entity_ID faceid, 
			     std::vector<AmanziGeometry::Point> *fcoords) const; 
    
  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions
    
  void cell_get_coordinates (const Entity_ID cellid, 
			     std::vector<AmanziGeometry::Point> *ccoords) const;
    
  // Modify the coordinates of a node
  void node_set_coordinates (const Entity_ID nodeid, const double *coords);

    
  //
  // Epetra maps
  //------------
    
    
  const Epetra_Map& cell_epetra_map (bool include_ghost) const;
    
  const Epetra_Map& face_epetra_map (bool include_ghost) const; 

  const Epetra_Map& node_epetra_map (bool include_ghost) const;
    
    
    
    
  //
  // Boundary Conditions or Sets
  //----------------------------
  //
    

  // Get number of entities of type 'category' in set

  unsigned int get_set_size (const Set_ID setid, 
			     const Entity_kind kind,
			     const Parallel_type ptype) const;

  unsigned int get_set_size (const Set_Name setname, 
			     const Entity_kind kind,
			     const Parallel_type ptype) const;


  unsigned int get_set_size (const char *setname, 
			     const Entity_kind kind,
			     const Parallel_type ptype) const;

  // Get list of entities of type 'category' in set

  void get_set_entities (const Set_ID setid, 
			 const Entity_kind kind, 
			 const Parallel_type ptype, 
			 Entity_ID_List *entids) const; 

  void get_set_entities (const Set_Name setname, 
			 const Entity_kind kind, 
			 const Parallel_type ptype, 
			 std::vector<Entity_ID> *entids) const; 


  void get_set_entities (const char *setname, 
			 const Entity_kind kind, 
			 const Parallel_type ptype, 
			 std::vector<Entity_ID> *entids) const; 


};



    
  inline Parallel_type Mesh_MSTK::entity_get_ptype(const Entity_kind kind, const Entity_ID entid) const {
    MEntity_ptr ment;
      
    switch(kind) {
    case CELL:
      ment = (MEntity_ptr) cell_id_to_handle[entid];
      break;
    case FACE:
      ment = (MEntity_ptr) face_id_to_handle[entid];
      break;
    case NODE:
      ment = (MEntity_ptr) vtx_id_to_handle[entid];
      break;
    }
      
    if (MEnt_PType(ment) == PGHOST)
      return GHOST;
    else
      return OWNED;
  }

} // End namespace AmanziMesh

} // End namespace Amanzi

#endif /* _AMANZI_MESH_MSTK_H_ */
