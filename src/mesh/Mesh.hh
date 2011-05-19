#ifndef _AMANZI_MESH_H_
#define _AMANZI_MESH_H_

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

#include <memory>

#include "MeshDefs.hh"
#include "Cell_topology.hh"
#include "Point.hh"

namespace Amanzi
{

namespace AmanziMesh
{

  class Mesh
  {

  private:

    unsigned int celldim, spacedim;
    bool geometry_precomputed;
    int precompute_geometric_quantities();
    std::vector<double> cell_volumes, face_areas;
    std::vector<AmanziGeometry::Point> cell_centroids, face_normals;

    Epetra_Comm *comm; // temporary until we get an amanzi communicator
    
  public:
      
    // constructor 

    Mesh()
      : spacedim(3), celldim(3), geometry_precomputed(false), comm(NULL)
    {
    }

    // destructor

    virtual ~Mesh() {};

    inline
    void set_space_dimension(const unsigned int dim) {
      spacedim = dim;
    }

    inline
    unsigned int space_dimension() const 
    {
      return spacedim;
    }

    inline
    void set_cell_dimension(const unsigned int dim) {
      celldim = dim;   // 3 is solid mesh, 2 is surface mesh
    }

    inline 
    unsigned int cell_dimension() const
    {
      return celldim;
    }



    // Get parallel type of eneity
    
    virtual
    Parallel_type entity_get_ptype(const Entity_kind kind, 
				   const Entity_ID entid) const = 0;




    // Get cell type
    
    virtual
    Cell_type cell_get_type(const Entity_ID cellid) const = 0;
        
    


    //
    // General mesh information
    // -------------------------
    //
    
    // Number of entities of any kind (cell, face, node) and in a
    // particular category (OWNED, GHOST, USED)
    
    virtual
    unsigned int num_entities (const Entity_kind kind,
			       const Parallel_type ptype) const = 0;
    
    
    // Global ID of any entity
    
    virtual
    unsigned int GID(const Entity_ID lid, const Entity_kind kind) const = 0;
    
    
    
    //
    // Mesh Entity Adjacencies 
    //-------------------------


    // Downward Adjacencies
    //---------------------
    
    // Get faces of a cell.

    // On a distributed mesh, this will return all the faces of the
    // cell, OWNED or GHOST. The faces will be returned in a standard
    // order according to Exodus II convention.
    
    virtual
    void cell_get_faces (const Entity_ID cellid, 
			 Entity_ID_List *faceids) const = 0;
    
    
    // Get directions in which a cell uses face
    // In 3D, direction is 1 if face normal points out of cell
    // and -1 if face normal points into cell
    // In 2D, direction is 1 if face/edge is defined in the same
    // direction as the cell polygon, and -1 otherwise
    
    virtual
    void cell_get_face_dirs (const Entity_ID cellid, 
			     std::vector<int> *face_dirs) const = 0;
    
    
    
    // Get nodes of cell 
    // On a distributed mesh, all nodes (OWNED or GHOST) of the cell 
    // are returned
    // Nodes are returned in a standard order (Exodus II convention)
    // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES in 3D
    // For a general polyhedron this will return the nodes in
    // arbitrary order
    // In 2D, the nodes of the polygon will be returned in ccw order 
    // consistent with the face normal
    
    virtual
    void cell_get_nodes (const Entity_ID cellid, 
			 Entity_ID_List *nodeids) const = 0;
    
    
    // Get nodes of face 
    // On a distributed mesh, all nodes (OWNED or GHOST) of the face 
    // are returned
    // In 3D, the nodes of the face are returned in ccw order consistent
    // with the face normal
    // In 2D, nfnodes is 2
    
    virtual
    void face_get_nodes (const Entity_ID faceid, 
			 Entity_ID_List *nodeids) const = 0;
    


    // Upward adjacencies
    //-------------------
    
    // Cells of type 'ptype' connected to a node
    
    virtual 
    void node_get_cells (const Entity_ID nodeid, 
			 const Parallel_type ptype,
			 Entity_ID_List *cellids) const = 0;
    
    // Faces of type 'ptype' connected to a node
    
    virtual
    void node_get_faces (const Entity_ID nodeid, 
			 const Parallel_type ptype,
			 Entity_ID_List *faceids) const = 0;
    
    // Get faces of ptype of a particular cell that are connected to the
    // given node
    
    virtual
    void node_get_cell_faces (const Entity_ID nodeid, 
			      const Entity_ID cellid,
			      const Parallel_type ptype,
			      Entity_ID_List *faceids) const = 0;    
    
    // Cells connected to a face
    
    virtual 
    void face_get_cells (const Entity_ID faceid, 
			 const Parallel_type ptype,
			 Entity_ID_List *cellids) const = 0;
    


    // Same level adjacencies
    //-----------------------

    // Face connected neighboring cells of given cell of a particular ptype
    // (e.g. a hex has 6 face neighbors)

    // The order in which the cellids are returned cannot be
    // guaranteed in general except when ptype = USED, in which case
    // the cellids will correcpond to cells across the respective
    // faces given by cell_get_faces

    virtual
    void cell_get_face_adj_cells(const Entity_ID cellid,
				 const Parallel_type ptype,
				 Entity_ID_List *fadj_cellids) const = 0;

    // Node connected neighboring cells of given cell
    // (a hex in a structured mesh has 26 node connected neighbors)
    // The cells are returned in no particular order

    virtual
    void cell_get_node_adj_cells(const Entity_ID cellid,
				 const Parallel_type ptype,
				 Entity_ID_List *nadj_cellids) const = 0;


    
    //
    // Mesh Topology for viz  
    //----------------------
    //
    // We need a special function because certain types of degenerate
    // hexes will not be recognized as any standard element type (hex,
    // pyramid, prism or tet). The original topology of this element 
    // without any collapsed nodes will be returned by this call.
    
    
    // Original cell type 
    
    virtual
    Cell_type cell_get_type_4viz(const Entity_ID cellid) const = 0;
    
    
    // See cell_get_nodes for details on node ordering
    
    virtual
    void cell_get_nodes_4viz (const Entity_ID cellid, 
			      Entity_ID_List *nodeids) const = 0;
    
    
    
    //
    // Mesh entity geometry
    //--------------
    //
    
    // Node coordinates - 3 in 3D and 2 in 2D
    
    virtual
    void node_get_coordinates (const Entity_ID nodeid, 
			       AmanziGeometry::Point *ncoord) const = 0;
    
    
    // Face coordinates - conventions same as face_to_nodes call 
    // Number of nodes is the vector size divided by number of spatial dimensions
    
    virtual
    void face_get_coordinates (const Entity_ID faceid, 
			       std::vector<AmanziGeometry::Point> *fcoords) const = 0; 
    
    // Coordinates of cells in standard order (Exodus II convention)
    // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
    // For a general polyhedron this will return the node coordinates in
    // arbitrary order
    // Number of nodes is vector size divided by number of spatial dimensions
    
    virtual 
    void cell_get_coordinates (const Entity_ID cellid, 
			       std::vector<AmanziGeometry::Point> *ccoords) const = 0;
    
    

    // Mesh entity geometry
    //--------------
    //
    
    
    // Volume/Area of cell
    
    inline
    double cell_volume (const Entity_ID cellid) const {
      assert (geometry_precomputed == true);
      return cell_volumes[cellid];
    }
    
    // Area/length of face

    inline
    double face_area(const Entity_ID faceid) const {
      assert (geometry_precomputed == true);
      return face_areas[faceid];
    }
    
    
    // Centroid of cell

    inline
    AmanziGeometry::Point cell_centroid (const Entity_ID cellid) const {
      assert (geometry_precomputed == true);
      return cell_centroids[cellid];
    }
    
    // Normal to face
    // The vector is normalized and then weighted by the area of the face
    
    inline
    AmanziGeometry::Point face_normal (const Entity_ID faceid) const {
      assert (geometry_precomputed == true);
      return face_normals[faceid];
    }


    //
    // Epetra maps
    //------------
    
    
    virtual
    const Epetra_Map& cell_epetra_map (const bool include_ghost) const = 0;
    
    virtual
    const Epetra_Map& face_epetra_map (const bool include_ghost) const = 0; 
    
    virtual
    const Epetra_Map& node_epetra_map (const bool include_ghost) const = 0;
    
    
    
    
    //
    // Boundary Conditions or Sets
    //----------------------------
    //
    
    // Number of sets containing entities of type 'kind' in mesh
    
    virtual
    unsigned int num_sets(const Entity_kind kind) const = 0;
    
    
    // Ids of sets containing entities of 'kind'

    virtual
    void get_set_ids (const Entity_kind kind, 
		      Set_ID_List *setids) const = 0;


    // Is this is a valid ID of a set containing entities of 'kind'

    virtual
    bool valid_set_id (const Set_ID setid, 
		       const Entity_kind kind) const = 0;


    // Get number of entities of type 'category' in set

    virtual
    unsigned int get_set_size (const Set_ID setid, 
			       const Entity_kind kind,
			       const Parallel_type ptype) const = 0;


    // Get list of entities of type 'category' in set

    virtual
    void get_set_entities (const Set_ID setid, 
			   const Entity_kind kind, 
			   const Parallel_type ptype, 
			   Entity_ID_List *entids) const = 0; 






    // communicator access
    // temporary until we set up an amanzi communicator

    inline
    const Epetra_Comm* get_comm() { return comm;};

    inline
    void set_comm(MPI_Comm incomm) {comm = new Epetra_MpiComm(incomm);};



    // Temporary routines for backward compatibility

    void cell_to_faces (unsigned int cell, 
			std::vector<unsigned int>::iterator begin, 
			std::vector<unsigned int>::iterator end) const;

    void cell_to_faces (unsigned int cell, 
			unsigned int* begin, unsigned int *end) const;


    void cell_to_face_dirs (unsigned int cell, 
			    std::vector<int>::iterator begin, 
			    std::vector<int>::iterator end) const;
    void cell_to_face_dirs (unsigned int cell, 
			    int * begin, int * end) const;
  
    

    void cell_to_nodes (unsigned int cell, 
			std::vector<unsigned int>::iterator begin, 
			std::vector<unsigned int>::iterator end) const;
    void cell_to_nodes (unsigned int cell, 
			unsigned int * begin, unsigned int * end) const;
      



    void face_to_nodes (unsigned int face, 
			std::vector<unsigned int>::iterator begin, 
			std::vector<unsigned int>::iterator end) const;
    void face_to_nodes (unsigned int face, 
			unsigned int * begin, unsigned int * end) const;



    void node_to_coordinates (unsigned int node, 
			      std::vector<double>::iterator begin, 
			      std::vector<double>::iterator end) const;
    void node_to_coordinates (unsigned int node, 
			      double * begin, 
			      double * end) const;

    void face_to_coordinates (unsigned int face, 
			      std::vector<double>::iterator begin, 
			      std::vector<double>::iterator end) const;
    void face_to_coordinates (unsigned int face, 
			      double * begin, 
			      double * end) const;
      
    void cell_to_coordinates (unsigned int cell, 
			      std::vector<double>::iterator begin,
			      std::vector<double>::iterator end) const;
    void cell_to_coordinates (unsigned int cell, 
			      double * begin,
			      double * end) const;
      

    const Epetra_Map& cell_map (bool include_ghost) const 
    {
      return cell_epetra_map (include_ghost);
    };
      
    const Epetra_Map& face_map (bool include_ghost) const 
    {
      return face_epetra_map (include_ghost);
    }; 
      
    const Epetra_Map& node_map (bool include_ghost) const 
    {
      return node_epetra_map (include_ghost);
    };

  
    unsigned int count_entities (Entity_kind kind,
				 Parallel_type ptype) const;

    // Unchanged in new interface
    // unsigned int num_sets(Entity_kind kind) const {};
      
    // Unchanged in new interface
    // unsigned int get_set_size (unsigned int set_id, 
    //			 Entity_kind kind,
    //			 Parallel_type ptype) const {};

    // Id numbers
    void get_set_ids (Entity_kind kind, 
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end) const;
    void get_set_ids (Entity_kind kind, 
		      unsigned int * begin, 
		      unsigned int * end) const;

    // Unchanged in new interface
    // bool valid_set_id (unsigned int id, Entity_kind kind) const {};
      
    void get_set (unsigned int set_id, Entity_kind kind, 
		  Parallel_type ptype, 
		  std::vector<unsigned int>::iterator begin, 
		  std::vector<unsigned int>::iterator end) const;
    void get_set (unsigned int set_id, Entity_kind kind, 
		  Parallel_type ptype, 
		  unsigned int * begin, 
		  unsigned int * end) const;

  }; // End class Mesh


} // end namespace AmanziMesh

} // end namespace Amanzi


 

#endif /* _AMANZI_MESH_H_ */
