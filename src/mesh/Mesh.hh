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
    
    public:
      
      // constructor 

      Mesh()
        : spacedim(3), celldim(3), geometry_precomputed(false), comm(NULL)
      {
      }

      // destructor

      ~Mesh() {};

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
			   Entity_ID_List *faceids) = 0;
    
    
      // Get directions in which a cell uses face
      // In 3D, direction is 1 if face normal points out of cell
      // and -1 if face normal points into cell
      // In 2D, direction is 1 if face/edge is defined in the same
      // direction as the cell polygon, and -1 otherwise
    
      virtual
      void cell_get_face_dirs (const Entity_ID cellid, 
			       std::vector<int> *face_dirs) = 0;
    
    
    
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
			   Entity_ID_List *nodeids) = 0;
    
    
      // Get nodes of face 
      // On a distributed mesh, all nodes (OWNED or GHOST) of the face 
      // are returned
      // In 3D, the nodes of the face are returned in ccw order consistent
      // with the face normal
      // In 2D, nfnodes is 2
    
      virtual
      void face_get_nodes (const Entity_ID faceid, 
			   Entity_ID_List *nodeids) = 0;
    


      // Upward adjacencies
      //-------------------
    
      // Cells of type 'ptype' connected to a node
    
      virtual 
      void node_get_cells (const Entity_ID nodeid, 
			   const Parallel_type ptype,
			   Entity_ID_List *cellids) = 0;
    
      // Faces of type 'ptype' connected to a node
    
      virtual
      void node_get_faces (const Entity_ID nodeid, 
			   const Parallel_type ptype,
			   Entity_ID_List *faceids) = 0;
    
      // Get faces of ptype of a particular cell that are connected to the
      // given node
    
      virtual
      void node_get_cell_faces (const Entity_ID nodeid, 
				const Entity_ID cellid,
				const Parallel_type ptype,
				Entity_ID_List *faceids) = 0;    
    
      // Cells connected to a face
    
      virtual 
      void face_get_cells (const Entity_ID faceid, 
			   const Parallel_type ptype,
			   Entity_ID_List *cellids) = 0;
    


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
				   Entity_ID_List *fadj_cellids) = 0;

      // Node connected neighboring cells of given cell
      // (a hex in a structured mesh has 26 node connected neighbors)
      // The cells are returned in no particular order

      virtual
      void cell_get_node_adj_cells(const Entity_ID cellid,
				   const Parallel_type ptype,
				   Entity_ID_List *nadj_cellids) = 0;


    
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
				Entity_ID_List *nodeids) = 0;
    
    
    
      //
      // Mesh entity geometry
      //--------------
      //
    
      // Node coordinates - 3 in 3D and 2 in 2D
    
      virtual
      void node_get_coordinates (const Entity_ID nodeid, 
			     AmanziGeometry::Point *ncoord) = 0;
    
    
      // Face coordinates - conventions same as face_to_nodes call 
      // Number of nodes is the vector size divided by number of spatial dimensions
    
      virtual
      void face_get_coordinates (const Entity_ID faceid, 
			     std::vector<AmanziGeometry::Point> *fcoords) = 0; 
    
      // Coordinates of cells in standard order (Exodus II convention)
      // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
      // For a general polyhedron this will return the node coordinates in
      // arbitrary order
      // Number of nodes is vector size divided by number of spatial dimensions
    
      virtual 
      void cell_get_coordinates (const Entity_ID cellid, 
			     std::vector<AmanziGeometry::Point> *ccoords) = 0;
    
    
      // Volume/Area of cell

      double cell_volume (const Entity_ID cellid);
    
      // Area/length of face

      double face_area(const Entity_ID faceid);
    
      // Centroid of cell

      AmanziGeometry::Point cell_centroid (const Entity_ID cell);
    
      // Normal to face
      // The vector is not normalized or in other words, this is an area
      // weighted normal

      AmanziGeometry::Point face_normal (const Entity_ID face);
    
    
    
      //
      // Epetra maps
      //------------
    
    
      virtual
      inline const Epetra_Map& cell_epetra_map (const bool include_ghost) const = 0;
    
      virtual
      inline const Epetra_Map& face_epetra_map (const bool include_ghost) const = 0; 
    
      virtual
      inline const Epetra_Map& node_epetra_map (const bool include_ghost) const = 0;
    
    
    
    
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
			Set_ID_List *setids) = 0;


      // Is this is a valid ID of a set containing entities of 'kind'

      virtual
      bool valid_set_id (const Set_ID setid, 
			 const Entity_kind kind) const = 0;


      // Get number of entities of type 'category' in set

      virtual
      unsigned int get_set_size (const Set_ID setid, 
				 const Entity_kind kind,
				 const Parallel_type ptype) = 0;


      // Get list of entities of type 'category' in set

      virtual
      void get_set_entities (const Set_ID setid, 
			     const Entity_kind kind, 
			     const Parallel_type ptype, 
			     Entity_ID_List *entids) = 0; 






      // communicator access
      // temporary until we set up an amanzi communicator

      inline
      const Epetra_Comm* get_comm() { return comm;};

      inline
      void set_comm(MPI_Comm incomm) {comm = new Epetra_MpiComm(incomm);};

      // this should be used with extreme caution: WHO USES THIS????
      // IF IT IS ONLY SIMPLE MESH, THEN SIMPLE MESH SHOULD IMPLEMENT IT
      // modify coordinates  
      // virtual void set_coordinate(Entity_ID local_node_id, double *coords) {};

    private:

      unsigned int celldim, spacedim;
      bool geometry_precomputed;
      int precompute_geometric_quantities();
      std::vector<double> cell_volumes, face_areas;
      std::vector<AmanziGeometry::Point> cell_centroids, face_normals;

      Epetra_Comm *comm; // temporary until we get an amanzi communicator

    }; // End class Mesh


  } // end namespace AmanziMesh

} // end namespace Amanzi

#endif /* _AMANZI_MESH_H_ */
