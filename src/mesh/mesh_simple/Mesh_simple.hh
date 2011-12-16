/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef _MESH_SIMPLE_H_
#define _MESH_SIMPLE_H_

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"
#include "Region.hh"

#include "GeometricModel.hh"

namespace Amanzi {

namespace AmanziMesh {

class GenerationSpec;

class Mesh_simple : public virtual Mesh
{
      
public:
      
  Mesh_simple (double x0, double y0, double z0,
	       double x1, double y1, double z1,
	       int nx, int ny, int nz, Epetra_Comm *communicator,
	       const AmanziGeometry::GeometricModelPtr &gm = (AmanziGeometry::GeometricModelPtr) NULL);
  
  Mesh_simple ( const GenerationSpec& gspec,
		Epetra_Comm *communicator,
                const AmanziGeometry::GeometricModelPtr &gm = (AmanziGeometry::GeometricModelPtr) NULL);

  Mesh_simple ( Teuchos::ParameterList &parameter_list,
		Epetra_Comm *communicator,
		const AmanziGeometry::GeometricModelPtr &gm = (AmanziGeometry::GeometricModelPtr) NULL);
  
  virtual ~Mesh_simple ();
  
  void update ();


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
		       std::vector<Entity_ID> *faceids) const;
    
    
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
		       std::vector<Entity_ID> *nodeids) const;
    
    
  // Get nodes of face 
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face 
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
    
  void face_get_nodes (const Entity_ID faceid, 
		       std::vector<Entity_ID> *nodeids) const;
    


  // Upward adjacencies
  //-------------------
    
  // Cells of type 'ptype' connected to a node
    
  void node_get_cells (const Entity_ID nodeid, 
		       const Parallel_type ptype,
		       std::vector<Entity_ID> *cellids) const;
    
  // Faces of type 'ptype' connected to a node
    
  void node_get_faces (const Entity_ID nodeid, 
		       const Parallel_type ptype,
		       std::vector<Entity_ID> *faceids) const;
    
  // Get faces of ptype of a particular cell that are connected to the
  // given node
    
  void node_get_cell_faces (const Entity_ID nodeid, 
			    const Entity_ID cellid,
			    const Parallel_type ptype,
			    std::vector<Entity_ID> *faceids) const;    
    
  // Cells connected to a face
    
  void face_get_cells (const Entity_ID faceid, 
		       const Parallel_type ptype,
		       std::vector<Entity_ID> *cellids) const;
    


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
			       std::vector<Entity_ID> *fadj_cellids) const;

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order

  void cell_get_node_adj_cells(const Entity_ID cellid,
			       const Parallel_type ptype,
			       std::vector<Entity_ID> *nadj_cellids) const;


    
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
			    std::vector<Entity_ID> *nodeids) const;
    
    
    
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
			 std::vector<Entity_ID> *entids) const; 

  void get_set_entities (const Set_Name setname, 
			 const Entity_kind kind, 
			 const Parallel_type ptype, 
			 std::vector<Entity_ID> *entids) const; 


  void get_set_entities (const char *setname, 
			 const Entity_kind kind, 
			 const Parallel_type ptype, 
			 std::vector<Entity_ID> *entids) const; 



  // this should be used with extreme caution:
  // modify coordinates
  void set_coordinate(unsigned int local_node_id,
		      double* source_begin, double* source_end);

private:
  void generate_(const GenerationSpec& g);
  void update_internals_();
  void clear_internals_();
  void build_maps_();


  Epetra_Map *cell_map_, *face_map_, *node_map_;


  std::vector<double> coordinates_;

  inline unsigned int node_index_(int i, int j, int k) const;
  inline unsigned int xyface_index_(int i, int j, int k) const;
  inline unsigned int yzface_index_(int i, int j, int k) const;
  inline unsigned int xzface_index_(int i, int j, int k) const;
  inline unsigned int cell_index_(int i, int j, int k) const;

  int nx_, ny_, nz_;  // number of cells in the three coordinate directions
  double x0_, x1_, y0_, y1_, z0_, z1_;  // coordinates of lower left front and upper right back of brick


  int num_cells_;
  int num_nodes_;
  int num_faces_;

  // Local-id tables of entities
  std::vector<unsigned int> cell_to_face_;
  std::vector<int> cell_to_face_dirs_;
  std::vector<unsigned int> cell_to_node_;
  std::vector<unsigned int> face_to_node_;
  std::vector<unsigned int> face_to_cell_;
  std::vector<unsigned int> node_to_face_;
  std::vector<unsigned int> node_to_cell_;

  // The following are mutable because they have to be modified 
  // after the class construction even though the class is instantiated
  // as a constant class

  mutable std::vector<std::vector<unsigned int> > side_sets_;
  mutable std::vector<std::vector<unsigned int> > element_blocks_;
  mutable std::vector<std::vector<unsigned int> > node_sets_;
  mutable std::vector<AmanziGeometry::RegionPtr> element_block_regions_;
  mutable std::vector<AmanziGeometry::RegionPtr> side_set_regions_;
  mutable std::vector<AmanziGeometry::RegionPtr> node_set_regions_;

  Epetra_Comm  *communicator_;

};

  // -------------------------
  // Template & inline members
  // ------------------------


  unsigned int Mesh_simple::node_index_(int i, int j, int k) const
  {
    return i + j*(nx_+1) + k*(nx_+1)*(ny_+1);
  }

  unsigned int Mesh_simple::cell_index_(int i, int j, int k) const
  {
    return i + j*nx_ + k*nx_*ny_;
  }

  unsigned int Mesh_simple::xyface_index_(int i, int j, int k) const
  {
    return i + j*nx_ + k*nx_*ny_;
  }

  unsigned int Mesh_simple::xzface_index_(int i, int j, int k) const
  {
    return i + j*nx_ + k*nx_*(ny_+1) +  xyface_index_(0,0,nz_+1);
  }

  unsigned int Mesh_simple::yzface_index_(int i, int j, int k) const
  {
    return i + j*(nx_+1) + k*(nx_+1)*ny_ + xzface_index_(0,0,nz_);
  }

} // close namespace AmanziMesh
} // close namespace Amanzi





#endif /* _MESH_MAPS_H_ */
