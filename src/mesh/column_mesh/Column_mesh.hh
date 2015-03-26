/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef _COLUMN_MESH_H_
#define _COLUMN_MESH_H_

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"
#include "Region.hh"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {

namespace AmanziMesh {

  // This is a mesh with a vertical column of prismatic cells (the
  // horizontal faces can be polygonal). Users of this class must note
  // several important assumptions and simplifications made in its
  // implementation.
  //
  // 1. A prerequisite of instantiation of this mesh class is that
  // columns should have been built in the base class
  //
  // 2. The base face of a column is considered to be perfectly
  // horizontal 
  //
  // 3. All cells in the column are assumed to have the
  // same prismatic topology (horizontal faces can have polygonal
  // topology, lateral faces are all quads)
  //
  // 4. The X-Y coordinates of nodes above the base face are adjusted
  // (only in the column mesh not the full 3D mesh) so that they are
  // perfectly stacked on top of the nodes of the base face.  The
  // Z-coordinates of the nodes of a face are adjusted so that they
  // match the Z-coordinate of the centroid of the original face
  //
  // 5. Adjustment of the node coordinates in the column results in
  // the face centroids and cell centroids of the column to line up
  // vertically; however, it will also result in small changes to the
  // cell volumes and face areas
  // 
  // 6. The lateral faces of the cells are ignored - so each cell is
  // considered to have only two faces, a bottom face and a top face.
  //
  // 7. The normals of the faces are pointing vertically down or up.
  //

class Column_mesh : public virtual Mesh
{
      
public:
      
  Column_mesh (Mesh& inmesh,
	       const int column_id, 
	       const VerboseObject *verbosity_obj = (VerboseObject *) NULL);

  ~Column_mesh () {};
  

  // Get parallel type of entity - does not make sense to have 
  // GHOST entities in columns
    
  Parallel_type entity_get_ptype(const Entity_kind kind, 
				 const Entity_ID entid) const {
    return OWNED; 
  }

  // Get parent entity 

  Entity_ID entity_get_parent(const Entity_kind kind,
                              const Entity_ID entid) const;

  // Get cell type
    
  Cell_type cell_get_type(const Entity_ID cellid) const {
    int parent_cellid = column_cells_[cellid];
    return parent_mesh_.cell_get_type(parent_cellid);
  }
   
  //
  // General mesh information
  // -------------------------
  //
    
  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, USED)
  //
  // Does not make sense to have GHOST entities in columns
    
  unsigned int num_entities (const Entity_kind kind,
			     const Parallel_type ptype) const {
    switch (kind) {
    case NODE: {
      int nfaces = column_faces_.size();
      return nfaces*nfnodes_;
    }

    case EDGE:
      return 0;

    case FACE:
      return (ptype == GHOST) ? 0 : column_faces_.size(); 

    case CELL:
      return (ptype == GHOST) ? 0 : column_cells_.size();
    }
  }
    
    
  // Global ID of any entity (same as local ID)
    
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const {
    return lid;
  }
    
    
    
  //
  // Mesh Entity Adjacencies 
  //-------------------------


  // Downward Adjacencies
  //---------------------
    
  // Get faces of a cell.

  // Get faces of a cell and directions in which the cell uses the face 

  // The Amanzi coding guidelines regarding function arguments is purposely
  // violated here to allow for a default input argument

  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order

  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise

  unsigned int cell_get_num_faces(const Entity_ID cellid) const {
    return 2;
  }

  void cell_get_faces (const Entity_ID cellid,
                       Entity_ID_List *faceids,
                       const bool ordered=false) const {
    faceids->resize(2);
    (*faceids)[0] = cellid;
    (*faceids)[1] = cellid+1;
  }


  void cell_get_faces_and_dirs (const Entity_ID cellid,
                                Entity_ID_List *faceids,
                                std::vector<int> *facedirs,
                                const bool ordered=false) const {

    // We don't want the base Mesh class to cache this data so we have
    // an direct call to cell_get_faces_and_dirs_internal instead of
    // letting the Mesh class call it

    cell_get_faces_and_dirs_internal(cellid, faceids, facedirs, ordered);
   }

  // Edges of a cell

  void cell_get_edges (const Entity_ID cellid,
                                Entity_ID_List *edgeids) const 
  { 
    Errors::Message mesg("Not implemented");
    Exceptions::amanzi_throw(mesg);
  }

   
    
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
		       std::vector<Entity_ID> *nodeids) const {

    // Cell i is bound by face i and face i+1
    // Each face has nfnodes_ nodes 

    int flo = cellid;
    int fhi = cellid+1;

    nodeids->resize(2*nfnodes_);

    for (int i = 0; i < nfnodes_; ++i) {
      (*nodeids)[i] = flo*nfnodes_ + i;
      (*nodeids)[i+nfnodes_] = fhi*nfnodes_ + i;
    }
  }
    
  // Get edges of a face and directions in which the face uses the edges 

  void face_get_edges_and_dirs (const Entity_ID faceid,
				Entity_ID_List *edgeids,
				std::vector<int> *edge_dirs,
				const bool ordered=false) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }

    
  // Get the local index of a face edge in a cell edge list

  void face_to_cell_edge_map(const Entity_ID faceid, const Entity_ID cellid,
			     std::vector<int> *map) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }


  // Get nodes of face 
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face 
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
    
  void face_get_nodes (const Entity_ID faceid, 
		       std::vector<Entity_ID> *nodeids) const {
    nodeids->resize(nfnodes_);
    for (int i = 0; i < nfnodes_; ++i)
      (*nodeids)[i] = faceid*nfnodes_ + i;
  }
    

  // Get nodes of edge

  void edge_get_nodes (const Entity_ID edgeid, Entity_ID *nodeid0,
		       Entity_ID *nodeid1) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }

  // Upward adjacencies
  //-------------------
    
  // Cells of type 'ptype' connected to a node
  // Since all cells of a column are on processor, ptype can be ignored
    
  void node_get_cells (const Entity_ID nodeid, 
		       const Parallel_type ptype,
		       std::vector<Entity_ID> *cellids) const {

    // get faceid from nodeid first

    int faceid = nodeid%nfnodes_;
    face_get_cells(faceid, ptype, cellids);

  }
    
  // Faces of type 'ptype' connected to a node
  // Since all faces of a column are on processor, ptype can be ignored
    
  void node_get_faces (const Entity_ID nodeid, 
		       const Parallel_type ptype,
		       std::vector<Entity_ID> *faceids) const {

    faceids->resize(1);
    (*faceids)[0] = nodeid%nfnodes_;

  }
    
  // Get faces of ptype of a particular cell that are connected to the
  // given node
    
  void node_get_cell_faces (const Entity_ID nodeid, 
			    const Entity_ID cellid,
			    const Parallel_type ptype,
			    std::vector<Entity_ID> *faceids) const {    

    faceids->resize(1);
    (*faceids)[0] = nodeid%nfnodes_;

  }

  // Cells connected to a face
    
  void face_get_cells (const Entity_ID faceid, 
                       const Parallel_type ptype,
                       std::vector<Entity_ID> *cellids) const {

    // We don't want the base Mesh class to cache this data so we have
    // an direct call to face_get_cells_internal instead of letting
    // the Mesh class call it

    face_get_cells_internal(faceid, ptype, cellids);

  }


  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell
  // Will return cell_below and cell_above in that order if they both exist
  // Otherwise it will only return the one that exists
  //
  // Since column meshes are all on processor, ptype is irrelevant

  void cell_get_face_adj_cells(const Entity_ID cellid,
			       const Parallel_type ptype,
			       std::vector<Entity_ID> *fadj_cellids) const {
    
    assert(cellid > 0 && cellid < column_cells_.size());

    if (cellid == 0) {
      fadj_cellids->resize(1);      
      (*fadj_cellids)[0] = cellid+1;
    }
    else if (cellid == column_cells_.size()-1) {
      fadj_cellids->resize(1);
      (*fadj_cellids)[0] = cellid-1;
    }
    else {
      fadj_cellids->resize(2);
      (*fadj_cellids)[0] = cellid-1;
      (*fadj_cellids)[1] = cellid+1;
    }

  }

  // Node connected neighboring cells of given cell
  // The cells are returned in no particular order
  //
  // Since column meshes are all on processor, ptype is irrelevant

  void cell_get_node_adj_cells(const Entity_ID cellid,
			       const Parallel_type ptype,
			       std::vector<Entity_ID> *nadj_cellids) const {

    // since its a column mesh face adjacent and node adjacent cell
    // list is the same

    cell_get_face_adj_cells(cellid, ptype, nadj_cellids);

  }

  // Special adjacency information for geological domains with a
  // semi-structured mesh. Will return -1 if there is no suitable cell
  // Since this a column mesh these just return the mesh information

  Entity_ID_List const & cells_of_column(const int columnID) const {
    amanzi_throw(Errors::Message("Meaningless to ask for column cells in a Column_mesh. Cell IDs go from 0 to ncells-1"));
  }

  Entity_ID_List const & faces_of_column(const int columnID) const {
    amanzi_throw(Errors::Message("Meaningless to ask for column faces in special Column_mesh. Face IDs go from 0 to nfaces-1"));
  }

  // Given a cell get its column ID

  int column_ID(const Entity_ID cellid) const {
    return 0; // only one column with ID 0
  }

  // Number of columns in mesh

  int num_columns() const {
    return 1;
  }
  
  Entity_ID cell_get_cell_above(const Entity_ID cellid) const {
    ASSERT(cellid >= 0 && cellid < column_cells_.size());
    return (cellid < column_cells_.size()-1 ? cellid+1 : -1);
  }


  Entity_ID cell_get_cell_below(const Entity_ID cellid) const {
    ASSERT(cellid >= 0 && cellid < column_cells_.size());
    return (cellid > 0 ? cellid-1 : -1);
  }

  Entity_ID node_get_node_above(const Entity_ID nodeid) const {
    ASSERT(nodeid >= 0 && nodeid < nfnodes_*column_faces_.size());
    int faceid = nodeid%nfnodes_;
    return (faceid < column_faces_.size()-1 ? (faceid+1)*nfnodes_+nodeid : -1);
  }

    
  //
  // Mesh entity geometry
  //--------------
  //
    
  // Node coordinates - 3 in 3D and 2 in 2D
    
  void node_get_coordinates (const Entity_ID nodeid, 
			     AmanziGeometry::Point *ncoord) const {
    *ncoord = node_coordinates_[nodeid];
  }
    
    
  // Face coordinates - conventions same as face_to_nodes call 
  // Number of nodes is the vector size divided by number of spatial dimensions
    
  void face_get_coordinates (const Entity_ID faceid, 
			     std::vector<AmanziGeometry::Point> *fcoords) const {
    fcoords->clear();
    for (int i = 0; i < nfnodes_; ++i)
      fcoords->push_back(node_coordinates_[faceid*nfnodes_+i]);
  }
    
  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions
    
  void cell_get_coordinates (const Entity_ID cellid, 
			     std::vector<AmanziGeometry::Point> *ccoords) const {
    ccoords->clear();
    for (int i = 0; i < nfnodes_; ++i)
      ccoords->push_back(node_coordinates_[cellid*nfnodes_+i]);
    for (int i = 0; i < nfnodes_; ++i)
      ccoords->push_back(node_coordinates_[(cellid+1)*nfnodes_+i]);    
  }

    
  // Modify the coordinates of a node

  void node_set_coordinates (const Entity_ID nodeid, const AmanziGeometry::Point coords) {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }

  void node_set_coordinates (const Entity_ID nodeid, const double *coords) {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }
    
  //
  // Epetra maps
  //------------
    
    
  const Epetra_Map& cell_map (bool include_ghost) const {
    return *cell_map_;
  }
    
  const Epetra_Map& face_map (bool include_ghost) const {
    return *face_map_;
  }

  const Epetra_Map& node_map (bool include_ghost) const {
    return *node_map_;
  }
    
  const Epetra_Map& exterior_face_map (void) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }
    
  // Epetra importer that will allow apps to import values from a
  // Epetra vector defined on all owned faces into an Epetra vector
  // defined only on exterior faces
  
  const Epetra_Import& exterior_face_importer (void) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }
    
    
  //
  // Boundary Conditions or Sets
  //----------------------------
  //
    

  // Get number of entities of type 'category' in set

  unsigned int get_set_size (const Set_ID setid, 
			     const Entity_kind kind,
			     const Parallel_type ptype) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }



  unsigned int get_set_size (const Set_Name setname, 
			     const Entity_kind kind,
			     const Parallel_type ptype) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }


  unsigned int get_set_size (const char *setname, 
			     const Entity_kind kind,
			     const Parallel_type ptype) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }


  // Get list of entities of type 'category' in set

  void get_set_entities (const Set_ID setid, 
			 const Entity_kind kind, 
			 const Parallel_type ptype, 
			 Entity_ID_List *entids) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }

  void get_set_entities (const Set_Name setname, 
			 const Entity_kind kind, 
			 const Parallel_type ptype, 
			 Entity_ID_List *entids) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }


  void get_set_entities (const char *setname, 
			 const Entity_kind kind, 
			 const Parallel_type ptype, 
			 Entity_ID_List *entids) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }


  // Miscellaneous

  void write_to_exodus_file (const std::string filename) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }



  // this should be used with extreme caution:
  // modify coordinates
  void set_coordinate(Entity_ID local_node_id,
		      double* source_begin, double* source_end) {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }


  // Deform a mesh so that cell volumes conform as closely as possible
  // to target volumes without dropping below the minimum volumes.  If
  // move_vertical = true, nodes will be allowed to move only in the
  // vertical direction (right now arbitrary node movement is not allowed)
  
  int deform(const std::vector<double>& target_cell_volumes_in, 
             const std::vector<double>& min_cell_volumes_in, 
             const Entity_ID_List& fixed_nodes,
             const bool move_vertical) {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }
  

private:
  Mesh const& parent_mesh_;
  int nfnodes_;
  int column_id_;
  Entity_ID_List const& column_cells_;
  Entity_ID_List const& column_faces_;
  std::vector< AmanziGeometry::Point > node_coordinates_;

  void build_epetra_maps_();
  void compute_special_node_coordinates_();

  // get faces and face dirs of a cell. This can be called by
  // cell_get_faces_and_dirs method of the base class and the data
  // cached or it can be called directly by the
  // cell_get_faces_and_dirs method of this class

  void cell_get_faces_and_dirs_internal (const Entity_ID cellid,
                                         Entity_ID_List *faceids,
                                         std::vector<int> *facedirs,
                                         const bool ordered=false) const {
    faceids->resize(2);
    facedirs->resize(2);
    
    // NOTE: the face directions with respect to the cell may be at
    // odds with how it is in the parent mesh but within this mesh its
    // consistent - so we think everything will work as it should


    (*faceids)[0] = cellid;
    (*facedirs)[0] = -1;
    (*faceids)[1] = cellid+1;
    (*facedirs)[1] = 1;

  }


  // Get cells attached to a face. This can be called by
  // face_get_cells method of the base class and the data cached or it
  // can be directly called by the face_get_cells method of this class
  

  void face_get_cells_internal (const Entity_ID faceid,
                                const Parallel_type ptype,
                                Entity_ID_List *cellids) const {

    if (faceid == 0) { // bottom face
      cellids->resize(1);
      (*cellids)[0] = faceid;
    }
    else if (faceid == column_faces_.size()-1) { // top face
      cellids->resize(1);
      (*cellids)[0] = faceid-1;
    }
    else {
      cellids->resize(2);
      (*cellids)[0] = faceid-1;
      (*cellids)[1] = faceid;
    }

  }

  void face_get_edges_and_dirs_internal (const Entity_ID faceid,
					 Entity_ID_List *edgeids,
					 std::vector<int> *edge_dirs,
					 const bool ordered=true) const {
    amanzi_throw(Errors::Message("Not implemented"));
  }

  void cell_get_edges_internal (const Entity_ID cellid,
				Entity_ID_List *edgeids) const {
    amanzi_throw(Errors::Message("Not implemented"));
  }

  Epetra_Map *cell_map_, *face_map_, *node_map_;

};


} // close namespace AmanziMesh
} // close namespace Amanzi





#endif /* _MESH_MAPS_H_ */
