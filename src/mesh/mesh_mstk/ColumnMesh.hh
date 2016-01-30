/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_COLUMN_MESH_H_
#define AMANZI_COLUMN_MESH_H_

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Mesh_MSTK.hh"
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

class ColumnMesh : public Mesh {
      
 public:
      
  ColumnMesh(const Mesh& inmesh,
             const int column_id, 
             const VerboseObject *verbosity_obj = (VerboseObject *) NULL);

  ~ColumnMesh ();
  
  // reference for vis.
  const Mesh& vis_mesh() const { return extracted_; }

  // Get parallel type of entity - OWNED, GHOST, USED (See MeshDefs.hh)
  virtual
  Parallel_type entity_get_ptype(const Entity_kind kind,
          const Entity_ID entid) const {
    return OWNED;
  }


  // Parent entity in the source mesh if mesh was derived from another mesh
  virtual
  Entity_ID entity_get_parent(const Entity_kind kind, const Entity_ID entid) const {
    Entity_ID ent;
    switch (kind) {
      case FACE:
        ent = extracted_.entity_get_parent(kind, column_faces_[entid]);
        break;

      default:
        ent = extracted_.entity_get_parent(kind, entid);
        break;
    }
    return ent;
  }



  // Get cell type - UNKNOWN, TRI, QUAD, POLYGON, TET, PRISM, PYRAMID, HEX, POLYHED 
  // See MeshDefs.hh

  virtual
  Cell_type cell_get_type(const Entity_ID cellid) const {
    return extracted_.cell_get_type(cellid);
  }


  //
  // General mesh information
  // -------------------------
  //

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, USED)

  virtual
  unsigned int num_entities (const Entity_kind kind,
                             const Parallel_type ptype) const {
    int count;
    switch (kind) {
      case FACE:
        count = (ptype == GHOST) ? 0 : column_faces_.size();
        break;

      case BOUNDARY_FACE:
        count = 2;
        break;

      default:
        count = extracted_.num_entities(kind, ptype);
        break;
    }
    return count;
  }


  // Global ID of any entity
  virtual
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const {
    return lid;
  }



  //
  // Mesh Entity Adjacencies
  //-------------------------


  // Downward Adjacencies
  //---------------------

  // Get nodes of a cell

  virtual
  void cell_get_nodes (const Entity_ID cellid, 
                       Entity_ID_List *nodeids) const {
    extracted_.cell_get_nodes(cellid, nodeids);
  }


  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2

  virtual
  void face_get_nodes (const Entity_ID faceid,
                       Entity_ID_List *nodeids) const {
    extracted_.face_get_nodes(column_faces_[faceid], nodeids);
  }


  // Get nodes of edge
  
  virtual
  void edge_get_nodes (const Entity_ID edgeid, 
                       Entity_ID *nodeid0, Entity_ID *nodeid1) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }


  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node - The order of cells is
  // not guaranteed to be the same for corresponding nodes on
  // different processors

  virtual
  void node_get_cells (const Entity_ID nodeid,
                       const Parallel_type ptype,
                       Entity_ID_List *cellids) const {
    extracted_.node_get_cells(nodeid, ptype, cellids);
  }


  // Faces of type 'ptype' connected to a node - The order of faces is
  // not guarnateed to be the same for corresponding nodes on
  // different processors


  virtual
  void node_get_faces (const Entity_ID nodeid,
                       const Parallel_type ptype,
                       Entity_ID_List *faceids) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }


  // Get faces of ptype of a particular cell that are connected to the
  // given node - The order of faces is not guarnateed to be the same
  // for corresponding nodes on different processors

  virtual
  void node_get_cell_faces (const Entity_ID nodeid,
                            const Entity_ID cellid,
                            const Parallel_type ptype,
                            Entity_ID_List *faceids) const {
    Errors::Message mesg("Not implemented");
    amanzi_throw(mesg);
  }

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
          Entity_ID_List *fadj_cellids) const {
    extracted_.cell_get_face_adj_cells(cellid, ptype, fadj_cellids);
  }
    

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order

  virtual
  void cell_get_node_adj_cells(const Entity_ID cellid,
          const Parallel_type ptype,
          Entity_ID_List *nadj_cellids) const {
    extracted_.cell_get_node_adj_cells(cellid, ptype, nadj_cellids);
  }
    


  //
  // Mesh entity geometry
  //--------------
  //

  // Node coordinates - 3 in 3D and 2 in 2D

  virtual
  void node_get_coordinates (const Entity_ID nodeid,
                             AmanziGeometry::Point *ncoord) const {
    extracted_.node_get_coordinates(nodeid, ncoord);
  }


  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size divided by number of spatial dimensions

  virtual
  void face_get_coordinates (const Entity_ID faceid,
                             std::vector<AmanziGeometry::Point> *fcoords) const {
    extracted_.face_get_coordinates(column_faces_[faceid], fcoords);
  }

  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions

  virtual
  void cell_get_coordinates (const Entity_ID cellid,
                             std::vector<AmanziGeometry::Point> *ccoords) const {
    extracted_.cell_get_coordinates(cellid, ccoords);
  }


  //
  // Mesh modification
  //-------------------

  // Set coordinates of node
  virtual
  void node_set_coordinates (const Entity_ID nodeid,
                             const AmanziGeometry::Point ncoord) {
    extracted_.node_set_coordinates(nodeid, ncoord);
  }


  virtual
  void node_set_coordinates (const Entity_ID nodeid,
                             const double *ncoord) {
    extracted_.node_set_coordinates(nodeid, ncoord);
  }


  // Deform a mesh so that cell volumes conform as closely as possible
  // to target volumes without dropping below the minimum volumes.  If
  // move_vertical = true, nodes will be allowed to move only in the
  // vertical direction (right now arbitrary node movement is not allowed)
  // Nodes in any set in the fixed_sets will not be permitted to move.
  virtual
  int deform(const std::vector<double>& target_cell_volumes_in,
             const std::vector<double>& min_cell_volumes_in,
             const Entity_ID_List& fixed_nodes,
             const bool move_vertical) {
    return extracted_.deform(target_cell_volumes_in, min_cell_volumes_in,
                       fixed_nodes, move_vertical);
  }

  //
  // Epetra maps
  //------------
  virtual
  const Epetra_Map& cell_map (const bool include_ghost) const {
    return extracted_.cell_map(include_ghost);
  }

  virtual
  const Epetra_Map& face_map (const bool include_ghost) const {
    return *face_map_;
  }

  // dummy implementation so that frameworks can skip or overwrite
  const Epetra_Map& edge_map (const bool include_ghost) const 
  {
    Errors::Message mesg("Edges not implemented in this framework");
    amanzi_throw(mesg);
  }; 

  virtual
  const Epetra_Map& node_map (const bool include_ghost) const {
    return extracted_.node_map(include_ghost);
  }

  virtual
  const Epetra_Map& exterior_face_map (void) const {
    return *exterior_face_map_;
  }
    

  // Epetra importer that will allow apps to import values from a
  // Epetra vector defined on all owned faces into an Epetra vector
  // defined only on exterior faces
  
  virtual
  const Epetra_Import& exterior_face_importer (void) const {
    return *exterior_face_importer_;
  }


  //
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------
  //

  // Get number of entities of type 'category' in set
  virtual
  unsigned int get_set_size (const Set_ID setid,
                             const Entity_kind kind,
                             const Parallel_type ptype) const {

    return get_set_size(set_name_from_id(setid), kind, ptype);
  }


  virtual
  unsigned int get_set_size (const Set_Name setname,
                             const Entity_kind kind,
                             const Parallel_type ptype) const {
    Entity_ID_List ents;
    get_set_entities(setname, kind, ptype, &ents);
    return ents.size();
  }

  virtual
  unsigned int get_set_size (const char *setname,
                             const Entity_kind kind,
                             const Parallel_type ptype) const {
    std::string setname1(setname);
    return get_set_size(setname1,kind,ptype);
  }


  // Get list of entities of type 'category' in set

  virtual
  void get_set_entities (const Set_ID setid,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         Entity_ID_List *entids) const {
    get_set_entities(set_name_from_id(setid), kind, ptype, entids);
  }

  virtual
  void get_set_entities (const Set_Name setname,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         Entity_ID_List *entids) const {
    switch (kind) {
      case FACE: {
        Entity_ID_List faces;
        extracted_.get_set_entities(setname, kind, ptype, &faces);

        for (Entity_ID_List::const_iterator f=faces.begin();
             f!=faces.end(); ++f) {
          if (face_in_column_[*f] >= 0) entids->push_back(face_in_column_[*f]);
        }
        break;
      }

      default: {
        extracted_.get_set_entities(setname, kind, ptype, entids);
        break;
      }
    }
  }
    

  virtual
  void get_set_entities (const char *setname,
                         const Entity_kind kind,
                         const Parallel_type ptype,
                         Entity_ID_List *entids) const {
    std::string setname1(setname);
    get_set_entities(setname1,kind,ptype,entids);
  }


  // Miscellaneous functions

  virtual
  void write_to_exodus_file(const std::string filename) const {
    extracted_.write_to_exodus_file(filename);
  }


 protected:

  // get faces and face dirs of a cell. This can be called by
  // cell_get_faces_and_dirs method of the base class and the data
  // cached or it can be called directly by the
  // cell_get_faces_and_dirs method of this class
  virtual
  void cell_get_faces_and_dirs_internal (const Entity_ID cellid,
          Entity_ID_List *faceids,
          std::vector<int> *face_dirs,
          const bool ordered=false) const {

    faceids->resize(2);
    face_dirs->resize(2);
    
    // NOTE: the face directions with respect to the cell may be at
    // odds with how it is in the parent mesh but within this mesh its
    // consistent - so we think everything will work as it should
    Entity_ID_List faceids_extracted;
    std::vector<int> face_dirs_extracted;
    extracted_.cell_get_faces_and_dirs(cellid, &faceids_extracted,
            &face_dirs_extracted, ordered);

    int count = 0;
    for (int i=0; i!=faceids_extracted.size(); ++i) {
      if (face_in_column_[faceids_extracted[i]] >= 0) {      
        (*faceids)[count] = face_in_column_[faceids_extracted[i]];
        (*face_dirs)[count] = face_dirs_extracted[i];
        count++;
      }
    }
  }
    

  // Cells connected to a face - this function is implemented in each
  // mesh framework. The results are cached in the base class
  virtual
  void face_get_cells_internal (const Entity_ID faceid,
          const Parallel_type ptype,
          Entity_ID_List *cellids) const {
    extracted_.face_get_cells(column_faces_[faceid], ptype, cellids);
  }


  // edges of a face - this function is implemented in each mesh
  // framework. The results are cached in the base class
  virtual
  void face_get_edges_and_dirs_internal (const Entity_ID faceid,
					 Entity_ID_List *edgeids,
					 std::vector<int> *edge_dirs,
          const bool ordered=true) const {
    Errors::Message mesg("Not implemented");
    Exceptions::amanzi_throw(mesg);
  }

  // edges of a cell - this function is implemented in each mesh
  // framework. The results are cached in the base class. 
  virtual
  void cell_get_edges_internal (const Entity_ID cellid,
          Entity_ID_List *edgeids) const {
    Errors::Message mesg("Not implemented");
    Exceptions::amanzi_throw(mesg);
  }


  // edges and directions of a 2D cell - this function is implemented
  // in each mesh framework. The results are cached in the base class.
  virtual
  void cell_2D_get_edges_and_dirs_internal (const Entity_ID cellid,
                                            Entity_ID_List *edgeids,
          std::vector<int> *edge_dirs) const {
    Errors::Message mesg("Not implemented");
    Exceptions::amanzi_throw(mesg);
  }


  void build_epetra_maps_();
  void compute_special_node_coordinates_();


 private:

  Mesh const& parent_mesh_;
  Mesh_MSTK extracted_;
  int nfnodes_;
  int column_id_;
  Entity_ID_List column_faces_;
  Entity_ID_List face_in_column_;

  Epetra_Map *face_map_;
  Epetra_Map *exterior_face_map_;
  Epetra_Import *exterior_face_importer_;




};


} // close namespace AmanziMesh
} // close namespace Amanzi





#endif /* _MESH_MAPS_H_ */
