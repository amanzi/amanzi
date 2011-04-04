#ifndef _MESH_MAPS_MSTK_H_
#define _MESH_MAPS_MSTK_H_


#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include "Entity_kind.hh"

#include "Mesh_maps_base.hh"

#include <memory>
#include <vector>
#include <sstream>

#include "MSTK.h"



class Mesh_maps_mstk : public virtual Mesh_maps_base
{

  private:

  MPI_Comm comm;
  int myprocid, numprocs;

  Mesh_ptr mesh;
  //    MSTKComm *mstkcomm; Not defined - will define if we need it

  int serial_run;
  
  int spacedim;
  int celldim; // Topological dimension of highest level entities
  int facedim; // Topological dimension of 2nd highest level entities
  
  Epetra_MpiComm *epcomm;



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


  // Local ID to MSTK handle map

  std::vector<MEntity_ptr> vtx_id_to_handle;
  std::vector<MEntity_ptr> face_id_to_handle;
  std::vector<MEntity_ptr> cell_id_to_handle;


  // Maps
  
  Epetra_Map *cell_map_wo_ghosts_, *face_map_wo_ghosts_, *node_map_wo_ghosts_;
  Epetra_Map *cell_map_w_ghosts_, *face_map_w_ghosts_, *node_map_w_ghosts_;
  
  
  // Sets (material sets, sidesets, nodesets)
  // We store the number of sets in the whole problem regardless of whether
  // they are represented on this processor or not
  // We also store the IDs of the sets and the dimension of entities 
  // in those sets
  
  // We could also store a single array of setids and another array of setdims
  // like we do for Mesh_Maps_Moab. Some code is easier this way and some code
  // easier the other way

  // Cannot use std::vector<int> because we cannot pass it into MPI routines

  int nmatsets, nsidesets, nnodesets;
  int *matset_ids, *sideset_ids, *nodeset_ids; 
                                              
  
  
  // flag whether to flip a face dir or not when returning nodes of a face
  
  bool *faceflip;
  
  // Private methods
  // ----------------------------
  
  bool valid_entity_kind_ (int kind) const;
  
  
  void clear_internals_();
  
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

public:
  
  Mesh_maps_mstk (const char *filename, MPI_Comm comm);
  ~Mesh_maps_mstk();
  
  void update ();
  
  // Local id interfaces.
  // --------------------
  
  void cell_to_faces (unsigned int cell, 
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end);
  void cell_to_faces (unsigned int cell, 
		      unsigned int * begin, unsigned int * end);

  void cell_to_face_dirs (unsigned int cell, 
			 std::vector<int>::iterator begin, 
			 std::vector<int>::iterator end);
  void cell_to_face_dirs (unsigned int cell, 
			  int * begin, int * end);


  void cell_to_nodes (unsigned int cell, 
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end);
  void cell_to_nodes (unsigned int cell, 
		      unsigned int * begin, unsigned int * end);


  void face_to_nodes (unsigned int face,
		      std::vector<unsigned int>::iterator begin, 
		      std::vector<unsigned int>::iterator end);
  void face_to_nodes (unsigned int face,
		      unsigned int * begin, unsigned int * end);
  

  void node_to_coordinates (unsigned int node, 
			    std::vector<double>::iterator begin,
			    std::vector<double>::iterator end);
  void node_to_coordinates (unsigned int node, 
			    double * begin,
			    double * end);
  
  void face_to_coordinates (unsigned int face,
			    std::vector<double>::iterator begin, 
			    std::vector<double>::iterator end);
  void face_to_coordinates (unsigned int face,
			    double * begin, 
			    double * end);
  
  void cell_to_coordinates (unsigned int cell, 
			    std::vector<double>::iterator begin,
			    std::vector<double>::iterator end);
  void cell_to_coordinates (unsigned int cell, 
			    double * begin,
			    double * end);
  
  unsigned int count_entities (Mesh_data::Entity_kind kind,
			       Element_Category category) const;

  unsigned int num_sets(Mesh_data::Entity_kind kind) const;
  
  unsigned int get_set_size (unsigned int set_id, 
			     Mesh_data::Entity_kind kind,
			     Element_Category category) const;

  unsigned int get_set_size (const char* name,
			     Mesh_data::Entity_kind kind,
			     Element_Category category) const;
  

  void get_set_ids (Mesh_data::Entity_kind kind, 
		    std::vector<unsigned int>::iterator begin, 
		    std::vector<unsigned int>::iterator end) const;
  void get_set_ids (Mesh_data::Entity_kind kind, 
		    unsigned int * begin, 
		    unsigned int * end) const;
 
  bool valid_set_id (unsigned int id, Mesh_data::Entity_kind kind) const;

  void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
		Element_Category category,
		std::vector<unsigned int>::iterator begin,
		std::vector<unsigned int>::iterator end) const;
  void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, 
		Element_Category category,
		unsigned int * begin,
		unsigned int * end) const;

  inline const Epetra_Comm * get_comm() { return epcomm; };

  unsigned int GID(unsigned int lid, Mesh_data::Entity_kind kind);

  // this should be used with extreme caution:
  // modify coordinates  
  void set_coordinate(unsigned int local_node_id, 
		      double* source_begin, double* source_end) 
  { throw std::exception(); };


  inline const Epetra_Map& cell_map (bool include_ghost) const {
    if (serial_run)
      return *cell_map_wo_ghosts_;
    else
      return (include_ghost ? *cell_map_w_ghosts_ : *cell_map_wo_ghosts_);
  }
	    

  inline const Epetra_Map& face_map (bool include_ghost) const {
    if (serial_run)
      return *face_map_wo_ghosts_;
    else
      return (include_ghost ? *face_map_w_ghosts_ : *face_map_wo_ghosts_);
  }

  inline const Epetra_Map& node_map (bool include_ghost) const {
    if (serial_run)
      return *node_map_wo_ghosts_;
    else
      return (include_ghost ? *node_map_w_ghosts_ : *node_map_wo_ghosts_);
  }
  
  
};


#endif /* _MESH_MAPS_MSTK_H_ */
