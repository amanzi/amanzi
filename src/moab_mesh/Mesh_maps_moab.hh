#ifndef _MESH_MAPS_MOAB_H_
#define _MESH_MAPS_MOAB_H_

// MOAB include files - for some reason they have to be first -
// if not the compiler cannot find MPI_COMM_WORLD

#include "MBmpi.h"
#include "MBTypes.h"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"
#include "MBReadUtilIface.hpp"
#include "MBParallelConventions.h"
#include "MBParallelComm.hpp"

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include "../mesh_data/Entity_kind.hh"

#include "../simple_mesh/Mesh_maps_base.hh"

#include <memory>
#include <vector>


class Mesh_maps_moab : public virtual Mesh_maps_base
{

  private:

    MBCore *mbcore;
    MBParallelComm *mbcomm;

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
    

    // These are MOAB's definitions of types of parallel mesh entities
    // These definitions are slightly different from what Amanzi has defined
    //
    // There are 2 types of cells - Owned and Ghost

    // There are 4 types of lower dimensional entities 
    //
    // 1. OWNED - owned by this processor
    //
    // 2. NOTOWNED - not owned by this processor
    //
    // 3. USED - connected to at least one cell owned by this
    // processor (may or may not be owned by this processor)
    //
    // 4. GHOST - neither the entity nor a cell connected to the
    // entity is owned by the processor

    // UNFORTUNATELY, THE TERMINOLOGY USED BY THE API USES GHOST TO
    // MEAN NOTOWNED

    // ALL = OWNED + NOTOWNED or USED + GHOST


    MBRange AllVerts, OwnedVerts, NotOwnedVerts, UsedVerts;

    MBRange AllFaces, OwnedFaces, NotOwnedFaces, UsedFaces;

    MBRange AllCells, OwnedCells, GhostCells;


    // tag handles

    MBTag lid_tag;        // Local ID
    MBTag gid_tag;        // Global ID
    MBTag mattag;         // Material tag
    MBTag sstag;          // Sideset tag
    MBTag nstag;          // Nodeset tag


    // Local ID to MOAB handle map

    std::vector<MBEntityHandle> vtx_id_to_handle;
    std::vector<MBEntityHandle> face_id_to_handle;
    std::vector<MBEntityHandle> cell_id_to_handle;


    // Maps

    Epetra_Map *cell_map_wo_ghosts_, *face_map_wo_ghosts_, *node_map_wo_ghosts_;
    Epetra_Map *cell_map_w_ghosts_, *face_map_w_ghosts_, *node_map_w_ghosts_;


    // Sets (material sets, sidesets, nodesets)
    // We store the number of sets in the whole problem regardless of whether
    // they are represented on this processor or not
    // We also store the IDs of the sets and the dimension of entities 
    // in those sets
  
    int nsets;
    int *setids, *setdims;
    
    // Private methods
    // ----------------------------

    bool valid_entity_kind_ (int kind) const;


    void clear_internals_();

    void init_pvert_lists();
    void init_pface_lists();
    void init_pcell_lists();

    void init_id_handle_maps();

    void init_cell_map();
    void init_face_map();
    void init_node_map();

    void init_set_info();

public:
  
  Mesh_maps_moab (const char *filename, MPI_Comm comm);
  
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


#endif /* _MESH_MAPS_MOAB_H_ */
