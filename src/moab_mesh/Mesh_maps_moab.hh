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

    
    // Maps, Accessors and setters.
    // ----------------------------

    bool valid_entity_kind_ (int kind) const;


    void clear_internals_();

    void init_pvert_sets();
    void init_pface_sets();
    void init_pcell_sets();

    void init_id_handle_maps();

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

  // this should be used with extreme caution:
  // modify coordinates  
  void set_coordinate(unsigned int local_node_id, 
		      double* source_begin, double* source_end) 
  { throw std::exception(); };


  inline const Epetra_Map& cell_map (bool include_ghost) const;

  inline const Epetra_Map& face_map (bool include_ghost) const;

  inline const Epetra_Map& node_map (bool include_ghost) const;
  
  
};


// Epetra map for nodes - basically a structure specifying the
// global IDs of vertices owned or used by this processor

inline const Epetra_Map& Mesh_maps_moab::cell_map (bool include_ghost) const
{
  int *cell_gids;
  int ncell;

  if (!serial_run) {
    if (include_ghost) {

      // Put in owned cells before the ghost cells

      cell_gids = new int[OwnedCells.size()+GhostCells.size()];

      mbcore->tag_get_data(gid_tag,OwnedCells,cell_gids);
      ncell = OwnedCells.size();

      mbcore->tag_get_data(gid_tag,GhostCells,&(cell_gids[ncell]));
      ncell += GhostCells.size();

    }
    else {

      cell_gids = new int[OwnedCells.size()];

      mbcore->tag_get_data(gid_tag,OwnedCells,cell_gids);
      ncell = OwnedCells.size();
    }
  }
  else {
    cell_gids = new int[AllCells.size()];

    mbcore->tag_get_data(gid_tag,AllCells,cell_gids);
    ncell = AllCells.size();
  }

  for (int i = 0; i < ncell; i++) cell_gids[i] -= 1;

  Epetra_Map *map = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);

  delete [] cell_gids;

  return *map;
}




// Epetra map for nodes - basically a structure specifying the
// global IDs of vertices owned or used by this processor
  
inline const Epetra_Map& Mesh_maps_moab::face_map (bool include_ghost) const
{
  int *face_gids;
  int nface;

  if (!serial_run) {
    if (include_ghost) {
	
      // Put in owned faces before not owned faces

      face_gids = new int[OwnedFaces.size()+NotOwnedFaces.size()];

      mbcore->tag_get_data(gid_tag,OwnedFaces,face_gids);
      nface = OwnedFaces.size();

      mbcore->tag_get_data(gid_tag,NotOwnedFaces,&(face_gids[nface]));
      nface += NotOwnedFaces.size();
    }
    else {

      face_gids = new int[OwnedFaces.size()];

      mbcore->tag_get_data(gid_tag,OwnedFaces,face_gids);
      nface = OwnedFaces.size();
	
    }
  }
  else {
      
    face_gids = new int[AllFaces.size()];

    mbcore->tag_get_data(gid_tag,AllFaces,face_gids);
    nface = AllFaces.size();
  }


  for (int i = 0; i < nface; i++) face_gids[i] -= 1;

  Epetra_Map *map = new Epetra_Map(-1,nface,face_gids,0,*epcomm);


  delete [] face_gids;

  return *map;
}




// Epetra map for nodes - basically a structure specifying the
// global IDs of vertices owned or used by this processor

inline const Epetra_Map& Mesh_maps_moab::node_map (bool include_ghost) const
{
  int *vert_gids;
  int nvert;

  if (!serial_run) {
    if (include_ghost) {
	
      // Put in owned vertices before not owned vertices

      vert_gids = new int[OwnedVerts.size()+NotOwnedVerts.size()];

      mbcore->tag_get_data(gid_tag,OwnedVerts,vert_gids);
      nvert = OwnedVerts.size();

      mbcore->tag_get_data(gid_tag,NotOwnedVerts,&(vert_gids[nvert]));
      nvert += NotOwnedVerts.size();
    }
    else {
      vert_gids = new int[OwnedVerts.size()];    

      mbcore->tag_get_data(gid_tag,OwnedVerts,vert_gids);
      nvert = OwnedVerts.size();
    }
  }
  else {
    vert_gids = new int[AllVerts.size()];

    mbcore->tag_get_data(gid_tag,AllVerts,vert_gids);
    nvert = AllVerts.size();
  }

  for (int i = 0; i < nvert; i++) vert_gids[i] -= 1;

  Epetra_Map *map = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);

  delete [] vert_gids;

  return *map;
}


#endif /* _MESH_MAPS_MOAB_H_ */
