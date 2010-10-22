#ifndef _MESH_MAPS_H_
#define _MESH_MAPS_H_

// MOAB include files - for some reason they have to be first -
// if not the compiler cannot find MPI_COMM_WORLD

#include "MBmpi.h"
#include "MBTypes.h"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBReadUtilIface.hpp"
#include "MBParallelConventions.h"
#include "MBParallelComm.hpp"

#include "Entity_kind.hh"
#include "Data_structures.hh"
#include "Element_category.hh"
#include "dbc.hh"

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCPDecl.hpp>

#include <memory>


namespace MOAB_mesh
{
  
  class Mesh_maps
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


  // THIS IS THE API THAT THE OUTSIDE WORLD SEES

  public:

    Mesh_maps (const char *filename, MPI_Comm comm);

    void update ();

    // Local id interfaces.
    // --------------------

    template <typename IT>
    void cell_to_faces (unsigned int cell, IT begin, IT end) const;

    template <typename IT>
    void cell_to_face_dirs (unsigned int cell, IT begin, IT end) const;

    template <typename IT>
    void cell_to_nodes (unsigned int cell, IT begin, IT end) const;

    template <typename IT>
    void face_to_nodes (unsigned int face, IT begin, IT end) const;

    template <typename IT>
    void node_to_coordinates (unsigned int node, IT begin, IT end) const;

    template <typename IT>
    void face_to_coordinates (unsigned int face, IT begin, IT end) const;

    template <typename IT>
    void cell_to_coordinates (unsigned int call, IT begin, IT end) const;

    inline const Epetra_Map& cell_map (bool include_ghost) const;
    inline const Epetra_Map& face_map (bool include_ghost) const;
    inline const Epetra_Map& node_map (bool include_ghost) const;

    unsigned int count_entities (Mesh_data::Entity_kind kind, Element_Category category) const;

    // Entity Sets (cell, side, node)
    // ------------------------------
    unsigned int num_sets     (int kind) const;
    bool         valid_set_id (int kind, unsigned int id) const;
    unsigned int set_size     (unsigned int set_id, int kind, 
                               Element_Category category) const;

    template <typename IT>
    void set_ids (int kind, IT begin, IT end) const;

    template <typename IT>
    void get_set (unsigned int set_id, int kind, Element_Category category,
                  IT begin, IT end) const;

  };




  // The template implementations must be here because when we call 
  // them from some driver code, the compiler must be able to compile the
  // code for that particular data type



  template <typename IT>
  void Mesh_maps::cell_to_faces (unsigned int cellid, IT destination_begin, IT destination_end) const
  {
    MBEntityHandle cell;
    MBRange cell_faces;
    int *cell_faceids;
    int nf;

    cell = cell_id_to_handle[cellid];
      
    mbcore->get_adjacencies(&cell, 1, facedim, true, cell_faces, 
			    MBInterface::INTERSECT);

    nf = cell_faces.size();
    assert ((unsigned int) (destination_end - destination_begin) >= nf);
    cell_faceids = new int[nf];

    mbcore->tag_get_data(lid_tag,cell_faces,cell_faceids);

    std::copy (cell_faceids, cell_faceids+nf, destination_begin);

    delete [] cell_faceids;
  }


  template <typename IT>
  void Mesh_maps::cell_to_face_dirs (unsigned int cellid, IT destination_begin, IT destination_end) const
  {
    MBEntityHandle cell;
    MBRange cell_faces;
    int *cell_facedirs;
    int j,nf;


    cell = cell_id_to_handle[cellid];

    mbcore->get_adjacencies(&cell, 1, facedim, true, cell_faces, 
			    MBInterface::INTERSECT);

    nf = cell_faces.size();
    assert ((unsigned int) (destination_end - destination_begin) >= nf);
    cell_facedirs = new int[nf];

    j = 0;
    for (MBRange::iterator i = cell_faces.begin(); i != cell_faces.end(); i++) {
      MBEntityHandle face = *i;
      int sidenum, offset;
      mbcore->side_number(cell,face,sidenum,(cell_facedirs[j]),offset);
      j++;
    }
    
    std::copy (cell_facedirs, cell_facedirs+nf, destination_begin);

    delete [] cell_facedirs;
  }



  template <typename IT>
  void Mesh_maps::cell_to_nodes (unsigned int cellid, IT destination_begin, IT destination_end) const
  {
    MBEntityHandle cell;
    std::vector<MBEntityHandle> cell_nodes;
    int *cell_nodeids;
    int nn;


    cell = cell_id_to_handle[cellid];
      
    mbcore->get_connectivity(&cell, 1, cell_nodes);

    nn = cell_nodes.size();
    assert ((unsigned int) (destination_end - destination_begin) >= nn);
    cell_nodeids = new int[nn];
    
    for (int i = 0; i < nn; i++)
      mbcore->tag_get_data(lid_tag,&(cell_nodes[i]),1,&(cell_nodeids[i]));

    std::copy (cell_nodeids, cell_nodeids+nn, destination_begin);

    delete [] cell_nodeids;
  }



  template <typename IT>
  void Mesh_maps::face_to_nodes (unsigned int faceid, IT destination_begin, IT destination_end) const
  {
    MBEntityHandle face;
    std::vector<MBEntityHandle> face_nodes;
    int *face_nodeids;
    int nn;

    face = face_id_to_handle[faceid];

    mbcore->get_connectivity(&face, 1, face_nodes, true);

    nn = face_nodes.size();
    assert ((unsigned int) (destination_end - destination_begin) >= nn);

    face_nodeids = new int[nn];
    for (int i = 0; i < nn; i++) 
      mbcore->tag_get_data(lid_tag,&(face_nodes[i]),1,&(face_nodeids[i]));

    std::copy (face_nodeids, face_nodeids+nn, destination_begin);

    delete [] face_nodeids;
  }
  


  
  template <typename IT>
  void Mesh_maps::node_to_coordinates (unsigned int node_id, IT begin, IT end) const
  {
    MBEntityHandle node;

    assert ((unsigned int) (end-begin) >= spacedim);
    
    double coords[3];

    node = vtx_id_to_handle[node_id];

    mbcore->get_coords(&node, 1, coords);

    std::copy (coords, coords+spacedim, begin);

  }



  template <typename IT>
  void Mesh_maps::cell_to_coordinates (unsigned int cellid, IT destination_begin, IT destination_end) const
  {
    MBEntityHandle cell;
    std::vector<MBEntityHandle> cell_nodes;
    double *coords;
    int nn;


    cell = cell_id_to_handle[cellid];
      
    mbcore->get_connectivity(&cell, 1, cell_nodes);

    nn = cell_nodes.size();
    assert ((unsigned int) (destination_end - destination_begin) >= nn);

    coords = new double[spacedim*nn];

    for (int i = 0; i < nn; i++)
      mbcore->get_coords(&(cell_nodes[i]),1,coords+spacedim*i);

    std::copy (coords, coords+spacedim*nn, destination_begin);

    delete [] coords;
  }



  template <typename IT>
  void Mesh_maps::face_to_coordinates (unsigned int faceid, IT destination_begin, IT destination_end) const
  {
    MBEntityHandle face;
    std::vector<MBEntityHandle> face_nodes;
    double *coords;
    int nn;

    face = face_id_to_handle[faceid];

    mbcore->get_connectivity(&face, 1, face_nodes, true);

    nn = face_nodes.size();
    assert ((unsigned int) (destination_end - destination_begin) >= nn);

    coords = new double[spacedim*nn];
    
    for (int i = 0; i < nn; i++)
      mbcore->get_coords(&(face_nodes[i]),1,coords+spacedim*i);

    std::copy (coords, coords+spacedim*nn, destination_begin);

    delete [] coords;
  }
  


}

#endif /* _MESH_MAPS_H_ */
