#include "Mesh_maps.hh"
#include "dbc.hh"

using namespace std;
using namespace moab;


namespace MOAB_mesh
{

  // Constructor - load up mesh from file

  Mesh_maps::Mesh_maps (const char *filename, MPI_Comm comm)
  {
    int result, rank;

    clear_internals_();

    
    // Core MOAB object
    
    mbcore = new MBCore();

    if (comm) {
      // MOAB's parallel communicator
      int mbcomm_id;
      mbcomm = new ParallelComm(mbcore,comm,&mbcomm_id);

      if (!mbcomm) {
	cerr << "Failed to initialize MOAB communicator\n";
	assert(mbcomm == 0);
      }
    }

    epcomm = new Epetra_MpiComm(comm);


    if (!mbcomm || mbcomm->size() == 1) 
      serial_run = true;
    else
      serial_run = false;

    if (!serial_run) {
      
      // Load partitioned mesh - serial read of mesh with partition
      // info, deletion of non-local entities, resolution of
      // interprocessor connections. If we need ghosts we have to add
      // the option "PARALLEL_GHOSTS=A.B.C" where A is usually the
      // topological dimension of the mesh cells, B is the bridge
      // dimension or the dimension of entities across which we want
      // ghosts (0 for vertex connected ghost cells) and C indicates
      // the number of layers of ghost cells we want

      result = 
	mbcore->load_file(filename,NULL,
			"PARALLEL=READ_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION",
			NULL,NULL,NULL);
      
      rank = mbcomm->rank();
      
    }
    else {

      // Load serial mesh

      result =
	mbcore->load_file(filename,NULL,NULL,NULL,NULL,NULL);

      rank = 0;
    }

    if (result != MB_SUCCESS) {
      std::cerr << "FAILED" << std::endl;
      std::cerr << "Failed to load " << filename << " on processor " << rank << std::endl;
      std::cerr << "MOAB error code " << result << std::endl;
      assert(result == MB_SUCCESS);
    }
      
      

    // Dimension of space, mesh cells, faces etc
    
    result = mbcore->get_dimension(spacedim);
    

    // Highest topological dimension
    
    int nent;
    result = mbcore->get_number_entities_by_dimension(0,3,nent,false);
    if (nent) {
      celldim = 3;
      facedim = 2;
    }
    else {
      result = mbcore->get_number_entities_by_dimension(0,2,nent,false);
      if (nent) {
	celldim = 2;
	facedim = 1;
      }
      else {
	std::cerr << "Flow code works only on 2D and 3D meshes" << std::endl;
	assert(nent > 0);
      }
    }
      


    

    // Get all vertices on this processor 
    
    result = mbcore->get_entities_by_dimension(0,0,AllVerts,false);
    if (result != MB_SUCCESS) {
      std::cerr << "Could not get vertices" << std::endl;
      assert(result == MB_SUCCESS);
    }
    


    // Get all cells (faces in 2D, regions in 3D) on this processor
    
    result = mbcore->get_entities_by_dimension(0,celldim,AllCells,false);
    if (result != MB_SUCCESS) {
      std::cerr << "Could not get cells" << std::endl;
      assert(result == MB_SUCCESS);
    }
    
    



    // We have to make MOAB create the missing 'faces' (faces in 3D,
    // edges in 2D). We do this by looping over the cells and asking
    // for their faces with the create_if_missing=true option


    for (MBRange::iterator it = AllCells.begin(); it != AllCells.end(); it++) {
      MBEntityHandle cell = *it;
      MBRange cfaces;
      
      mbcore->get_adjacencies(&cell,1,facedim,true,cfaces,MBCore::UNION);
      if (result != MB_SUCCESS) {
	std::cerr << "Could not get faces of cell" << cell << std::endl;
	assert(result == MB_SUCCESS);
      }
    }
    

    // Get all "faces" (edges in 2D, faces in 3D) on this processor
    
    result = mbcore->get_entities_by_dimension(0,facedim,AllFaces,false);
    if (result != MB_SUCCESS) {
      std::cerr << "Could not get 'faces'" << std::endl;
      assert(result == MB_SUCCESS);
    }




    // Assign local IDs to entities


    result = mbcore->tag_create("LOCAL_ID",sizeof(unsigned int),
				MB_TAG_DENSE,MB_TYPE_INTEGER,lid_tag,
				0,true);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag handle for LOCAL_ID" << std::endl;
      assert(result == MB_SUCCESS);
    }
      
    nent = AllVerts.size();
    assert(nent > 0);
    int *lids = new int[nent];
    for (int i = 0; i < nent; i++) lids[i] = i;

    result = mbcore->tag_set_data(lid_tag,AllVerts,lids);

    delete [] lids;


    nent = AllFaces.size();
    assert(nent > 0);
    lids = new int[nent];
    for (int i = 0; i < nent; i++) lids[i] = i;

    result = mbcore->tag_set_data(lid_tag,AllFaces,lids);

    delete [] lids;


    nent = AllCells.size();
    assert(nent > 0);
    lids = new int[nent];
    for (int i = 0; i < nent; i++) lids[i] = i;

    result = mbcore->tag_set_data(lid_tag,AllCells,lids);

    delete [] lids;
      


    
    
    if (!serial_run) {
      // Ask Parallel Communicator to assign global IDs to entities

      bool largest_dim_only=false;
      int start_id=1;
      int largest_dim=celldim;
      result = mbcomm->assign_global_ids(0,largest_dim,start_id,largest_dim_only);
      if (result != MB_SUCCESS) {
	std::cerr << "Problem assigning global IDS" << std::endl;
	assert(result == MB_SUCCESS);
      }
    

      // Exchange global IDs across all processors

      result = mbcore->tag_get_handle("GLOBAL_ID",gid_tag);
      if (result != MB_SUCCESS) {
	std::cerr << "Could not get tag handle for GLOBAL_ID data" << std::endl;
	assert(result == MB_SUCCESS);
      }

      mbcomm->exchange_tags(gid_tag,AllVerts);
      mbcomm->exchange_tags(gid_tag,AllFaces);
      mbcomm->exchange_tags(gid_tag,AllCells);

      init_pvert_sets();
      init_pface_sets();
      init_pcell_sets();


      for (MBRange::iterator it = OwnedFaces.begin(); it != OwnedFaces.end(); ++it) {
	MBEntityHandle face = *it;
	int fgid;
	mbcore->tag_get_data(gid_tag,&face,1,&fgid);
      }
      
    }
    else {
      // Serial case - we assign global IDs ourselves

      result = mbcore->tag_create("GLOBAL_ID",sizeof(unsigned int),
				  MB_TAG_DENSE,MB_TYPE_INTEGER,gid_tag,
				  0,true);
      
      nent = AllVerts.size();
      int *gids = new int[nent];
      for (int i = 0; i < nent; i++) gids[i] = i+1;

      result = mbcore->tag_set_data(gid_tag,AllVerts,gids);

      delete [] gids;


      nent = AllFaces.size();
      gids = new int[nent];
      for (int i = 0; i < nent; i++) gids[i] = i+1;

      result = mbcore->tag_set_data(gid_tag,AllFaces,gids);

      delete [] gids;


      nent = AllCells.size();
      gids = new int[nent];
      for (int i = 0; i < nent; i++) gids[i] = i+1;

      result = mbcore->tag_set_data(gid_tag,AllCells,gids);

      delete [] gids;
      
    }


    // Create maps from IDs to MOAB entity handles

    init_id_handle_maps();
  

  }


  // Some initializations

  void Mesh_maps::clear_internals_ () 
  { 
    mbcore = NULL;
    mbcomm = NULL;

    AllVerts.clear();
    OwnedVerts.clear();
    NotOwnedVerts.clear();
    UsedVerts.clear();
    AllFaces.clear();
    OwnedFaces.clear();
    NotOwnedFaces.clear();
    UsedFaces.clear();
    AllCells.clear();
    OwnedCells.clear();
    GhostCells.clear();
    
    gid_tag = 0;

    spacedim = 3;
    celldim = -1;
    facedim = -1;
  }


  void Mesh_maps::init_id_handle_maps() {
    int i, nv, nf, nc;

    nv = AllVerts.size();

    vtx_id_to_handle.reserve(nv);

    i = 0;
    for (MBRange::iterator it = AllVerts.begin(); it != AllVerts.end(); ++it) {
      MBEntityHandle vtx = *it;
      vtx_id_to_handle[i++] = vtx;
    }
    


    nf = AllFaces.size();

    face_id_to_handle.reserve(nf);

    i = 0;
    for (MBRange::iterator it = AllFaces.begin(); it != AllFaces.end(); ++it) {
      MBEntityHandle face = *it;
      face_id_to_handle[i++] = face;
    }
    


    nc = AllCells.size();

    cell_id_to_handle.reserve(nc);

    i = 0;
    for (MBRange::iterator it = AllCells.begin(); it != AllCells.end(); ++it) {
      MBEntityHandle cell = *it;
      cell_id_to_handle[i++] = cell;
    }
    

  }


  void Mesh_maps::init_pvert_sets() {

    // Get not owned vertices

    mbcomm->get_pstatus_entities(0,PSTATUS_NOT_OWNED,NotOwnedVerts);

    // Subtract from all vertices on processor to get owned vertices only

    OwnedVerts = AllVerts;  // I think we DO want a data copy here
    OwnedVerts -= NotOwnedVerts;

    // Get ghost vertices (in MOAB, ghost vertices are not owned
    // by a processor and not connected to an element owned by
    // the processor

    MBRange GhostVerts;
    
    mbcomm->get_pstatus_entities(0,PSTATUS_GHOST,GhostVerts);
    
    UsedVerts = AllVerts; // I think we DO want a data copy here
    UsedVerts -= GhostVerts;
    
  }


  void Mesh_maps::init_pface_sets() {

    // Get not owned faces   

    mbcomm->get_pstatus_entities(facedim,PSTATUS_NOT_OWNED,NotOwnedFaces);
    
    // Subtract from all faces on processor to get owned faces only
    
    OwnedFaces = AllFaces;  // I think we DO want a data copy here
    OwnedFaces -= NotOwnedFaces;
    
    // Get ghost faces (in MOAB, ghost faces are not owned
    // by a processor and not connected to an element owned by
    // the processor
    
    MBRange GhostFaces;
    
    mbcomm->get_pstatus_entities(facedim,PSTATUS_GHOST,GhostFaces);
    
    UsedFaces = AllFaces; // I think we DO want a data copy here
    UsedFaces -= GhostFaces;
    
  }


  void Mesh_maps::init_pcell_sets() {
    // Get not owned cells (which is the same as ghost cells)

    mbcomm->get_pstatus_entities(facedim,PSTATUS_GHOST,GhostCells);
    
    // Subtract from all cells on processor to get owned cells only
    
    OwnedCells = AllCells;  // I think we DO want a data copy here
    OwnedCells -= GhostCells;
  }



  // Number of OWNED, GHOST or USED entities of different types

    
  unsigned int Mesh_maps::count_entities (Mesh_data::Entity_kind kind, 
					  Element_Category category) const
  {
    const int rank = (int) kind;
    const int index = ((int) category) - 1;

    switch (kind) {
    case Mesh_data::NODE:

      switch (category) {
      case OWNED:
	return !serial_run ? OwnedVerts.size() : AllVerts.size();
	break;
      case GHOST:
	return !serial_run ? NotOwnedVerts.size() : 0;
	break;
      case USED:
	return !serial_run ? UsedVerts.size() : AllVerts.size();
	break;
      default:
	return 0;
      }
      break;


    case Mesh_data::FACE:
      switch (category) {
      case OWNED:
	return !serial_run ? OwnedFaces.size() : AllFaces.size();
	break;
      case GHOST:
	return !serial_run ? NotOwnedFaces.size() : 0;
	break;
      case USED:
	return !serial_run ? UsedFaces.size() : AllFaces.size();
	break;
      default:
	return 0;
      }
      break;


    case Mesh_data::CELL:

      switch (category) {
      case OWNED:
	return !serial_run ? OwnedCells.size() : AllCells.size();
	break;
      case GHOST:
	return !serial_run ? GhostCells.size() : 0;
	break;
      case USED:
	return !serial_run ? OwnedCells.size() : AllCells.size();
	break;
      default:
	return 0;
      }

      break;
    default:
      std::cerr << "Count requested for unknown entity type" << std::endl;
    }
  }

  // Epetra map for nodes - basically a structure specifying the
  // global IDs of vertices owned or used by this processor

  const Epetra_Map& Mesh_maps::cell_map (bool include_ghost) const
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

    Epetra_Map *map = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);

    delete [] cell_gids;

    return *map;
  }




  // Epetra map for nodes - basically a structure specifying the
  // global IDs of vertices owned or used by this processor
  
  const Epetra_Map& Mesh_maps::face_map (bool include_ghost) const
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

	cout << "Allocated array of size " << OwnedFaces.size() << " on processor " << Mesh_maps::mbcomm->rank() << std::endl;

	mbcore->tag_get_data(gid_tag,OwnedFaces,face_gids);
	nface = OwnedFaces.size();
	
      }
    }
    else {
      
      face_gids = new int[AllFaces.size()];

      mbcore->tag_get_data(gid_tag,AllFaces,face_gids);
      nface = AllFaces.size();
    }


    cout << "Creating Epetra_map of size " << nface << " on processor " << epcomm->MyPID() << std::endl;

    for (int i = 0; i < nface; i++)
      cout << " " << face_gids[i] << " " << std::endl;

    Epetra_Map *map = new Epetra_Map(-1,nface,face_gids,0,*epcomm);

    cout << "FINISHED creating Epetra_map of size " << nface << " on processor " << epcomm->MyPID() << std::endl;

    

    delete [] face_gids;

    return *map;
  }




  // Epetra map for nodes - basically a structure specifying the
  // global IDs of vertices owned or used by this processor

  const Epetra_Map& Mesh_maps::node_map (bool include_ghost) const
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

    Epetra_Map *map = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);

    delete [] vert_gids;

    return *map;
  }

  
}

