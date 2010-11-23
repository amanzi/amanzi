#include "Mesh_maps_moab.hh"
//#include <Teuchos_RCP.hpp>

using namespace std;
using namespace moab;


  // Constructor - load up mesh from file

Mesh_maps_moab::Mesh_maps_moab (const char *filename, MPI_Comm comm)
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
    // the option "PARALLEL_GHOSTS=A.B.C.D" where A is usually the
    // topological dimension of the mesh cells, B is the bridge
    // dimension or the dimension of entities across which we want
    // ghosts (0 for vertex connected ghost cells) and C indicates
    // the number of layers of ghost cells we want, D indicates if
    // we want edges/faces bounding the ghost cells as well (1 for
    // edges, 2 for faces, 3 for faces and edges)

    // In the specification for the Ghosts we made the assumption 
    // that we are dealing with 3D meshes only

    result = 
      mbcore->load_file(filename,NULL,
			"PARALLEL=READ_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;PARALLEL_GHOSTS=3.0.1.2",
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
      



  { // Keep together and in this order 

    init_pvert_lists();
    init_pcell_lists(); // cells MUST be initialized before faces
    init_pface_lists();

    
    // Create maps from local IDs to MOAB entity handles (must be after
    // the various init_p*_list calls)
    
    init_id_handle_maps();


  }
    
    
  init_global_ids();


  init_pface_dirs();


  // Create Epetra_maps
  
  init_cell_map();
  init_face_map();
  init_node_map();
  


  // Initialize some info about the global number of sets, global set
  // IDs and set types

  init_set_info();

}


// Some initializations

void Mesh_maps_moab::clear_internals_ () 
{ 
  mbcore = NULL;
  mbcomm = NULL;

  AllVerts.clear();
  OwnedVerts.clear();
  NotOwnedVerts.clear();
  AllFaces.clear();
  OwnedFaces.clear();
  NotOwnedFaces.clear();
  AllCells.clear();
  OwnedCells.clear();
  GhostCells.clear();
    
  lid_tag = 0;
  gid_tag = 0;
  mattag = 0;
  sstag = 0;
  nstag = 0;

  spacedim = 3;
  celldim = -1;
  facedim = -1;

  faceflip = NULL;

  cell_map_w_ghosts_ = cell_map_wo_ghosts_ = NULL;
  face_map_w_ghosts_ = face_map_wo_ghosts_ = NULL;
  node_map_w_ghosts_ = node_map_wo_ghosts_ = NULL;

  nsets = 0;
  setids = setdims = NULL;
}


void Mesh_maps_moab::init_id_handle_maps() {
  int i, nv, nf, nc;
  int result;

  // Assign local IDs to entities


  result = mbcore->tag_create("LOCAL_ID",sizeof(unsigned int),
			      MB_TAG_DENSE,MB_TYPE_INTEGER,lid_tag,
			      0,true);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting tag handle for LOCAL_ID" << std::endl;
    assert(result == MB_SUCCESS);
  }
      


  nv = AllVerts.size();

  vtx_id_to_handle.reserve(nv);

  i = 0;
  for (MBRange::iterator it = OwnedVerts.begin(); it != OwnedVerts.end(); ++it) {
    MBEntityHandle vtx = *it;
    mbcore->tag_set_data(lid_tag,&vtx,1,&i);
    vtx_id_to_handle[i++] = vtx;
  }    
  for (MBRange::iterator it = NotOwnedVerts.begin(); it != NotOwnedVerts.end(); ++it) {
    MBEntityHandle vtx = *it;
    mbcore->tag_set_data(lid_tag,&vtx,1,&i);
    vtx_id_to_handle[i++] = vtx;
  }
    


  nf = AllFaces.size();

  face_id_to_handle.reserve(nf);

  i = 0;
  for (MBRange::iterator it = OwnedFaces.begin(); it != OwnedFaces.end(); ++it) {
    MBEntityHandle face = *it;
    mbcore->tag_set_data(lid_tag,&face,1,&i);
    face_id_to_handle[i++] = face;
  }
  for (MBRange::iterator it = NotOwnedFaces.begin(); it != NotOwnedFaces.end(); ++it) {
    MBEntityHandle face = *it;
    mbcore->tag_set_data(lid_tag,&face,1,&i);
    face_id_to_handle[i++] = face;
  }
    


  nc = AllCells.size();

  cell_id_to_handle.reserve(nc);

  i = 0;
  for (MBRange::iterator it = OwnedCells.begin(); it != OwnedCells.end(); ++it) {
    MBEntityHandle cell = *it;
    mbcore->tag_set_data(lid_tag,&cell,1,&i);
    cell_id_to_handle[i++] = cell;
  }
  for (MBRange::iterator it = GhostCells.begin(); it != GhostCells.end(); ++it) {
    MBEntityHandle cell = *it;
    mbcore->tag_set_data(lid_tag,&cell,1,&i);
    cell_id_to_handle[i++] = cell;
  }
    

}



void Mesh_maps_moab::init_global_ids() {
  int result;
    
  if (!serial_run) {
    // Ask Parallel Communicator to assign global IDs to entities

    bool largest_dim_only=false;
    int start_id=0;
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

  }
  else {
    // Serial case - we assign global IDs ourselves

    result = mbcore->tag_create("GLOBAL_ID",sizeof(unsigned int),
				MB_TAG_DENSE,MB_TYPE_INTEGER,gid_tag,
				0,true);
      
    int nent = AllVerts.size();
    int *gids = new int[nent];
    for (int i = 0; i < nent; i++) gids[i] = i;

    result = mbcore->tag_set_data(gid_tag,AllVerts,gids);

    delete [] gids;


    nent = AllFaces.size();
    gids = new int[nent];
    for (int i = 0; i < nent; i++) gids[i] = i;

    result = mbcore->tag_set_data(gid_tag,AllFaces,gids);

    delete [] gids;


    nent = AllCells.size();
    gids = new int[nent];
    for (int i = 0; i < nent; i++) gids[i] = i;

    result = mbcore->tag_set_data(gid_tag,AllCells,gids);

    delete [] gids;
      
  }


}


void Mesh_maps_moab::init_pvert_lists() {
  int result;

  // Get all vertices on this processor 
    
  result = mbcore->get_entities_by_dimension(0,0,AllVerts,false);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get vertices" << std::endl;
    assert(result == MB_SUCCESS);
  }    


  // Get not owned vertices

  result = mbcomm->get_pstatus_entities(0,PSTATUS_NOT_OWNED,NotOwnedVerts);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get NotOwned vertices" << std::endl;
    assert(result == MB_SUCCESS);
  }

  // Subtract from all vertices on processor to get owned vertices only

  OwnedVerts = AllVerts;  // I think we DO want a data copy here
  OwnedVerts -= NotOwnedVerts;

}


// init_pface_lists is more complicated than init_pvert_lists and
// init_pcell_lists because of the way MOAB is setting up shared
// entities and ghost entities. When we ask MOAB to resolve shared
// entities, then MOAB sets up faces on interprocessor boundaries and
// assigns each of them to some processor. Therefore, the pstatus tags
// on these faces are correctly set. On the other hand when we ask for
// ghost cells, MOAB does not automatically create ghost faces. Also,
// when we go through ghost cells and create their faces, MOAB does
// not tag them as ghost faces, tagging them as owned faces
// instead. So we have to process them specially.


void Mesh_maps_moab::init_pface_lists() {
  int result;


  // Make MOAB create the missing 'faces' (faces in 3D, edges in
  // 2D). We do this by looping over the cells and asking for their
  // faces with the create_if_missing=true option


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


  // Get not owned faces 

  result = mbcomm->get_pstatus_entities(facedim,PSTATUS_NOT_OWNED,NotOwnedFaces);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get NotOwned 'faces'" << std::endl;
    assert(result == MB_SUCCESS);
  }


  // Subtract from all faces on processor to get owned faces only
    
  OwnedFaces = AllFaces;  // I think we DO want a data copy here
  OwnedFaces -= NotOwnedFaces;
    

}


void Mesh_maps_moab::init_pface_dirs() {
  int result, zero=0, minus1=-1;
  int face_lid, face_gid, cell_gid;
  int sidenum, offset, facedir;
  MBEntityHandle cell, face;
  int DebugWait=1;

  // Do some additional processing to see if ghost faces and their masters
  // are oriented the same way; if not, turn on flag to flip them


  /* In this code, we increment local values of global IDs by 1 so
     that we can distinguish between the lowest gid and no data */

  //  result = mbcore->tag_create("TMP_FC_TAG",sizeof(MBEntityHandle),
  //			      MB_TAG_DENSE, MB_TYPE_HANDLE, tmp_fc_tag,
  //			      &zero,true);

  MBTag tmp_fc0_tag, tmp_fc1_tag;
  result = mbcore->tag_create("TMP_FC0_TAG",sizeof(int),
			      MB_TAG_DENSE, MB_TYPE_INTEGER, tmp_fc0_tag,
			      &zero,true);
  result = mbcore->tag_create("TMP_FC1_TAG",sizeof(int),
			      MB_TAG_DENSE, MB_TYPE_INTEGER, tmp_fc1_tag,
			      &zero,true);

  
  for (MBRange::iterator it = OwnedFaces.begin(); it != OwnedFaces.end(); it++) {
    MBRange fcells;
    face = *it;

    result = mbcore->get_adjacencies(&face,1,celldim,false,fcells,MBCore::UNION);
    if (result != MB_SUCCESS) {
      cout << "Could not get cells of face" << std::endl;
      assert(result == MB_SUCCESS);
    }

    mbcore->tag_set_data(tmp_fc0_tag,&face,1,&zero);
    mbcore->tag_set_data(tmp_fc1_tag,&face,1,&zero);


    for (MBRange::iterator jt = fcells.begin(); jt != fcells.end(); ++jt) {

      cell = *jt;
      result = mbcore->side_number(cell,face,sidenum,facedir,offset);
      if (result != MB_SUCCESS) {
	cout << "Could not get face dir w.r.t. cell" << std::endl;
	assert(result == MB_SUCCESS);
      }

      result = mbcore->tag_get_data(gid_tag,&cell,1,&cell_gid);
      cell_gid += 1; 

      if (facedir == 1)
	result = mbcore->tag_set_data(tmp_fc0_tag,&face,1,&cell_gid);
      else
	result = mbcore->tag_set_data(tmp_fc1_tag,&face,1,&cell_gid);
      
    }

  }



  mbcomm->exchange_tags(tmp_fc0_tag,AllFaces);
  if (result != MB_SUCCESS) {
    cout << "Could not get exchange tag data successfully" << std::endl;
    assert(result == MB_SUCCESS);
  }

  mbcomm->exchange_tags(tmp_fc1_tag,AllFaces);
  if (result != MB_SUCCESS) {
    cout << "Could not get exchange tag data successfully" << std::endl;
    assert(result == MB_SUCCESS);
  }



  faceflip = new bool[AllFaces.size()];
  for (int i = 0; i < AllFaces.size(); i++) faceflip[i] = false;

  for (MBRange::iterator it = NotOwnedFaces.begin(); it != NotOwnedFaces.end(); it++) {
    MBRange fcells;
    int ghost_cell0_gid = 0, ghost_cell1_gid = 0;
    int master_cell0_gid = 0, master_cell1_gid = 0;

    face = *it;


    mbcore->tag_get_data(tmp_fc0_tag,&face,1,&master_cell0_gid);
    if (result != MB_SUCCESS) {
      cout << "Could not get face tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
    mbcore->tag_get_data(tmp_fc1_tag,&face,1,&master_cell1_gid);
    if (result != MB_SUCCESS) {
      cout << "Could not get face tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }


    mbcore->get_adjacencies(&face,1,celldim,false,fcells,MBCore::UNION);
    if (result != MB_SUCCESS) {
      cout << "Could not get cells of face" << std::endl;
      assert(result == MB_SUCCESS);
    }


    for (MBRange::iterator jt = fcells.begin(); jt != fcells.end(); ++jt) {
      cell = *jt;

      result = mbcore->side_number(cell,face,sidenum,facedir,offset);
      if (result != MB_SUCCESS) {
	cout << "Could not get face dir w.r.t. cell" << std::endl;
	assert(result == MB_SUCCESS);
      }

      if (facedir == 1) {
	result = mbcore->tag_get_data(gid_tag,&cell,1,&ghost_cell0_gid);
	ghost_cell0_gid += 1;
      }
      else {
	result = mbcore->tag_get_data(gid_tag,&cell,1,&ghost_cell1_gid);
	ghost_cell1_gid += 1;
      }
    }

    if (ghost_cell0_gid == master_cell1_gid || 
	ghost_cell1_gid == master_cell0_gid) {

      // Both cells don't have to match because a ghost face may 
      // not have the cell on the other side
      
      result = mbcore->tag_get_data(lid_tag,&face,1,&face_lid);
      if (result != MB_SUCCESS) {
	cout << "Could not get face tag data" << std::endl;
	assert(result == MB_SUCCESS);
      }
      faceflip[face_lid] = true;

    }
    else { // Sanity check
      if (ghost_cell0_gid != master_cell0_gid &&
	  ghost_cell1_gid != master_cell1_gid) {

	// Problem if there is no match at all

	//
	mbcore->tag_get_data(gid_tag,&face,1,&face_gid);
	//

	cout << "Face cells mismatch between master and ghost (processor " << mbcomm->rank() << ")" << std::endl;
	cout << " Face " << face_gid << std::endl;
	cout << "Master cells " << master_cell0_gid << " " << master_cell1_gid << std::endl;
	cout << "Ghost cells " << ghost_cell0_gid << " " << ghost_cell1_gid << std::endl;
      }
    }
  }

}


void Mesh_maps_moab::init_pcell_lists() {
  int result;

  // Get all cells (faces in 2D, regions in 3D) on this processor
    
  result = mbcore->get_entities_by_dimension(0,celldim,AllCells,false);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get cells" << std::endl;
    assert(result == MB_SUCCESS);
  }
    
  // Get not owned cells (which is the same as ghost cells)

  result = mbcomm->get_pstatus_entities(celldim,PSTATUS_GHOST,GhostCells);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get ghost cells" << std::endl;
    assert(result == MB_SUCCESS);
  }
    
  // Subtract from all cells on processor to get owned cells only
    
  OwnedCells = AllCells;  // I think we DO want a data copy here
  OwnedCells -= GhostCells;
}



void Mesh_maps_moab::init_set_info() {
  int maxnsets, result;



  // Get material, sideset and nodeset tags

  result = mbcore->tag_get_handle(MATERIAL_SET_TAG_NAME,mattag);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get tag for material sets" << std::endl;
    assert(result == MB_SUCCESS);
  }
  mbcore->tag_get_handle(NEUMANN_SET_TAG_NAME,sstag);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get tag for side sets" << std::endl;
    assert(result == MB_SUCCESS);
  }
  mbcore->tag_get_handle(DIRICHLET_SET_TAG_NAME,nstag);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get tag for node sets" << std::endl;
    assert(result == MB_SUCCESS);
  }


  std::vector<MBTag> tag_handles;
  mbcore->tag_get_tags(tag_handles);



  nsets = 0;
  for (int i = 0; i < tag_handles.size(); i++) {
    MBTag tag = tag_handles[i];
    int n;
    
    if (tag != mattag && tag != sstag && tag != nstag) continue;

    mbcore->get_number_entities_by_type_and_tag(0,MBENTITYSET,&tag,0,1,n);
    
    nsets +=  n;
  }
    


  if (serial_run) {
    maxnsets = nsets;
  }
  else {
    // Figure out the maximum number of sets on any processor,

    MPI_Allreduce(&nsets,&maxnsets,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  }


  // We'll pad maxnsets so that we can try and avoid (but not
  // guarantee) reallocation

  maxnsets *= 2;

  setids = new int[maxnsets];
  setdims = new int[maxnsets];

  nsets = 0;
  for (int i = 0; i < tag_handles.size(); i++) {
    MBTag tag = tag_handles[i];
    MBRange tagsets;
    MBRange tagset;
    int n;
    
    if (tag != mattag && tag != sstag && tag != nstag) continue;

    mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&tag,0,1,tagsets);

    for (MBRange::iterator it = tagsets.begin(); it != tagsets.end(); ++it) {
      MBEntityHandle set = *it;
      int val;
      
      mbcore->tag_get_data(tag,&set,1,&val);
      setids[nsets] = val;

      MBRange setents;
      mbcore->get_entities_by_handle(set,setents,false);
      int dim = mbcore->dimension_from_handle(*(setents.begin()));

      setdims[nsets] = dim;

      nsets++;
    }

  }

  for (int i = nsets; i < maxnsets; i++) {
    setids[i] = 0;
    setdims[i] = -1;
  }


  if (serial_run) return;


  // In a parallel run, we want the global number of sets and the
  // global list of set ids 

  int rank = mbcomm->rank();
  int nprocs = mbcomm->size();


  // Collate all the sets across processors

  
  int *allsetids = new int[nprocs*maxnsets];
  int *allsetdims = new int[nprocs*maxnsets];
  int *allnsets = new int[nprocs];
  
  
  MPI_Allgather(setids,maxnsets,MPI_INT,allsetids,maxnsets,MPI_INT,
		MPI_COMM_WORLD);
  MPI_Allgather(setdims,maxnsets,MPI_INT,allsetdims,maxnsets,MPI_INT,
		MPI_COMM_WORLD);
  
  
  
  // Make a list of all setids/entity dimensions across all processors

  int nsetsall;
  if (rank == 0) {

    // Just add to the list of setids on processor 0 Since we doubled
    // the computed maxnsets and allocated allsetids to be of length
    // nprocs*maxnsets, we are guaranteed that we will not overwrite
    // processor 1 entries if we add unique entries from processor 1
    // to the end of processor 0's list of sets and so on

    nsetsall = nsets;

    for (int k = celldim; k >= 0; k--) {
      for (int ip = 1; ip < nprocs; ip++) {
	for (int i = 0; i < maxnsets; i++) {
	  if (allsetdims[ip*maxnsets+i] != k) continue;
	
	  int found = 0;
	  for (int j = 0; j < nsetsall; j++) {
	    if (allsetids[j] == allsetids[ip*maxnsets+i]) {
	      found = 1;
	      break;
	    }
	  }
	  
	  if (!found) {
	    allsetids[nsetsall] = allsetids[ip*maxnsets+i];
	    allsetdims[nsetsall] = allsetdims[ip*maxnsets+i];
	    nsetsall++;
	  }
	}
      }
    }

    if (nsetsall > maxnsets) {
      // We have to resize the allsetdims array

      int *tmpsetids = new int[nsetsall];
      int *tmpsetdims = new int[nsetsall];

      // is it possible to use memcpy on arrays allocated by new

      memcpy(tmpsetids, allsetids, nsetsall*sizeof(allsetids[0]));
      memcpy(tmpsetdims, allsetdims, nsetsall*sizeof(allsetdims[0]));

      delete [] allsetids;
      delete [] allsetdims;

      maxnsets = nsetsall;
      allsetids = new int[maxnsets];
      allsetdims = new int[maxnsets];

      memcpy(allsetids, tmpsetids, nsetsall*sizeof(tmpsetids[0]));
      memcpy(allsetdims, tmpsetdims, nsetsall*sizeof(tmpsetdims[0]));
    }



    for (int i = 1; i < nprocs; i++) {
      memcpy(allsetids+maxnsets*i,allsetids,nsetsall*sizeof(allsetids[0]));
      memcpy(allsetdims+maxnsets*i,allsetdims,nsetsall*sizeof(allsetdims[0]));
    }

    for (int i = 0; i < nprocs; i++)
      allnsets[i] = nsetsall;

  } // if rank == 0;



  
  MPI_Scatter(allnsets,1,MPI_INT,&nsets,1,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Scatter(allsetids,maxnsets,MPI_INT,setids,maxnsets,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Scatter(allsetdims,maxnsets,MPI_INT,setdims,maxnsets,MPI_INT,0,MPI_COMM_WORLD);
  
  delete [] allsetids;
  delete [] allsetdims;
  
}




// Number of OWNED, GHOST or USED entities of different types

    
unsigned int Mesh_maps_moab::count_entities (Mesh_data::Entity_kind kind, 
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
      return AllVerts.size();
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
      return AllFaces.size();
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
      return AllCells.size();
      break;
    default:
      return 0;
    }

    break;
  default:
    std::cerr << "Count requested for unknown entity type" << std::endl;
  }
}


void Mesh_maps_moab::cell_to_faces (unsigned int cellid, 
				    std::vector<unsigned int>::iterator begin, 
				    std::vector<unsigned int>::iterator end)
{
  cell_to_faces(cellid, &(*begin), &(*end));
}

void Mesh_maps_moab::cell_to_faces (unsigned int cellid, 
				    unsigned int *begin, 
				    unsigned int *end)
{
  MBEntityHandle cell;
  MBRange cell_faces;
  std::vector<MBEntityHandle> cell_nodes, face_nodes;
  int *cell_faceids;
  int nf;
  int cfstd[6][4] = {{0,1,5,4},    // Expected cell-face-node pattern
		     {1,2,6,5},
		     {2,3,7,6},
		     {0,4,7,3},
		     {0,3,2,1},
		     {4,5,6,7}};


  cell = cell_id_to_handle[cellid];
      
  mbcore->get_adjacencies(&cell, 1, facedim, true, cell_faces, 
			  MBInterface::INTERSECT);
  nf = cell_faces.size();
  assert ((unsigned int) (end - begin) >= nf);

  cell_faceids = new int[nf];			


  // Have to re-sort the faces according a specific template for hexes


  if (nf == 6) { // Hex

    MBEntityHandle *ordfaces, face;

    ordfaces = new MBEntityHandle[6];

    mbcore->get_connectivity(&cell, 1, cell_nodes);

    for (int i = 0; i < nf; i++) {
  
      // Search for a face that has all the expected nodes

      bool found = false;
      int j;
      for (j = 0; j < nf; j++) {

	face = cell_faces[j];
	mbcore->get_connectivity(&face, 1, face_nodes);

	// Check if this face has all the expected nodes

	bool all_present = true;

	for (int k = 0; k < 4; k++) {
	  unsigned int node = cell_nodes[cfstd[i][k]];

	  if (face_nodes[0] != node && face_nodes[1] != node &&
	      face_nodes[2] != node && face_nodes[3] != node) {
	    all_present = false;
	    break;
	  }
	}

	if (all_present) {
	  found = true;
	  break;
	}
      }

      assert(found);

      if (found)
	ordfaces[i] = face;
    }


    mbcore->tag_get_data(lid_tag,ordfaces,6,cell_faceids);

    delete [] ordfaces;

  }
  else {
    mbcore->tag_get_data(lid_tag,cell_faces,cell_faceids);
  }

  std::copy (cell_faceids, cell_faceids+nf, begin);

  delete [] cell_faceids;
}




void Mesh_maps_moab::cell_to_face_dirs (unsigned int cellid, 
				       std::vector<int>::iterator begin, 
				       std::vector<int>::iterator end){
  cell_to_face_dirs(cellid, &(*begin), &(*end));
}

void Mesh_maps_moab::cell_to_face_dirs (unsigned int cellid, 
				       int *begin, 
				       int *end)
{
  MBEntityHandle cell;
  unsigned int *cell_faces;
  int *cell_facedirs;
  int j,nf, result;
  

  nf = 6; // HEX SPECIFIC - THIS NUMBER SHOULD COME FROM cell_to_faces

  cell = cell_id_to_handle[cellid];

  cell_faces = new unsigned int[nf];

  cell_to_faces(cellid,cell_faces,cell_faces+nf);

  assert ((unsigned int) (end - begin) >= nf);
  cell_facedirs = new int[nf];

  for (int i = 0; i < nf; i++) {
    MBEntityHandle face = face_id_to_handle[cell_faces[i]];
    int sidenum, offset;

    result = mbcore->side_number(cell,face,sidenum,cell_facedirs[i],offset);

    if (result != MB_SUCCESS) {
      cerr << "Could not find face dir in cell" << std::endl;
      assert(result == MB_SUCCESS);
    }
    
    // If this is a ghost face and the master has the opposite direction
    // we are supposed to flip it

    if (faceflip[cell_faces[i]]) cell_facedirs[i] *= -1;
  }
    
  std::copy (cell_facedirs, cell_facedirs+nf, begin);

  delete [] cell_faces;
  delete [] cell_facedirs;
}



void Mesh_maps_moab::cell_to_nodes (unsigned int cellid, 
				    std::vector<unsigned int>::iterator begin, 
				    std::vector<unsigned int>::iterator end)
{
  cell_to_nodes(cellid, &(*begin), &(*end));
}



void Mesh_maps_moab::cell_to_nodes (unsigned int cellid, unsigned int *begin, 
				    unsigned int *end)
{
  MBEntityHandle cell;
  std::vector<MBEntityHandle> cell_nodes;
  int *cell_nodeids;
  int nn;


  cell = cell_id_to_handle[cellid];
      
  mbcore->get_connectivity(&cell, 1, cell_nodes);

  nn = cell_nodes.size();
  assert ((unsigned int) (end - begin) >= nn);
  cell_nodeids = new int[nn];
    
  for (int i = 0; i < nn; i++)
    mbcore->tag_get_data(lid_tag,&(cell_nodes[i]),1,&(cell_nodeids[i]));

  std::copy (cell_nodeids, cell_nodeids+nn, begin);

  delete [] cell_nodeids;
}




void Mesh_maps_moab::face_to_nodes (unsigned int faceid, 
				    std::vector<unsigned int>::iterator begin, 
				    std::vector<unsigned int>::iterator end)
{
  face_to_nodes(faceid, &(*begin), &(*end));
}

  
void Mesh_maps_moab::face_to_nodes (unsigned int faceid, unsigned int *begin, 
				    unsigned int *end) {
  MBEntityHandle face;
  std::vector<MBEntityHandle> face_nodes;
  int *face_nodeids;
  int nn;

  face = face_id_to_handle[faceid];

  mbcore->get_connectivity(&face, 1, face_nodes, true);

  nn = face_nodes.size();
  assert ((unsigned int) (end - begin) >= nn);

  face_nodeids = new int[nn];
  if (faceflip[faceid]) {
    for (int i = nn-1; i >= 0; i--) 
      mbcore->tag_get_data(lid_tag,&(face_nodes[i]),1,&(face_nodeids[nn-i-1]));
  }
  else {
    for (int i = 0; i < nn; i++) 
      mbcore->tag_get_data(lid_tag,&(face_nodes[i]),1,&(face_nodeids[i]));
  }

  std::copy (face_nodeids, face_nodeids+nn, begin);

  delete [] face_nodeids;
}
  


void Mesh_maps_moab::node_to_coordinates (unsigned int node_id, 
					  std::vector<double>::iterator begin, 
					  std::vector<double>::iterator end)
{
  node_to_coordinates(node_id, &(*begin), &(*end));
}

void Mesh_maps_moab::node_to_coordinates (unsigned int node_id, 
					  double *begin, 
					  double *end)
{
  MBEntityHandle node;

  assert ((unsigned int) (end-begin) >= spacedim);
    
  double coords[3];

  node = vtx_id_to_handle[node_id];

  mbcore->get_coords(&node, 1, coords);

  std::copy (coords, coords+spacedim, begin);

}



void Mesh_maps_moab::cell_to_coordinates (unsigned int cellid, 
					  std::vector<double>::iterator begin, 
					  std::vector<double>::iterator end)
{
  cell_to_coordinates(cellid, &(*begin), &(*end));
}


void Mesh_maps_moab::cell_to_coordinates (unsigned int cellid, 
					  double *begin, double *end)
{
  MBEntityHandle cell;
  std::vector<MBEntityHandle> cell_nodes;
  double *coords;
  int nn;


  cell = cell_id_to_handle[cellid];
      
  mbcore->get_connectivity(&cell, 1, cell_nodes);

  nn = cell_nodes.size();
  assert ((unsigned int) (end - begin) >= nn);

  coords = new double[spacedim*nn];

  for (int i = 0; i < nn; i++)
    mbcore->get_coords(&(cell_nodes[i]),1,coords+spacedim*i);

  std::copy (coords, coords+spacedim*nn, begin);

  delete [] coords;
}




void Mesh_maps_moab::face_to_coordinates (unsigned int faceid, 
					  std::vector<double>::iterator begin, 
					  std::vector<double>::iterator end)
{
  face_to_coordinates(faceid, &(*begin), &(*end));
}
  
void Mesh_maps_moab::face_to_coordinates (unsigned int faceid, 
					  double *begin, 
					  double *end)
{
    MBEntityHandle face;
    std::vector<MBEntityHandle> face_nodes;
    double *coords;
    int nn;

    face = face_id_to_handle[faceid];

    mbcore->get_connectivity(&face, 1, face_nodes, true);

    nn = face_nodes.size();
    assert ((unsigned int) (end - begin) >= nn);

    coords = new double[spacedim*nn];
    
    if (faceflip[faceid]) {
      for (int i = nn-1; i >=0; i--)
	mbcore->get_coords(&(face_nodes[i]),1,coords+spacedim*(nn-i-1));
    }
    else {
      for (int i = 0; i < nn; i++)
	mbcore->get_coords(&(face_nodes[i]),1,coords+spacedim*i);
    }

    std::copy (coords, coords+spacedim*nn, begin);

    delete [] coords;
}
  


void Mesh_maps_moab::get_set_ids (Mesh_data::Entity_kind kind, 
			      std::vector<unsigned int>::iterator begin, 
			      std::vector<unsigned int>::iterator end) const {
  get_set_ids(kind, &(*begin), &(*end));
}


void Mesh_maps_moab::get_set_ids (Mesh_data::Entity_kind kind, 
			      unsigned int *begin, unsigned int *end) const {
  int sids[nsets];
  int i, ns=0;

  switch (kind) {
  case Mesh_data::CELL: {

    for (i = 0; i < nsets; i++)
      if (setdims[i] == celldim)
	sids[ns++] = setids[i];

    assert ((unsigned int) (end - begin) >= ns);

    std::copy (sids, sids + ns, begin);

    return;
    break;
  }

  case Mesh_data::FACE: {

    for (i = 0; i < nsets; i++)
      if (setdims[i] == facedim)
	sids[ns++] = setids[i];

    assert ((unsigned int) (end - begin) >= ns);

    std::copy (sids, sids + ns, begin);

    return;
    break;
  }

  case Mesh_data::NODE: {

    for (i = 0; i < nsets; i++)
      if (setdims[i] == 0)
	sids[ns++] = setids[i];

    assert ((unsigned int) (end - begin) >= ns);

    std::copy (sids, sids + ns, begin);

    return;
    break;
  }

  default:
    return;
  }
    
}



void Mesh_maps_moab::get_set (unsigned int set_id, 
			      Mesh_data::Entity_kind kind, 
			      Element_Category category,
			      std::vector<unsigned int>::iterator begin, 
			      std::vector<unsigned int>::iterator end) const {
  get_set(set_id, kind, category, &(*begin), &(*end));
}

void Mesh_maps_moab::get_set (unsigned int set_id, 
			      Mesh_data::Entity_kind kind, 
			      Element_Category category,
			      unsigned int *begin, unsigned int *end) const {

  int i;
  const void *valarr[1] = {&set_id};

  switch (kind) {
  case Mesh_data::CELL: {
    MBRange matsets, cellrange;
    int *cells;
    int ncells;

    mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&mattag,valarr,1,matsets);

    if (matsets.size() == 0) return;

    MBEntityHandle matset = *(matsets.begin());

    mbcore->get_entities_by_dimension(matset,celldim,cellrange);
    if (cellrange.size()) {
      if (category == OWNED)
	cellrange -= GhostCells;
      else if (category == GHOST)
	cellrange -= OwnedCells;
    }

    ncells = cellrange.size();

    assert ((unsigned int) (end - begin) >= ncells);

    cells = new int[ncells];
    i = 0;
    for (MBRange::iterator jt = cellrange.begin(); jt != cellrange.end(); ++jt) {
      MBEntityHandle ent = *jt;

      mbcore->tag_get_data(lid_tag,&ent,1,&(cells[i++]));
    }

    std::copy(cells, cells + ncells, begin);

    delete [] cells;

    return;
    break;
  }
  case Mesh_data::FACE: {
    MBRange sidesets, facerange;
    int *faces;
    int nfaces;

    mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&sstag,valarr,1,sidesets);

    if (sidesets.size() == 0) return;

    MBEntityHandle sideset = *(sidesets.begin());

    mbcore->get_entities_by_dimension(sideset,facedim,facerange); 
    
    if (facerange.size()) {
      if (category == OWNED)
	facerange -= NotOwnedFaces;
      else if (category == GHOST)
	facerange -= OwnedFaces;
    }
    nfaces = facerange.size();

    assert ((unsigned int) (end - begin) >= nfaces);

    faces = new int[nfaces];
    i = 0;
    for (MBRange::iterator jt = facerange.begin(); jt != facerange.end(); ++jt) {
      MBEntityHandle ent = *jt;

      mbcore->tag_get_data(lid_tag,&ent,1,&(faces[i++]));
    }

    std::copy(faces, faces + nfaces, begin);
      
    delete [] faces;

    return;
    break;
  }
  case Mesh_data::NODE: {
    MBRange nodesets, noderange;
    int *nodes;
    int nnodes;

    mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&nstag,valarr,1,nodesets);

    if (nodesets.size() == 0) return;

    MBEntityHandle nodeset = *(nodesets.begin());

    mbcore->get_entities_by_dimension(nodeset,0,noderange);
    if (noderange.size()) {
      if (category == OWNED)
	noderange -= NotOwnedVerts;
      else if (category == GHOST)
	noderange -= OwnedVerts;
    }

    nnodes = noderange.size();

    assert ((unsigned int) (end - begin) >= nnodes);

    nodes = new int[nnodes];
    i = 0;
    for (MBRange::iterator jt = noderange.begin(); jt != noderange.end(); ++jt) {
      MBEntityHandle ent = *jt;

      mbcore->tag_get_data(lid_tag,&ent,1,&(nodes[i++]));
    }

    std::copy(nodes, nodes + nnodes, begin);
      
    delete [] nodes;

    return;
    break;
  }
  default:
    return;
  }
    
}


unsigned int Mesh_maps_moab::num_sets(Mesh_data::Entity_kind kind) const {
  int n = 0;
    
  switch (kind) {
  case Mesh_data::CELL:
    for (int i = 0; i < nsets; i++)
      if (setdims[i] == celldim) n++;
    return n;
  case Mesh_data::FACE:
    for (int i = 0; i < nsets; i++)
      if (setdims[i] == facedim) n++;
    return n;
  case Mesh_data::NODE:
    for (int i = 0; i < nsets; i++)
      if (setdims[i] == 0) n++;
    return n;
  default:
    return 0;
  }
}

bool Mesh_maps_moab::valid_set_id (unsigned int id, Mesh_data::Entity_kind kind) const {
  int n = 0;

  switch (kind) {
  case Mesh_data::CELL:
    for (int i = 0; i < nsets; i++)
      if (setdims[i] == celldim && setids[i] == id) return true;
    break;
  case Mesh_data::FACE:
    for (int i = 0; i < nsets; i++)
      if (setdims[i] == facedim && setids[i] == id) return true;
    break;
  case Mesh_data::NODE:
    for (int i = 0; i < nsets; i++)
      if (setdims[i] == 0 && setids[i] == id) return true;
    break;
  default:
    return false;
  }    

  return false;
}

unsigned int Mesh_maps_moab::get_set_size(unsigned int set_id, 
				      Mesh_data::Entity_kind kind, 
				      Element_Category category) const {   

  const int loc_setid = set_id;
  const void* valarr[1];
  valarr[0] =  &loc_setid;

  switch (kind) {
  case Mesh_data::CELL: {
    MBRange matsets, cellrange;
    int ncells;

    mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&mattag,valarr,1,matsets);

    if (matsets.size() == 0) return 0;

    MBEntityHandle matset = *(matsets.begin());

    mbcore->get_entities_by_dimension(matset,celldim,cellrange);
    if (category == OWNED)
      cellrange -= GhostCells;
    else if (category == GHOST)
      cellrange -= OwnedCells;

    return cellrange.size();
    break;
  }

  case Mesh_data::FACE: {
    MBRange sidesets, facerange;
    int nfaces;

    mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&sstag,valarr,1,sidesets);
    if (sidesets.size() == 0) return 0;

    MBEntityHandle sideset = *(sidesets.begin());

    mbcore->get_entities_by_dimension(sideset,facedim,facerange);
    if (category == OWNED)
      facerange -= NotOwnedFaces;
    else if (category == GHOST)
      facerange -= OwnedFaces;

    return facerange.size();
    break;
  }

  case Mesh_data::NODE: {
    MBRange nodesets, noderange;
    int nnodes;

    mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&nstag,valarr,1,nodesets);
    if (nodesets.size() == 0) return 0;

    MBEntityHandle nodeset = *(nodesets.begin());

    mbcore->get_entities_by_dimension(nodeset,0,noderange);
    if (category == OWNED)
      noderange -= NotOwnedVerts;
    else if (category == GHOST)
      noderange -= OwnedVerts;

    return noderange.size();
    break;
  }

  default:
    return 0;
  }
    
}


// Epetra map for cells - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_maps_moab::init_cell_map ()
{
  int *cell_gids;
  int ncell;

  if (!serial_run) {

    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells

    
    cell_gids = new int[OwnedCells.size()+GhostCells.size()];
    
    mbcore->tag_get_data(gid_tag,OwnedCells,cell_gids);
    ncell = OwnedCells.size();
    

    cell_map_wo_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);
    



    mbcore->tag_get_data(gid_tag,GhostCells,&(cell_gids[ncell]));
    
    ncell += GhostCells.size();

    cell_map_w_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);

  }
  else {
    cell_gids = new int[AllCells.size()];

    mbcore->tag_get_data(gid_tag,AllCells,cell_gids);
    ncell = AllCells.size();

    cell_map_wo_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);
  }

  delete [] cell_gids;

}




// Epetra map for faces - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_maps_moab::init_face_map ()
{
  int *face_gids;
  int nface;

  if (!serial_run) {

    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells

    
    face_gids = new int[OwnedFaces.size()+NotOwnedFaces.size()];
    
    mbcore->tag_get_data(gid_tag,OwnedFaces,face_gids);
    nface = OwnedFaces.size();
    
    face_map_wo_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*epcomm);


    mbcore->tag_get_data(gid_tag,NotOwnedFaces,&(face_gids[nface]));
    
    nface += NotOwnedFaces.size();

    face_map_w_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*epcomm);

  }
  else {
    face_gids = new int[AllFaces.size()];

    mbcore->tag_get_data(gid_tag,AllFaces,face_gids);
    nface = AllFaces.size();

    face_map_wo_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*epcomm);
  }

  delete [] face_gids;

}




// Epetra map for nodes - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_maps_moab::init_node_map ()
{
  int *vert_gids;
  int nvert;

  if (!serial_run) {

    // For parallel runs create map without and with ghost verts included
    // Also, put in owned cells before the ghost verts

    
    vert_gids = new int[OwnedVerts.size()+NotOwnedVerts.size()];
    
    mbcore->tag_get_data(gid_tag,OwnedVerts,vert_gids);
    nvert = OwnedVerts.size();
    
    node_map_wo_ghosts_ = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);
    



    mbcore->tag_get_data(gid_tag,NotOwnedVerts,&(vert_gids[nvert]));
    
    nvert += NotOwnedVerts.size();

    node_map_w_ghosts_ = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);

  }
  else {
    vert_gids = new int[AllVerts.size()];

    mbcore->tag_get_data(gid_tag,AllVerts,vert_gids);
    nvert = AllVerts.size();

    node_map_wo_ghosts_ = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);
  }

  delete [] vert_gids;

}


unsigned int Mesh_maps_moab::GID(unsigned int lid, Mesh_data::Entity_kind kind) {
  MBEntityHandle ent;
  unsigned int gid;

  switch (kind) {
  case Mesh_data::NODE:
    ent = vtx_id_to_handle[lid];
    break;

  case Mesh_data::FACE:
    ent = face_id_to_handle[lid];
    break;

  case Mesh_data::CELL:
    ent = cell_id_to_handle[lid];
    break;
  default:
    std::cerr << "Global ID requested for unknown entity type" << std::endl;
  }

  mbcore->tag_get_data(gid_tag,&ent,1,&gid);

  return gid;
}



