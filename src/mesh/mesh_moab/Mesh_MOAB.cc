#include "Mesh_MOAB.hh"
//#include <Teuchos_RCP.hpp>

using namespace std;
using namespace moab;


namespace Amanzi
{
namespace AmanziMesh
{

  // Constructor - load up mesh from file

  Mesh_MOAB::Mesh_MOAB (const char *filename, const Epetra_MpiComm *comm, 
			const AmanziGeometry::GeometricModelPtr& gm) 
{
  int result, rank;
  
  clear_internals_();

    
  // Core MOAB object
    
  mbcore = new MBCore();

  if (comm) {
    // MOAB's parallel communicator
    int mbcomm_id;
    mbcomm = new ParallelComm(mbcore,comm->GetMpiComm(),&mbcomm_id);

    if (!mbcomm) {
      cerr << "Failed to initialize MOAB communicator\n";
      assert(mbcomm == 0);
    }
  }

  Mesh::set_comm(comm);


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
			NULL,NULL,0);
      
    rank = mbcomm->rank();
      
  }
  else {

    // Load serial mesh

    result =
      mbcore->load_file(filename,NULL,NULL,NULL,NULL,0);

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
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting number of entities of dim 3" << std::endl;
    assert(result == MB_SUCCESS);
  }
  if (nent) {
    celldim = 3;
    facedim = 2;
  }
  else {
    result = mbcore->get_number_entities_by_dimension(0,2,nent,false);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting number of entities of dim 2" << std::endl;
      assert(result == MB_SUCCESS);
    }
    if (nent) {
      celldim = 2;
      facedim = 1;
    }
    else {
      std::cerr << "Flow code works only on 2D and 3D meshes" << std::endl;
      assert(nent > 0);
    }
  }
      


  // Set the geometric model that this mesh is related to

  set_geometric_model(gm);



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


Mesh_MOAB::~Mesh_MOAB() {
  delete cell_map_wo_ghosts_;
  delete cell_map_w_ghosts_;
  delete face_map_wo_ghosts_;
  delete face_map_w_ghosts_;
  delete node_map_wo_ghosts_;
  delete node_map_w_ghosts_;
  delete [] setids_;
  delete [] setdims_;
  delete [] faceflip;
  delete mbcore;
}


// Some initializations

void Mesh_MOAB::clear_internals_ () 
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
  setids_ = setdims_ = NULL;

  Mesh::set_geometric_model((Amanzi::AmanziGeometry::GeometricModelPtr) NULL);
}


void Mesh_MOAB::init_id_handle_maps() {
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
    result = mbcore->tag_set_data(lid_tag,&vtx,1,&i);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting local ID for vertex" << std::endl;
      assert(result == MB_SUCCESS);
    }
    vtx_id_to_handle[i++] = vtx;
  }    
  for (MBRange::iterator it = NotOwnedVerts.begin(); it != NotOwnedVerts.end(); ++it) {
    MBEntityHandle vtx = *it;
    result = mbcore->tag_set_data(lid_tag,&vtx,1,&i);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting local ID for vertex" << std::endl;
      assert(result == MB_SUCCESS);
    }
    vtx_id_to_handle[i++] = vtx;
  }
    


  nf = AllFaces.size();

  face_id_to_handle.reserve(nf);

  i = 0;
  for (MBRange::iterator it = OwnedFaces.begin(); it != OwnedFaces.end(); ++it) {
    MBEntityHandle face = *it;
    result = mbcore->tag_set_data(lid_tag,&face,1,&i);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting local ID for face" << std::endl;
      assert(result == MB_SUCCESS);
    }
    face_id_to_handle[i++] = face;
  }
  for (MBRange::iterator it = NotOwnedFaces.begin(); it != NotOwnedFaces.end(); ++it) {
    MBEntityHandle face = *it;
    result = mbcore->tag_set_data(lid_tag,&face,1,&i);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting local ID for face" << std::endl;
      assert(result == MB_SUCCESS);
    }
    face_id_to_handle[i++] = face;
  }
    


  nc = AllCells.size();

  cell_id_to_handle.reserve(nc);

  i = 0;
  for (MBRange::iterator it = OwnedCells.begin(); it != OwnedCells.end(); ++it) {
    MBEntityHandle cell = *it;
    result = mbcore->tag_set_data(lid_tag,&cell,1,&i);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting local ID for cell" << std::endl;
      assert(result == MB_SUCCESS);
    }
    cell_id_to_handle[i++] = cell;
  }
  for (MBRange::iterator it = GhostCells.begin(); it != GhostCells.end(); ++it) {
    MBEntityHandle cell = *it;
    result = mbcore->tag_set_data(lid_tag,&cell,1,&i);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting local ID for cell" << std::endl;
      assert(result == MB_SUCCESS);
    }
    cell_id_to_handle[i++] = cell;
  }
    

}



void Mesh_MOAB::init_global_ids() {
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
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag handle for GLOBAL_ID" << std::endl;
      assert(result == MB_SUCCESS);
    }
      
    int nent = AllVerts.size();
    int *gids = new int[nent];
    for (int i = 0; i < nent; i++) gids[i] = i;

    result = mbcore->tag_set_data(gid_tag,AllVerts,gids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem setting global IDs for vertices" << std::endl;
      assert(result == MB_SUCCESS);
    }

    delete [] gids;


    nent = AllFaces.size();
    gids = new int[nent];
    for (int i = 0; i < nent; i++) gids[i] = i;

    result = mbcore->tag_set_data(gid_tag,AllFaces,gids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem setting global IDs for faces" << std::endl;
      assert(result == MB_SUCCESS);
    }

    delete [] gids;


    nent = AllCells.size();
    gids = new int[nent];
    for (int i = 0; i < nent; i++) gids[i] = i;

    result = mbcore->tag_set_data(gid_tag,AllCells,gids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem setting global IDs for cells" << std::endl;
      assert(result == MB_SUCCESS);
    }

    delete [] gids;
  }


}


void Mesh_MOAB::init_pvert_lists() {
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


void Mesh_MOAB::init_pface_lists() {
  int result;


  // Make MOAB create the missing 'faces' (faces in 3D, edges in
  // 2D). We do this by looping over the cells and asking for their
  // faces with the create_if_missing=true option


  for (MBRange::iterator it = AllCells.begin(); it != AllCells.end(); it++) {
    MBEntityHandle cell = *it;
    MBRange cfaces;
      
    result = mbcore->get_adjacencies(&cell,1,facedim,true,cfaces,MBCore::UNION);
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


void Mesh_MOAB::init_pface_dirs() {
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
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting new tag handle" << std::endl;
    assert(result == MB_SUCCESS);
  }
  result = mbcore->tag_create("TMP_FC1_TAG",sizeof(int),
			      MB_TAG_DENSE, MB_TYPE_INTEGER, tmp_fc1_tag,
			      &zero,true);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting new tag handle" << std::endl;
    assert(result == MB_SUCCESS);
  }

  
  for (MBRange::iterator it = OwnedFaces.begin(); it != OwnedFaces.end(); it++) {
    MBRange fcells;
    face = *it;

    result = mbcore->get_adjacencies(&face,1,celldim,false,fcells,MBCore::UNION);
    if (result != MB_SUCCESS) {
      cout << "Could not get cells of face" << std::endl;
      assert(result == MB_SUCCESS);
    }

    result = mbcore->tag_set_data(tmp_fc0_tag,&face,1,&zero);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem setting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
    result = mbcore->tag_set_data(tmp_fc1_tag,&face,1,&zero);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem setting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }


    for (MBRange::iterator jt = fcells.begin(); jt != fcells.end(); ++jt) {

      cell = *jt;
      result = mbcore->side_number(cell,face,sidenum,facedir,offset);
      if (result != MB_SUCCESS) {
	cout << "Could not get face dir w.r.t. cell" << std::endl;
	assert(result == MB_SUCCESS);
      }

      result = mbcore->tag_get_data(gid_tag,&cell,1,&cell_gid);
      if (result != MB_SUCCESS) {
	std::cerr << "Problem getting tag data" << std::endl;
	assert(result == MB_SUCCESS);
      }
      cell_gid += 1; 

      if (facedir == 1) {
	result = mbcore->tag_set_data(tmp_fc0_tag,&face,1,&cell_gid);
	if (result != MB_SUCCESS) {
	  std::cerr << "Problem setting tag data" << std::endl;
	  assert(result == MB_SUCCESS);
	}
      }
      else {
	result = mbcore->tag_set_data(tmp_fc1_tag,&face,1,&cell_gid);
	if (result != MB_SUCCESS) {
	  std::cerr << "Problem setting tag data" << std::endl;
	  assert(result == MB_SUCCESS);
	}
      }
    }

  }



  result = mbcomm->exchange_tags(tmp_fc0_tag,AllFaces);
  if (result != MB_SUCCESS) {
    cout << "Could not get exchange tag data successfully" << std::endl;
    assert(result == MB_SUCCESS);
  }

  result = mbcomm->exchange_tags(tmp_fc1_tag,AllFaces);
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


    result = mbcore->tag_get_data(tmp_fc0_tag,&face,1,&master_cell0_gid);
    if (result != MB_SUCCESS) {
      cout << "Could not get face tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
    result = mbcore->tag_get_data(tmp_fc1_tag,&face,1,&master_cell1_gid);
    if (result != MB_SUCCESS) {
      cout << "Could not get face tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }


    result = mbcore->get_adjacencies(&face,1,celldim,false,fcells,MBCore::UNION);
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
	if (result != MB_SUCCESS) {
	  std::cerr << "Problem getting tag data" << std::endl;
	  assert(result == MB_SUCCESS);
	}
	ghost_cell0_gid += 1;
      }
      else {
	result = mbcore->tag_get_data(gid_tag,&cell,1,&ghost_cell1_gid);
	if (result != MB_SUCCESS) {
	  std::cerr << "Problem getting tag data" << std::endl;
	  assert(result == MB_SUCCESS);
	}
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
	result = mbcore->tag_get_data(gid_tag,&face,1,&face_gid);
	if (result != MB_SUCCESS) {
	  std::cerr << "Problem getting tag data" << std::endl;
	  assert(result == MB_SUCCESS);
	}
	//

	cout << "Face cells mismatch between master and ghost (processor " << mbcomm->rank() << ")" << std::endl;
	cout << " Face " << face_gid << std::endl;
	cout << "Master cells " << master_cell0_gid << " " << master_cell1_gid << std::endl;
	cout << "Ghost cells " << ghost_cell0_gid << " " << ghost_cell1_gid << std::endl;
      }
    }
  }

}


void Mesh_MOAB::init_pcell_lists() {
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



void Mesh_MOAB::init_set_info() {
  int maxnsets, result;



  // Get material, sideset and nodeset tags

  result = mbcore->tag_get_handle(MATERIAL_SET_TAG_NAME,mattag);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get tag for material sets" << std::endl;
    assert(result == MB_SUCCESS);
  }
  result = mbcore->tag_get_handle(NEUMANN_SET_TAG_NAME,sstag);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get tag for side sets" << std::endl;
    assert(result == MB_SUCCESS);
  }
  result = mbcore->tag_get_handle(DIRICHLET_SET_TAG_NAME,nstag);
  if (result != MB_SUCCESS) {
    std::cerr << "Could not get tag for node sets" << std::endl;
    assert(result == MB_SUCCESS);
  }


  std::vector<MBTag> tag_handles;
  result = mbcore->tag_get_tags(tag_handles);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting all tags" << std::endl;
    assert(result == MB_SUCCESS);
  }



  nsets = 0;
  for (int i = 0; i < tag_handles.size(); i++) {
    MBTag tag = tag_handles[i];
    int n;
    
    if (tag != mattag && tag != sstag && tag != nstag) continue;

    result = mbcore->get_number_entities_by_type_and_tag(0,MBENTITYSET,&tag,0,1,n);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag sets" << std::endl;
      assert(result == MB_SUCCESS);
    }
    
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

  setids_ = new int[maxnsets];
  setdims_ = new int[maxnsets];

  nsets = 0;
  for (int i = 0; i < tag_handles.size(); i++) {
    MBTag tag = tag_handles[i];
    MBRange tagsets;
    MBRange tagset;
    int n;
    
    if (tag != mattag && tag != sstag && tag != nstag) continue;

    result = mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&tag,0,1,tagsets);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entities by type and tag" << std::endl;
      assert(result == MB_SUCCESS);
    }

    for (MBRange::iterator it = tagsets.begin(); it != tagsets.end(); ++it) {
      MBEntityHandle set = *it;
      int val;
      
      result = mbcore->tag_get_data(tag,&set,1,&val);
      setids_[nsets] = val;

      MBRange setents;
      result = mbcore->get_entities_by_handle(set,setents,false);
      if (result != MB_SUCCESS) {
	std::cerr << "Problem getting entities by handle" << std::endl;
	assert(result == MB_SUCCESS);
      }
      int dim = mbcore->dimension_from_handle(*(setents.begin()));

      setdims_[nsets] = dim;

      nsets++;
    }

  }

  for (int i = nsets; i < maxnsets; i++) {
    setids_[i] = 0;
    setdims_[i] = -1;
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
  
  
  MPI_Allgather(setids_,maxnsets,MPI_INT,allsetids,maxnsets,MPI_INT,
		MPI_COMM_WORLD);
  MPI_Allgather(setdims_,maxnsets,MPI_INT,allsetdims,maxnsets,MPI_INT,
		MPI_COMM_WORLD);
  
  
  
  // Make a list of all setids/entity dimensions across all processors

  if (rank == 0) {

    // Just add to the list of setids on processor 0. Since we doubled
    // the computed maxnsets and allocated allsetids to be of length
    // nprocs*maxnsets, we are guaranteed that we will not overwrite
    // processor 1 entries if we add unique entries from processor 1
    // to the end of processor 0's list of sets and so on

    for (int k = celldim; k >= 0; k--) {
      for (int ip = 1; ip < nprocs; ip++) {
	for (int i = 0; i < maxnsets; i++) {
	  if (allsetdims[ip*maxnsets+i] != k) continue;
	
	  int found = 0;
	  for (int j = 0; j < nsets; j++) {
	    if (allsetids[j] == allsetids[ip*maxnsets+i]) {
	      found = 1;
	      break;
	    }
	  }
	  
	  if (!found) {
	    allsetids[nsets] = allsetids[ip*maxnsets+i];
	    allsetdims[nsets] = allsetdims[ip*maxnsets+i];
	    nsets++;
	  }
	}
      }
    }

    if (nsets > maxnsets) {
      // We have to resize the setids and setdims array

      delete [] setids_;
      delete [] setdims_;

      maxnsets = nsets;
      setids_ = new int[maxnsets];
      setdims_ = new int[maxnsets];
    }
    
    memcpy(setids_, allsetids, nsets*sizeof(allsetids[0]));
    memcpy(setdims_, allsetdims, nsets*sizeof(allsetdims[0]));

  } // if rank == 0


  // Tell all the other processors about number of sets to expect

  MPI_Bcast(&nsets,1,MPI_INT,0,MPI_COMM_WORLD);

  
  // Make sure the other processors have enough space allocated to receive
  // the expanded list of sets

  if (rank != 0) {					       
    if (nsets > maxnsets) { 
      delete [] setids_;
      delete [] setdims_;

      maxnsets = nsets;
      setids_ = new int[maxnsets];
      setdims_ = new int[maxnsets];
    }    
  }


  // Tell all the other processors about the complete list of sets
  
  MPI_Bcast(setids_,nsets,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(setdims_,nsets,MPI_INT,0,MPI_COMM_WORLD);
  

  delete [] allnsets;
  delete [] allsetids;
  delete [] allsetdims;
  
}




// Number of OWNED, GHOST or USED entities of different types

    
unsigned int Mesh_MOAB::num_entities (Entity_kind kind, 
					     Parallel_type ptype) const
{
  const int rank = (int) kind;
  const int index = ((int) ptype) - 1;

  switch (kind) {
  case NODE:

    switch (ptype) {
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


  case FACE:
    switch (ptype) {
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


  case CELL:

    switch (ptype) {
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


void Mesh_MOAB::cell_get_faces (Entity_ID cellid, Entity_ID_List *faceids) const
{
  MBEntityHandle cell;
  MBRange cell_faces;
  std::vector<MBEntityHandle> cell_nodes, face_nodes;
  int *cell_faceids;
  int nf, result;
  int cfstd[6][4] = {{0,1,5,4},    // Expected cell-face-node pattern
		     {1,2,6,5},
		     {2,3,7,6},
		     {0,4,7,3},
		     {0,3,2,1},
		     {4,5,6,7}};

  cell = cell_id_to_handle[cellid];

  faceids->clear();

      
  result = mbcore->get_adjacencies(&cell, 1, facedim, true, cell_faces, 
			  MBInterface::INTERSECT);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting faces of cell" << std::endl;
    assert(result == MB_SUCCESS);
  }
  nf = cell_faces.size();

  cell_faceids = new int[nf];			


  // Have to re-sort the faces according a specific template for hexes


  if (nf == 6) { // Hex

    MBEntityHandle *ordfaces, face;

    ordfaces = new MBEntityHandle[6];

    result = mbcore->get_connectivity(&cell, 1, cell_nodes);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting nodes of cell" << std::endl;
      assert(result == MB_SUCCESS);
    }

    for (int i = 0; i < nf; i++) {
  
      // Search for a face that has all the expected nodes

      bool found = false;
      int j;
      for (j = 0; j < nf; j++) {

	face = cell_faces[j];
	result = mbcore->get_connectivity(&face, 1, face_nodes);
	if (result != MB_SUCCESS) {
	  std::cerr << "Problem getting nodes of face" << std::endl;
	  assert(result == MB_SUCCESS);
	}


	// Check if this face has all the expected nodes

	bool all_present = true;

	for (int k = 0; k < 4; k++) {
	  Entity_ID node = cell_nodes[cfstd[i][k]];

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


    result = mbcore->tag_get_data(lid_tag,ordfaces,6,cell_faceids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
    
    delete [] ordfaces;

  }
  else {
    result = mbcore->tag_get_data(lid_tag,cell_faces,cell_faceids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
  }

  for (int i = 0; i < nf; i++)
    faceids->push_back(cell_faceids[i]);

  delete [] cell_faceids;
}




void Mesh_MOAB::cell_get_face_dirs (Entity_ID cellid, std::vector<int> *facedirs) const
{
  MBEntityHandle cell;
  Entity_ID *cell_faces;
  int *cell_facedirs;
  int j,nf, result;
  

  nf = 6; // HEX SPECIFIC - THIS NUMBER SHOULD COME FROM cell_to_faces

  cell = cell_id_to_handle[cellid];

  facedirs->clear();

  cell_faces = new Entity_ID[nf];

  cell_to_faces(cellid,cell_faces,cell_faces+nf);

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

  for (int i = 0; i < nf; i++)
    facedirs->push_back(cell_facedirs[i]);

  delete [] cell_faces;
  delete [] cell_facedirs;
}


  // Get faces of a cell and directions in which the cell uses the face 

  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order

  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise

 
void Mesh_MOAB::cell_get_faces_and_dirs (const Entity_ID cellid,
                                         Entity_ID_List *faceids,
                                         std::vector<int> *face_dirs,
					 const bool ordered) const
{
  
  MBEntityHandle cell;
  MBRange cell_faces;
  std::vector<MBEntityHandle> cell_nodes, face_nodes;
  int *cell_faceids, *cell_facedirs;
  int nf, result;
  int cfstd[6][4] = {{0,1,5,4},    // Expected cell-face-node pattern
		     {1,2,6,5},
		     {2,3,7,6},
		     {0,4,7,3},
		     {0,3,2,1},
		     {4,5,6,7}};

  cell = cell_id_to_handle[cellid];

  faceids->clear();
  face_dirs->clear();

      
  result = mbcore->get_adjacencies(&cell, 1, facedim, true, cell_faces, 
			  MBInterface::INTERSECT);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting faces of cell" << std::endl;
    assert(result == MB_SUCCESS);
  }
  nf = cell_faces.size();

  cell_faceids = new int[nf];			
  cell_facedirs = new int[nf];


  // Have to re-sort the faces according a specific template for hexes


  if (ordered && nf == 6) { // Hex

    MBEntityHandle *ordfaces, face;

    ordfaces = new MBEntityHandle[6];

    result = mbcore->get_connectivity(&cell, 1, cell_nodes);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting nodes of cell" << std::endl;
      assert(result == MB_SUCCESS);
    }

    for (int i = 0; i < nf; i++) {
  
      // Search for a face that has all the expected nodes

      bool found = false;
      int j;
      for (j = 0; j < nf; j++) {

	face = cell_faces[j];
	result = mbcore->get_connectivity(&face, 1, face_nodes);
	if (result != MB_SUCCESS) {
	  std::cerr << "Problem getting nodes of face" << std::endl;
	  assert(result == MB_SUCCESS);
	}


	// Check if this face has all the expected nodes

	bool all_present = true;

	for (int k = 0; k < 4; k++) {
	  Entity_ID node = cell_nodes[cfstd[i][k]];

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


    result = mbcore->tag_get_data(lid_tag,ordfaces,6,cell_faceids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }

    for (int i = 0; i < nf; i++) {
      MBEntityHandle face = ordfaces[i];
      int sidenum, offset;

      result = mbcore->side_number(cell,face,sidenum,cell_facedirs[i],offset);
      if (result != MB_SUCCESS) {
        cerr << "Could not find face dir in cell" << std::endl;
        assert(result == MB_SUCCESS);
      }
    
      // If this is a ghost face and the master has the opposite direction
      // we are supposed to flip it

      if (faceflip[cell_faceids[i]]) cell_facedirs[i] *= -1;
    }

    delete [] ordfaces;

  }
  else {
    result = mbcore->tag_get_data(lid_tag,cell_faces,cell_faceids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }

    for (int i = 0; i < nf; i++) {
      MBEntityHandle face = cell_faces[i];
      int sidenum, offset;

      result = mbcore->side_number(cell,face,sidenum,cell_facedirs[i],offset);
      if (result != MB_SUCCESS) {
        cerr << "Could not find face dir in cell" << std::endl;
        assert(result == MB_SUCCESS);
      }
    
      // If this is a ghost face and the master has the opposite direction
      // we are supposed to flip it

      if (faceflip[cell_faceids[i]]) cell_facedirs[i] *= -1;
    }
  }

  for (int i = 0; i < nf; i++) {
    faceids->push_back(cell_faceids[i]);
    face_dirs->push_back(cell_facedirs[i]);
  }

  delete [] cell_faceids;
  delete [] cell_facedirs;
}


void Mesh_MOAB::cell_get_nodes (Entity_ID cellid, Entity_ID_List *cnodes) const
{
  MBEntityHandle cell;
  std::vector<MBEntityHandle> cell_nodes;
  int *cell_nodeids;
  int nn, result;


  cell = cell_id_to_handle[cellid];
      
  result = mbcore->get_connectivity(&cell, 1, cell_nodes);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting nodes of cell" << std::endl;
    assert(result == MB_SUCCESS);
  }

  nn = cell_nodes.size();
  cell_nodeids = new int[nn];

  cnodes->clear();
    
  for (int i = 0; i < nn; i++) {
    result = mbcore->tag_get_data(lid_tag,&(cell_nodes[i]),1,&(cell_nodeids[i]));
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
  }

  for (int i = 0; i < nn; i++)
    cnodes->push_back(cell_nodeids[i]);

  delete [] cell_nodeids;
}




void Mesh_MOAB::face_get_nodes (Entity_ID faceid, Entity_ID_List *fnodes) const 
{
  MBEntityHandle face;
  std::vector<MBEntityHandle> face_nodes;
  int *face_nodeids;
  int nn, result;

  face = face_id_to_handle[faceid];

  fnodes->clear();

  result = mbcore->get_connectivity(&face, 1, face_nodes, true);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting nodes of face" << std::endl;
    assert(result == MB_SUCCESS);
  }

  nn = face_nodes.size();

  face_nodeids = new int[nn];
  if (faceflip[faceid]) {
    for (int i = nn-1; i >= 0; i--) {
      result = mbcore->tag_get_data(lid_tag,&(face_nodes[i]),1,&(face_nodeids[nn-i-1]));
      if (result != MB_SUCCESS) {
	std::cerr << "Problem getting tag data" << std::endl;
	assert(result == MB_SUCCESS);
      }
    }
  }
  else {
    for (int i = 0; i < nn; i++) {
      result = mbcore->tag_get_data(lid_tag,&(face_nodes[i]),1,&(face_nodeids[i]));
      if (result != MB_SUCCESS) {
	std::cerr << "Problem getting tag data" << std::endl;
	assert(result == MB_SUCCESS);
      }
    }
  }

  for (int i = 0; i < nn; i++)
    fnodes->push_back(face_nodeids[i]);

  delete [] face_nodeids;
}
  


void Mesh_MOAB::node_get_coordinates (Entity_ID node_id, AmanziGeometry::Point *ncoord) const 
{
  MBEntityHandle node;
  double coords[3];

  node = vtx_id_to_handle[node_id];

  int result = mbcore->get_coords(&node, 1, coords);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting node coordinates" << std::endl;
    assert(result == MB_SUCCESS);
  }

  ncoord->init(spacedim);
  ncoord->set(coords);

}

// Modify a node's coordinates

void Mesh_MOAB::node_set_coordinates(const AmanziMesh::Entity_ID nodeid, 
                                    const double *coords) {
  MBEntityHandle v = vtx_id_to_handle[nodeid];

  int result = mbcore->set_coords(&v, 1, coords);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem setting node coordinates" << std::endl;
    assert(result == MB_SUCCESS);
  }

}

void Mesh_MOAB::node_set_coordinates(const AmanziMesh::Entity_ID nodeid, 
                                     const AmanziGeometry::Point coords) {
  MBEntityHandle v = vtx_id_to_handle[nodeid];

  double coordarray[3] = {0.0,0.0,0.0};

  for (int i = 0; i < spacedim; i++)
    coordarray[i] = coords[i];

  int result = mbcore->set_coords(&v, 1, coordarray);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem setting node coordinates" << std::endl;
    assert(result == MB_SUCCESS);
  }

}





void Mesh_MOAB::cell_get_coordinates (Entity_ID cellid, std::vector<AmanziGeometry::Point> *ccoords) const
{
  MBEntityHandle cell;
  std::vector<MBEntityHandle> cell_nodes;
  double *coords;
  int nn, result;


  cell = cell_id_to_handle[cellid];

  ccoords->clear();
      
  result = mbcore->get_connectivity(&cell, 1, cell_nodes);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting nodes of a cell" << std::endl;
    assert(result == MB_SUCCESS);
  }

  nn = cell_nodes.size();

  coords = new double[spacedim];

  for (int i = 0; i < nn; i++) {
    result = mbcore->get_coords(&(cell_nodes[i]),1,coords);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting coordinates of a node" << std::endl;
      assert(result == MB_SUCCESS);
    }

    AmanziGeometry::Point xyz(spacedim);
    xyz.set(coords);
    ccoords->push_back(xyz);
  }

  delete [] coords;
}




void Mesh_MOAB::face_get_coordinates (Entity_ID faceid, std::vector<AmanziGeometry::Point> *fcoords) const
{
    MBEntityHandle face;
    std::vector<MBEntityHandle> face_nodes;
    double *coords;
    int nn, result;

    face = face_id_to_handle[faceid];

    fcoords->clear();

    result = mbcore->get_connectivity(&face, 1, face_nodes, true);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting nodes of face" << std::endl;
      assert(result == MB_SUCCESS);
    }

    nn = face_nodes.size();

    coords = new double[spacedim];
    
    if (faceflip[faceid]) {
      for (int i = nn-1; i >=0; i--) {
	result = mbcore->get_coords(&(face_nodes[i]),1,coords);
	if (result != MB_SUCCESS) {
	  std::cerr << "Problem getting coordinates of node" << std::endl;
	  assert(result == MB_SUCCESS);
	}

	AmanziGeometry::Point xyz(spacedim);
	xyz.set(coords);
	fcoords->push_back(xyz);
      }
    }
    else {
      for (int i = 0; i < nn; i++) {
	result = mbcore->get_coords(&(face_nodes[i]),1,coords);
	if (result != MB_SUCCESS) {
	  std::cerr << "Problem getting tag data" << std::endl;
	  assert(result == MB_SUCCESS);
	}

	AmanziGeometry::Point xyz(spacedim);
	xyz.set(coords);
	fcoords->push_back(xyz);
      }
    }

    delete [] coords;
}
  

void Mesh_MOAB::get_set_entities (const Set_Name set_name, 
                                  const Entity_kind kind, 
                                  const Parallel_type ptype,
                                  Entity_ID_List *setents) const {
}

void Mesh_MOAB::get_set_entities (const char *set_name, 
                                  const Entity_kind kind, 
                                  const Parallel_type ptype,
                                  Entity_ID_List *setents) const {
}


void Mesh_MOAB::get_set_entities (const Set_ID set_id, 
                                  const Entity_kind kind, 
                                  const Parallel_type ptype,
                                  Entity_ID_List *setents) const {

  int result;
  const void *valarr[1] = {&set_id};

  setents->clear();

  switch (kind) {
  case CELL: {
    MBRange matsets, cellrange;

    result = mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&mattag,valarr,1,matsets);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entity sets" << std::endl;
      assert(result == MB_SUCCESS);
    }


    if (matsets.size() == 0) return;

    MBEntityHandle matset = *(matsets.begin());

    result = mbcore->get_entities_by_dimension(matset,celldim,cellrange);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entities from entity set" << std::endl;
      assert(result == MB_SUCCESS);
    }

    if (cellrange.size()) {
      if (ptype == OWNED)
	cellrange -= GhostCells;
      else if (ptype == GHOST)
	cellrange -= OwnedCells;
    }


    for (MBRange::iterator jt = cellrange.begin(); jt != cellrange.end(); ++jt) {
      MBEntityHandle ent = *jt;
      int cellid;

      result = mbcore->tag_get_data(lid_tag,&ent,1,&(cellid));
      if (result != MB_SUCCESS) {
	std::cerr << "Problem getting tag data" << std::endl;
	assert(result == MB_SUCCESS);
      }
      
      setents->push_back(cellid);
    }

    return;
    break;
  }
  case FACE: {
    MBRange sidesets, facerange;

    result = mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&sstag,valarr,1,sidesets);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entity sets" << std::endl;
      assert(result == MB_SUCCESS);
    }

    if (sidesets.size() == 0) return;

    MBEntityHandle sideset = *(sidesets.begin());

    result = mbcore->get_entities_by_dimension(sideset,facedim,facerange); 
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entities by dimension" << std::endl;
      assert(result == MB_SUCCESS);
    }
    
    if (facerange.size()) {
      if (ptype == OWNED)
	facerange -= NotOwnedFaces;
      else if (ptype == GHOST)
	facerange -= OwnedFaces;
    }

    for (MBRange::iterator jt = facerange.begin(); jt != facerange.end(); ++jt) {
      MBEntityHandle ent = *jt;
      int face;

      result = mbcore->tag_get_data(lid_tag,&ent,1,&(face));
      if (result != MB_SUCCESS) {
	std::cerr << "Problem getting tag data" << std::endl;
	assert(result == MB_SUCCESS);
      }

      setents->push_back(face);
    }

    return;
    break;
  }
  case NODE: {
    MBRange nodesets, noderange;

    result = mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&nstag,valarr,1,nodesets);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem gettin entity set info" << std::endl;
      assert(result == MB_SUCCESS);
    }

    if (nodesets.size() == 0) return;

    MBEntityHandle nodeset = *(nodesets.begin());

    result = mbcore->get_entities_by_dimension(nodeset,0,noderange);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entities by dimension" << std::endl;
      assert(result == MB_SUCCESS);
    }

    if (noderange.size()) {
      if (ptype == OWNED)
	noderange -= NotOwnedVerts;
      else if (ptype == GHOST)
	noderange -= OwnedVerts;
    }

    for (MBRange::iterator jt = noderange.begin(); jt != noderange.end(); ++jt) {
      MBEntityHandle ent = *jt;
      int node;

      result = mbcore->tag_get_data(lid_tag,&ent,1,&(node));
      if (result != MB_SUCCESS) {
	std::cerr << "Problem getting tag data" << std::endl;
	assert(result == MB_SUCCESS);
      }

      setents->push_back(node);
    }

    return;
    break;
  }
  default:
    return;
  }
    
}

unsigned int Mesh_MOAB::get_set_size(const Set_Name set_name, 
                                     const Entity_kind kind, 
                                     const Parallel_type ptype) const {   
}
  
unsigned int Mesh_MOAB::get_set_size(const char *set_name, 
                                     const Entity_kind kind, 
                                     const Parallel_type ptype) const {   
}  

unsigned int Mesh_MOAB::get_set_size(const Set_ID set_id, 
                                     const Entity_kind kind, 
                                     const Parallel_type ptype) const {   
  
  const int loc_setid = set_id;
  const void* valarr[1];
  valarr[0] =  &loc_setid;
  int result;

  switch (kind) {
  case CELL: {
    MBRange matsets, cellrange;
    int ncells;

    result = mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&mattag,valarr,1,matsets);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting an entity set" << std::endl;
      assert(result == MB_SUCCESS);
    }

    if (matsets.size() == 0) return 0;

    MBEntityHandle matset = *(matsets.begin());

    result = mbcore->get_entities_by_dimension(matset,celldim,cellrange);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entities by dimension" << std::endl;
      assert(result == MB_SUCCESS);
	}

    if (ptype == OWNED)
      cellrange -= GhostCells;
    else if (ptype == GHOST)
      cellrange -= OwnedCells;

    return cellrange.size();
    break;
  }

  case FACE: {
    MBRange sidesets, facerange;
    int nfaces;

    result = mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&sstag,valarr,1,sidesets);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entity sets" << std::endl;
      assert(result == MB_SUCCESS);
    }

    if (sidesets.size() == 0) return 0;

    MBEntityHandle sideset = *(sidesets.begin());

    result = mbcore->get_entities_by_dimension(sideset,facedim,facerange);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entities by dimension" << std::endl;
      assert(result == MB_SUCCESS);
    }

    if (ptype == OWNED)
      facerange -= NotOwnedFaces;
    else if (ptype == GHOST)
      facerange -= OwnedFaces;

    return facerange.size();
    break;
  }

  case NODE: {
    MBRange nodesets, noderange;
    int nnodes;

    result = mbcore->get_entities_by_type_and_tag(0,MBENTITYSET,&nstag,valarr,1,nodesets);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting entity sets" << std::endl;
      assert(result == MB_SUCCESS);
    }

    if (nodesets.size() == 0) return 0;

    MBEntityHandle nodeset = *(nodesets.begin());

    result = mbcore->get_entities_by_dimension(nodeset,0,noderange);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }

    if (ptype == OWNED)
      noderange -= NotOwnedVerts;
    else if (ptype == GHOST)
      noderange -= OwnedVerts;

    return noderange.size();
    break;
  }

  default:
    return 0;
  }
    
}


// Upward adjacencies
//-------------------

// Cells of type 'ptype' connected to a node

void Mesh_MOAB::node_get_cells (const Entity_ID nodeid, 
				const Parallel_type ptype,
				Entity_ID_List *cellids) const 
{
  throw std::exception();
}
    
// Faces of type 'ptype' connected to a node

void Mesh_MOAB::node_get_faces (const Entity_ID nodeid, 
				const Parallel_type ptype,
				Entity_ID_List *faceids) const
{
  throw std::exception();
}
    
// Get faces of ptype of a particular cell that are connected to the
// given node

void Mesh_MOAB::node_get_cell_faces (const Entity_ID nodeid, 
				     const Entity_ID cellid,
				     const Parallel_type ptype,
				     Entity_ID_List *faceids) const
{
  throw std::exception();
}
    
// Cells connected to a face

void Mesh_MOAB::face_get_cells (const Entity_ID faceid, 
				const Parallel_type ptype,
				Entity_ID_List *cellids) const
{
  throw std::exception();
}
    


// Same level adjacencies
//-----------------------

// Face connected neighboring cells of given cell of a particular ptype
// (e.g. a hex has 6 face neighbors)

// The order in which the cellids are returned cannot be
// guaranteed in general except when ptype = USED, in which case
// the cellids will correcpond to cells across the respective
// faces given by cell_get_faces

void Mesh_MOAB::cell_get_face_adj_cells(const Entity_ID cellid,
					const Parallel_type ptype,
					Entity_ID_List *fadj_cellids) const
{
  throw std::exception();
}

// Node connected neighboring cells of given cell
// (a hex in a structured mesh has 26 node connected neighbors)
// The cells are returned in no particular order

void Mesh_MOAB::cell_get_node_adj_cells(const Entity_ID cellid,
					const Parallel_type ptype,
					Entity_ID_List *nadj_cellids) const
{
  throw std::exception();
}


    
//
// Mesh Topology for viz  
//----------------------
//
// We need a special function because certain types of degenerate
// hexes will not be recognized as any standard element type (hex,
// pyramid, prism or tet). The original topology of this element 
// without any collapsed nodes will be returned by this call.


// Original cell type 

Cell_type Mesh_MOAB::cell_get_type_4viz(const Entity_ID cellid) const 
{
  return HEX;
}
    
    
// See cell_get_nodes for details on node ordering

void Mesh_MOAB::cell_get_nodes_4viz (const Entity_ID cellid, 
				     Entity_ID_List *nodeids) const
{
  cell_get_nodes(cellid, nodeids);
}
    
    



// Epetra map for cells - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_MOAB::init_cell_map ()
{
  int *cell_gids;
  int ncell, result;
  const Epetra_Comm *epcomm = Mesh::get_comm();

  if (!serial_run) {

    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells

    
    cell_gids = new int[OwnedCells.size()+GhostCells.size()];
    
    result = mbcore->tag_get_data(gid_tag,OwnedCells,cell_gids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
    ncell = OwnedCells.size();
    

    cell_map_wo_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);
    



    result = mbcore->tag_get_data(gid_tag,GhostCells,&(cell_gids[ncell]));
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
    
    ncell += GhostCells.size();

    cell_map_w_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);

  }
  else {
    cell_gids = new int[AllCells.size()];

    result = mbcore->tag_get_data(gid_tag,AllCells,cell_gids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }

    ncell = AllCells.size();

    cell_map_wo_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);
  }

  delete [] cell_gids;

}




// Epetra map for faces - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_MOAB::init_face_map ()
{
  int *face_gids;
  int nface, result;
  const Epetra_Comm *epcomm = Mesh::get_comm();

  if (!serial_run) {

    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells

    
    face_gids = new int[OwnedFaces.size()+NotOwnedFaces.size()];
    
    result = mbcore->tag_get_data(gid_tag,OwnedFaces,face_gids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }

    nface = OwnedFaces.size();
    
    face_map_wo_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*epcomm);


    result = mbcore->tag_get_data(gid_tag,NotOwnedFaces,&(face_gids[nface]));
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
    
    nface += NotOwnedFaces.size();

    face_map_w_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*epcomm);

  }
  else {
    face_gids = new int[AllFaces.size()];

    result = mbcore->tag_get_data(gid_tag,AllFaces,face_gids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }

    nface = AllFaces.size();

    face_map_wo_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*epcomm);
  }

  delete [] face_gids;

}




// Epetra map for nodes - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_MOAB::init_node_map ()
{
  int *vert_gids;
  int nvert, result;
  const Epetra_Comm *epcomm = Mesh::get_comm();

  if (!serial_run) {

    // For parallel runs create map without and with ghost verts included
    // Also, put in owned cells before the ghost verts

    
    vert_gids = new int[OwnedVerts.size()+NotOwnedVerts.size()];
    
    result = mbcore->tag_get_data(gid_tag,OwnedVerts,vert_gids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }

    nvert = OwnedVerts.size();
    
    node_map_wo_ghosts_ = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);
    



    result = mbcore->tag_get_data(gid_tag,NotOwnedVerts,&(vert_gids[nvert]));
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }
    
    nvert += NotOwnedVerts.size();

    node_map_w_ghosts_ = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);

  }
  else {
    vert_gids = new int[AllVerts.size()];

    result = mbcore->tag_get_data(gid_tag,AllVerts,vert_gids);
    if (result != MB_SUCCESS) {
      std::cerr << "Problem getting tag data" << std::endl;
      assert(result == MB_SUCCESS);
    }

    nvert = AllVerts.size();

    node_map_wo_ghosts_ = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);
  }

  delete [] vert_gids;

}


unsigned int Mesh_MOAB::GID(Entity_ID lid, Entity_kind kind) const {
  MBEntityHandle ent;
  Entity_ID gid;

  switch (kind) {
  case NODE:
    ent = vtx_id_to_handle[lid];
    break;

  case FACE:
    ent = face_id_to_handle[lid];
    break;

  case CELL:
    ent = cell_id_to_handle[lid];
    break;
  default:
    std::cerr << "Global ID requested for unknown entity type" << std::endl;
  }

  int result = mbcore->tag_get_data(gid_tag,&ent,1,&gid);
  if (result != MB_SUCCESS) {
    std::cerr << "Problem getting tag data" << std::endl;
    assert(result == MB_SUCCESS);
  }

  return gid;
}



inline const Epetra_Map& Mesh_MOAB::cell_epetra_map (bool include_ghost) const {
  if (serial_run)
    return *cell_map_wo_ghosts_;
  else
    return (include_ghost ? *cell_map_w_ghosts_ : *cell_map_wo_ghosts_);
}


inline const Epetra_Map& Mesh_MOAB::face_epetra_map (bool include_ghost) const {
  if (serial_run)
    return *face_map_wo_ghosts_;
  else
    return (include_ghost ? *face_map_w_ghosts_ : *face_map_wo_ghosts_);
}

inline const Epetra_Map& Mesh_MOAB::node_epetra_map (bool include_ghost) const {
  if (serial_run)
    return *node_map_wo_ghosts_;
  else
    return (include_ghost ? *node_map_w_ghosts_ : *node_map_wo_ghosts_);
}


  // Get parallel type of eneity
  
Parallel_type Mesh_MOAB::entity_get_ptype(const Entity_kind kind, 
				 const Entity_ID entid) const
  {
    throw std::exception(); // Not implemented
  }
  
  
  
  
  // Get cell type
  
Cell_type Mesh_MOAB::cell_get_type(const Entity_ID cellid) const
  {
    return HEX;
  }
    
        
    
  

} // close namespace AmanziMesh
} // close namespace Amanzi
