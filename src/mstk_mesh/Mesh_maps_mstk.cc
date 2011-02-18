#include "Mesh_maps_mstk.hh"
//#include <Teuchos_RCP.hpp>

using namespace std;


  // Constructor - load up mesh from file

Mesh_maps_mstk::Mesh_maps_mstk (const char *filename, MPI_Comm incomm)
{
  int ok;
  int ring = 1; // One layer of ghost cells in parallel meshes
  int with_attr = 1;  // update of attributes in parallel meshes

  clear_internals_();

  MSTK_Init();


  comm = incomm;
  epcomm = new Epetra_MpiComm(comm);
  MPI_Comm_rank(comm,&myprocid);
  MPI_Comm_size(comm,&numprocs);

  serial_run =  (!comm || numprocs == 1) ? true : false;

  if (serial_run) {

    // Load serial mesh

    mesh = MESH_New(F1);
    ok = MESH_ImportFromExodusII(mesh,filename);

    myprocid = 0;
  }
  else {
    Mesh_ptr globalmesh;
    int dim=3;
    
    if (myprocid == 0) {
      globalmesh = MESH_New(F1);
      ok = MESH_ImportFromExodusII(globalmesh,filename);
      
      dim = MESH_Num_Regions(globalmesh) ? 3 : 2;
      
      mesh = globalmesh;
    }
    else {
      mesh = MESH_New(UNKNOWN_REP);
      ok = 1;
    }

    int DebugWait=0;
    while (DebugWait);
    
    ok = ok & MSTK_Mesh_Distribute(&mesh,dim,ring,with_attr,myprocid,
					numprocs,comm);

    if (myprocid == 0)
      MESH_Delete(globalmesh);
  }

  if (!ok) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << "Failed to load " << filename << " on processor " << myprocid << std::endl;
    assert(ok);
  }


  // Dimension of space - always 3 for MSTK
      
  spacedim = 3; 

  // Highest topological dimension
    
  if (MESH_Num_Regions(mesh)) {
    celldim = 3;
    facedim = 2;
  }
  else if (MESH_Num_Faces(mesh)) {
    celldim = 2;
    facedim = 1;
  }
  else {
    std::cerr << "Flow code works only on 2D and 3D meshes" << std::endl;
    assert(MESH_Num_Faces(mesh) > 0);
  }
      



  { // Keep together and in this order 

    init_pvert_lists();
    init_pcell_lists(); // cells MUST be initialized before faces
    init_pface_lists();

    
    // Create maps from local IDs to MSTK entity handles (must be after
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


Mesh_maps_mstk::~Mesh_maps_mstk() {
  delete cell_map_wo_ghosts_;
  delete cell_map_w_ghosts_;
  delete face_map_wo_ghosts_;
  delete face_map_w_ghosts_;
  delete node_map_wo_ghosts_;
  delete node_map_w_ghosts_;
  delete [] matset_ids;
  delete [] sideset_ids;
  delete [] nodeset_ids;
  delete [] faceflip;
  delete epcomm;

  if (AllVerts) MSet_Delete(AllVerts);
  if (OwnedVerts) MSet_Delete(OwnedVerts);
  if (NotOwnedVerts) MSet_Delete(NotOwnedVerts);
  if (AllFaces) MSet_Delete(AllFaces);
  if (OwnedFaces) MSet_Delete(OwnedFaces);
  if (NotOwnedFaces) MSet_Delete(NotOwnedFaces);
  if (AllCells) MSet_Delete(AllCells);
  if (OwnedCells) MSet_Delete(OwnedCells);
  if (GhostCells) MSet_Delete(GhostCells);
    
  // MESH_Delete(mesh);
}


// Some initializations

void Mesh_maps_mstk::clear_internals_ () 
{ 
  spacedim = 3;
  celldim = -1;
  facedim = -1;

  faceflip = NULL;

  cell_map_w_ghosts_ = cell_map_wo_ghosts_ = NULL;
  face_map_w_ghosts_ = face_map_wo_ghosts_ = NULL;
  node_map_w_ghosts_ = node_map_wo_ghosts_ = NULL;

  nmatsets = 0;
  nsidesets = 0;
  nnodesets = 0;
  matset_ids = sideset_ids = nodeset_ids = NULL;

}


void Mesh_maps_mstk::init_id_handle_maps() {
  int i, lid, nv, nf, nc, idx;
  MVertex_ptr vtx;
  MEntity_ptr face, cell;

  // If the mesh is dynamic, then this code has to be revisited
  
  // Amanzi has IDs starting from 0, MSTK has IDs starting from 1

  nv = List_Num_Entries(AllVerts);

  vtx_id_to_handle.reserve(nv);

  idx = 0; lid = 1;
  while ((vtx = MSet_Next_Entry(OwnedVerts,&idx))) {
    MEnt_Set_ID(vtx,lid);
    vtx_id_to_handle[lid-1] = vtx;
    lid++;
  }
    
  idx = 0;
  while ((vtx = MSet_Next_Entry(NotOwnedVerts,&idx))) {
    MEnt_Set_ID(vtx,lid);
    vtx_id_to_handle[lid-1] = vtx;
    lid++;
  }
    

  nf = MSet_Num_Entries(AllFaces);

  face_id_to_handle.reserve(nf);

  idx = 0; lid = 1;
  while ((face = MSet_Next_Entry(OwnedFaces,&idx))) {
    MEnt_Set_ID(face,lid);
    face_id_to_handle[lid-1] = face;
    lid++;
  }
  
  idx = 0;
  while ((face = MSet_Next_Entry(NotOwnedFaces,&idx))) {
    MEnt_Set_ID(face,lid);
    face_id_to_handle[lid-1] = face;
    lid++;
  }
    


  nc = MSet_Num_Entries(AllCells);

  cell_id_to_handle.reserve(nc);

  idx = 0; lid = 1;
  while ((cell = MSet_Next_Entry(OwnedCells,&idx))) {
    MEnt_Set_ID(cell,lid);
    cell_id_to_handle[lid-1] = cell;
    lid++;
  }
    
  idx = 0;
  while ((cell = MSet_Next_Entry(GhostCells,&idx))) {
    MEnt_Set_ID(cell,lid);
    cell_id_to_handle[lid-1] = cell;
    lid++;
  }
    
}



void Mesh_maps_mstk::init_global_ids() {
    
}


void Mesh_maps_mstk::init_pvert_lists() {
  int idx = 0;
  MVertex_ptr vtx;

  // Get all vertices on this processor 

  AllVerts = MSet_New(mesh,"AllVerts",MVERTEX);
  NotOwnedVerts = MSet_New(mesh,"NotOwnedVerts",MVERTEX);
  OwnedVerts = MSet_New(mesh,"OwnedVerts",MVERTEX);

  idx = 0;
  while ((vtx = MESH_Next_Vertex(mesh,&idx))) {
    MSet_Add(AllVerts,vtx);
    if (MV_PType(vtx) == PGHOST)
      MSet_Add(NotOwnedVerts,vtx);
    else
      MSet_Add(OwnedVerts,vtx);
  }

}


void Mesh_maps_mstk::init_pface_lists() {
  int idx = 0;
  MFace_ptr face;
  MEdge_ptr edge;

  // Get all faces on this processor 

  if (facedim == 2) {
    AllFaces = MSet_New(mesh,"AllFaces",MFACE);
    NotOwnedFaces = MSet_New(mesh,"NotOwnedFaces",MFACE);
    OwnedFaces = MSet_New(mesh,"OwnedFaces",MFACE);

    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      MSet_Add(AllFaces,face);
      if (MF_PType(face) == PGHOST)
	MSet_Add(NotOwnedFaces,face);
      else
	MSet_Add(OwnedFaces,face);
    }
  }
  else {
    std::cerr << "Not implemented for 2D" << std::endl;
  }
}


void Mesh_maps_mstk::init_pface_dirs() {
    MRegion_ptr region0, region1;
    MFace_ptr face;
    MAttrib_ptr attfc0, attfc1;
    int idx;
    int local_regid0, local_regid1;
    int remote_regid0, remote_regid1;


    if (serial_run) {
      faceflip = new bool[MSet_Num_Entries(AllFaces)];
      for (int i = 0; i < MSet_Num_Entries(AllFaces); i++) faceflip[i] = false;
    }
    else {
      // Do some additional processing to see if ghost faces and their masters
      // are oriented the same way; if not, turn on flag to flip them

      if (facedim == 2) {
	attfc0 = MAttrib_New(mesh,"TMP_FC0_ATT",INT,MFACE);
	attfc1 = MAttrib_New(mesh,"TMP_FC1_ATT",INT,MFACE);
    
	idx = 0;
	while ((face = MESH_Next_Face(mesh,&idx))) {
	  if (MF_PType(face) != PINTERIOR) {
	    region0 = MF_Region(face,0);
	    if (region0)
	      MEnt_Set_AttVal(face,attfc0,MEnt_GlobalID(region0),0.0,NULL);

	    region1 = MF_Region(face,1);
	    if (region1)
	      MEnt_Set_AttVal(face,attfc1,MEnt_GlobalID(region1),0.0,NULL);
	  }
	}    
      }



      MSTK_UpdateAttr(mesh, myprocid, numprocs, comm);



      faceflip = new bool[MSet_Num_Entries(AllFaces)];
      for (int i = 0; i < MSet_Num_Entries(AllFaces); i++) faceflip[i] = false;
    
      if (facedim == 2) {
	double rval;
	void *pval;

	idx = 0;
	while ((face = MSet_Next_Entry(NotOwnedFaces,&idx))) {
      
	  MEnt_Get_AttVal(face,attfc0,&remote_regid0,&rval,&pval);
	  MEnt_Get_AttVal(face,attfc1,&remote_regid1,&rval,&pval);
      
	  region0 = MF_Region(face,0);
	  local_regid0 = region0 ? MEnt_GlobalID(region0) : 0;
	  region1 = MF_Region(face,1);
	  local_regid1 = region1 ? MEnt_GlobalID(region1) : 0;
      
	  if (remote_regid1 == local_regid0 || 
	      remote_regid0 == local_regid1) {
	    int lid = MEnt_ID(face);
	    faceflip[lid-1] = true;
	  }
	  else { // Sanity Check
	
	    if (remote_regid1 != local_regid1 &&
		remote_regid0 != local_regid0) {
	  
	      cout << "Face cells mismatch between master and ghost (processor " << myprocid << ")" << std::endl;
	      cout << " Face " << MEnt_GlobalID(face) << std::endl;
	      cout << "Remote cells " << remote_regid0 << " " << remote_regid1 << std::endl;
	      cout << "Local cells " << local_regid0 << " " << local_regid1 << std::endl;
	    }
	  }
	}
      }
      else {
	cout << "init_pface_dir not implemented for 2D yet" << std::endl;
	assert(facedim != 2);
      }
    }    
}


void Mesh_maps_mstk::init_pcell_lists() {
  int idx = 0;

  if (celldim == 3) {
    MRegion_ptr region;

    AllCells = MSet_New(mesh,"AllCells",MREGION);
    OwnedCells = MSet_New(mesh,"OwnedCells",MREGION);
    GhostCells = MSet_New(mesh,"GhostCells",MREGION);

    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {
      MSet_Add(AllCells,region);
      if (MR_PType(region) == PGHOST)
	MSet_Add(GhostCells,region);
      else
	MSet_Add(OwnedCells,region);
    }
  }
  else {
    std::cerr << "Not implemented for 2D" << std::endl;
    assert (celldim == 3);
  }

}



void Mesh_maps_mstk::init_set_info() {
  int idx, i, j, ip, gid, found, setid;
  int maxmatsets, maxsidesets, maxnodesets;
  std::vector<int> loc_matset_ids, loc_sideset_ids, loc_nodeset_ids;
  MRegion_ptr mr;
  MFace_ptr mf;
  MVertex_ptr mv;
  MSet_ptr mset;
  char msetname[256];

  // Get material, sideset and nodeset tags

  idx = 0;
  while ((mset = MESH_Next_MSet(mesh,&idx))) {
    MSet_Name(mset,msetname);
    if (strncmp(msetname,"matset_",6) == 0) {
      sscanf(&(msetname[7]),"%d",&setid);
      loc_matset_ids.push_back(setid);
    }
  }
  nmatsets = loc_matset_ids.size();



  idx = 0;
  while ((mset = MESH_Next_MSet(mesh,&idx))) {
    MSet_Name(mset,msetname);
    if (strncmp(msetname,"sideset_",6) == 0) {
      sscanf(&(msetname[8]),"%d",&setid);
      loc_sideset_ids.push_back(setid);
    }
  }  
  nsidesets = loc_sideset_ids.size();



  idx = 0;
  while ((mset = MESH_Next_MSet(mesh,&idx))) {
    MSet_Name(mset,msetname);
    if (strncmp(msetname,"nodeset_",6) == 0) {
      sscanf(&(msetname[8]),"%d",&setid);
      loc_nodeset_ids.push_back(setid);
    }
  }
  nnodesets = loc_nodeset_ids.size();



  if (!serial_run) {

    // Figure out the maximum number of sets on any processor,

    MPI_Allreduce(&nmatsets,&maxmatsets,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&nsidesets,&maxsidesets,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&nnodesets,&maxnodesets,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);



    // We'll pad maxnsets so that we can try and avoid (but not
    // guarantee) reallocation
    
    maxmatsets *= 2;
    maxsidesets *= 2;
    maxnodesets *= 2;
    
    matset_ids = new int[maxmatsets];
    sideset_ids = new int[maxsidesets];
    nodeset_ids = new int[maxnodesets];

    for (i = 0; i < maxmatsets; i++)
      matset_ids[i] = (i < nmatsets) ? loc_matset_ids[i] : 0;
    for (i = 0; i < maxsidesets; i++)
      sideset_ids[i] = (i < nsidesets) ? loc_sideset_ids[i] : 0;
    for (i = 0; i < maxnodesets; i++)
      nodeset_ids[i] = (i < nnodesets) ? loc_nodeset_ids[i] : 0;


    // Collate all the sets across processors

  
    int *allmatset_ids = new int[numprocs*maxmatsets];
    int *allsideset_ids = new int[numprocs*maxsidesets];
    int *allnodeset_ids = new int[numprocs*maxnodesets];
  
  
    MPI_Allgather(matset_ids,maxmatsets,MPI_INT,allmatset_ids,maxmatsets,
		  MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(sideset_ids,maxsidesets,MPI_INT,allsideset_ids,maxsidesets,
		  MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(nodeset_ids,maxnodesets,MPI_INT,allnodeset_ids,maxnodesets,
		  MPI_INT,MPI_COMM_WORLD);
  
  
  
  // Make a list of all setids/entity dimensions across all processors

  if (myprocid == 0) {

    // Just add to the list of setids on processor 0. Since we doubled
    // the computed maxnsets and allocated allsetids to be of length
    // nprocs*maxnsets, we are guaranteed that we will not overwrite
    // processor 1 entries if we add unique entries from processor 1
    // to the end of processor 0's list of sets and so on

    // Mat sets

    for (ip = 1; ip < numprocs; ip++) {
      for (i = 0; i < maxmatsets; i++) {
	if (!allmatset_ids[ip*maxmatsets+i]) continue;

	found = 0;
	for (j = 0; j < nmatsets; j++) {
	  if (allmatset_ids[j] == allmatset_ids[ip*maxmatsets+i]) {
	    found = 1;
	    break;
	  }
	}
	  
	if (!found) {
	  allmatset_ids[nmatsets] = allmatset_ids[ip*maxmatsets+i];
	  nmatsets++;
	}
      }
    }

    if (nmatsets > maxmatsets) {
      // We have to resize the matset_ids array
      delete [] matset_ids;
      maxmatsets = nmatsets;
      matset_ids = new int[maxmatsets];
    }
    
    memcpy(matset_ids, allmatset_ids, nmatsets*sizeof(allmatset_ids[0]));


    // Sidesets


    for (ip = 1; ip < numprocs; ip++) {
      for (i = 0; i < maxsidesets; i++) {
	found = 0;
	for (j = 0; j < nsidesets; j++) {
	  if (allsideset_ids[j] == allsideset_ids[ip*maxsidesets+i]) {
	    found = 1;
	    break;
	  }
	}
	  
	if (!found) {
	  allsideset_ids[nsidesets] = allsideset_ids[ip*maxsidesets+i];
	  nsidesets++;
	}
      }
    }

    if (nsidesets > maxsidesets) {
      // We have to resize the sideset_ids array
      delete [] sideset_ids;
      maxsidesets = nsidesets;
      sideset_ids = new int[maxsidesets];
    }
    
    memcpy(sideset_ids, allsideset_ids, nsidesets*sizeof(allsideset_ids[0]));


    // Nodesets

    for (ip = 1; ip < numprocs; ip++) {
      for (i = 0; i < maxnodesets; i++) {
	found = 0;
	for (j = 0; j < nnodesets; j++) {
	  if (allnodeset_ids[j] == allnodeset_ids[ip*maxnodesets+i]) {
	    found = 1;
	    break;
	  }
	}
	  
	if (!found) {
	  allnodeset_ids[nnodesets] = allnodeset_ids[ip*maxnodesets+i];
	  nnodesets++;
	}
      }
    }

    if (nnodesets > maxnodesets) {
      // We have to resize the nodeset_ids array
      delete [] nodeset_ids;
      maxnodesets = nnodesets;
      nodeset_ids = new int[maxnodesets];
    }
    
    memcpy(nodeset_ids, allnodeset_ids, nnodesets*sizeof(allnodeset_ids[0]));

  } // if myprocid == 0


  // Tell all the other processors about number of sets to expect

  MPI_Bcast(&nmatsets,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nsidesets,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nnodesets,1,MPI_INT,0,MPI_COMM_WORLD);

  
  // Make sure the other processors have enough space allocated to receive
  // the expanded list of sets

  if (myprocid != 0) {					       
    if (nmatsets > maxmatsets) { 
      delete [] matset_ids;
      maxmatsets = nmatsets;
      matset_ids = new int[maxmatsets];
    }    

    if (nsidesets > maxsidesets) { 
      delete [] sideset_ids;
      maxsidesets = nsidesets;
      sideset_ids = new int[maxsidesets];
    }    
  }


  // Tell all the other processors about the complete list of sets
  
  MPI_Bcast(matset_ids,nmatsets,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(sideset_ids,nsidesets,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(nodeset_ids,nnodesets,MPI_INT,0,MPI_COMM_WORLD);
  

  delete [] allmatset_ids;
  delete [] allsideset_ids;
  delete [] allnodeset_ids;


  }
  else { /* serial run */

    matset_ids = new int[nmatsets];
    for (i = 0; i < nmatsets; i++)
      matset_ids[i] = loc_matset_ids[i];

    sideset_ids = new int[nsidesets];
    for (i = 0; i < nsidesets; i++)
      sideset_ids[i] = loc_sideset_ids[i];

    nodeset_ids = new int[nnodesets];
    for (i = 0; i < nnodesets; i++)
      nodeset_ids[i] = loc_nodeset_ids[i];

  }


  // Now create the sets (including empty sets)
  
  for (i = 0; i < nmatsets; i++) {
    stringstream setname;
  
    setname << "matset_" << matset_ids[i];
    mset = MESH_MSetByName(mesh,setname.str().c_str());
    if (!mset) { 
      if (celldim == 3) { // create an empty set
	mset = MSet_New(mesh,setname.str().c_str(),MREGION);
      }
      else if (celldim == 2) {
	std::cerr << "Not implemented for 2D" << std::endl;
	assert (celldim == 3);
      }
    }
  }
  

  for (i = 0; i < nsidesets; i++) {
    stringstream setname;
  
    setname << "sideset_" << sideset_ids[i];
    mset = MESH_MSetByName(mesh,setname.str().c_str());
    if (!mset) {
      if (facedim == 2) {
	mset = MSet_New(mesh,setname.str().c_str(),MFACE);
      }
      else if (facedim == 1) {
	std::cerr << "Not implemented for 2D" << std::endl;
	assert (facedim == 2);
      }
    }
  }



  for (i = 0; i < nnodesets; i++) {
    stringstream setname;
  
    setname << "nodeset_" << nodeset_ids[i];
    mset = MESH_MSetByName(mesh,setname.str().c_str());
    if (!mset) {
      mset = MSet_New(mesh,setname.str().c_str(),MVERTEX);
    }
  }

}




// Number of OWNED, GHOST or USED entities of different types

    
unsigned int Mesh_maps_mstk::count_entities (Mesh_data::Entity_kind kind, 
					     Element_Category category) const
{
  const int rank = (int) kind;
  const int index = ((int) category) - 1;

  switch (kind) {
  case Mesh_data::NODE:

    switch (category) {
    case OWNED:
      return !serial_run ? MSet_Num_Entries(OwnedVerts) : MSet_Num_Entries(AllVerts);
      break;
    case GHOST:
      return !serial_run ? MSet_Num_Entries(NotOwnedVerts) : 0;
      break;
    case USED:
      return MSet_Num_Entries(AllVerts);
      break;
    default:
      return 0;
    }
    break;


  case Mesh_data::FACE:
    switch (category) {
    case OWNED:
      return !serial_run ? MSet_Num_Entries(OwnedFaces) : MSet_Num_Entries(AllFaces);
      break;
    case GHOST:
      return !serial_run ? MSet_Num_Entries(NotOwnedFaces) : 0;
      break;
    case USED:
      return MSet_Num_Entries(AllFaces);
      break;
    default:
      return 0;
    }
    break;


  case Mesh_data::CELL:

    switch (category) {
    case OWNED:
      return !serial_run ? MSet_Num_Entries(OwnedCells) : MSet_Num_Entries(AllCells);
      break;
    case GHOST:
      return !serial_run ? MSet_Num_Entries(GhostCells) : 0;
      break;
    case USED:
      return MSet_Num_Entries(AllCells);
      break;
    default:
      return 0;
    }

    break;
  default:
    std::cerr << "Count requested for unknown entity type" << std::endl;
  }
}


void Mesh_maps_mstk::cell_to_faces (unsigned int cellid, 
				    std::vector<unsigned int>::iterator begin, 
				    std::vector<unsigned int>::iterator end)
{
  cell_to_faces(cellid, &(*begin), &(*end));
}

void Mesh_maps_mstk::cell_to_faces (unsigned int cellid, 
				    unsigned int *begin, 
				    unsigned int *end)
{
  int nf;
  MEntity_ptr cell;
  int *cell_faceids;
  int cfstd[6][4] = {{0,1,5,4},    // Expected cell-face-node pattern
		     {1,2,6,5},
		     {2,3,7,6},
		     {0,4,7,3},
		     {0,3,2,1},
		     {4,5,6,7}};

  cell = cell_id_to_handle[cellid];
      
  if (celldim == 3) {
    
    List_ptr rverts = MR_Vertices((MRegion_ptr)cell);
    List_ptr rfaces = MR_Faces((MRegion_ptr)cell);

    nf = List_Num_Entries(rfaces);
    cell_faceids = new int[nf];

    if (nf == 6) {
      // Have to re-sort the faces according a specific template for hexes
      
      for (int i = 0; i < 6; i++) {
	
	// Search for a face that has all the expected nodes
	
	bool found = false;
	  
	for (int j = 0; j < 6; j++) {
	  
	  bool all_present = true;
	  MFace_ptr face = List_Entry(rfaces,j);	  
	  List_ptr fverts = MF_Vertices(face,1,0);
	  
	  // Check if this face has all the expected nodes
	  
	  for (int k = 0; k < 4; k++) {
	    MVertex_ptr vtx = List_Entry(rverts,cfstd[i][k]);
	    
	    if (!List_Contains(fverts,vtx)) {
	      all_present = false;
	      break;
	    }
	  }

	  if (all_present) {
	    int lid = MEnt_ID(face);
	    cell_faceids[i] = lid-1;
	    found = true;
	    break;
	  }

	  List_Delete(fverts);
	  
	} // for (int j = 0; j < 6; j++) 

	if (!found) {
	  std::cerr << "Could not find face in hex" << std::endl;
	  assert(found);
	}	   

      } // for (int i = 0; i < 6; i++) 

    } // If (List_Num_Entries(rfaces) == 6)

    List_Delete(rfaces);
    List_Delete(rverts);

  } // If celldim == 3
  else {
    std::cerr << "Not implemented" << std::endl;
    assert (celldim != 3);
  }
  
  std::copy (cell_faceids, cell_faceids+nf, begin);

  delete [] cell_faceids;
}




void Mesh_maps_mstk::cell_to_face_dirs (unsigned int cellid, 
				       std::vector<int>::iterator begin, 
				       std::vector<int>::iterator end){
  cell_to_face_dirs(cellid, &(*begin), &(*end));
}

void Mesh_maps_mstk::cell_to_face_dirs (unsigned int cellid, 
				       int *begin, 
				       int *end)
{
  MEntity_ptr cell;
  unsigned int *cell_faces;
  int *cell_facedirs;
  int j,nf, result;
  

  if (celldim == 3) {
    nf = 6; // HEX SPECIFIC - THIS NUMBER SHOULD COME FROM cell_to_faces

    cell = cell_id_to_handle[cellid];

    cell_faces = new unsigned int[nf];

    cell_to_faces(cellid,cell_faces,cell_faces+nf);

    assert ((unsigned int) (end - begin) >= nf);
    cell_facedirs = new int[nf];

    for (int i = 0; i < nf; i++) {
      MFace_ptr face = face_id_to_handle[cell_faces[i]];

      cell_facedirs[i] = MR_FaceDir(cell,face) == 1 ? 1 : -1;

      // If this is a ghost face and the master has the opposite direction
      // we are supposed to flip it

      if (faceflip[cell_faces[i]]) cell_facedirs[i] *= -1;
    }
    
    std::copy (cell_facedirs, cell_facedirs+nf, begin);

    delete [] cell_faces;
    delete [] cell_facedirs;
  }
  else {
    std::cerr << "Not implemented in 2D" << std::endl;
    assert(celldim == 3);
  }
}



void Mesh_maps_mstk::cell_to_nodes (unsigned int cellid, 
				    std::vector<unsigned int>::iterator begin, 
				    std::vector<unsigned int>::iterator end)
{
  cell_to_nodes(cellid, &(*begin), &(*end));
}



void Mesh_maps_mstk::cell_to_nodes (unsigned int cellid, unsigned int *begin, 
				    unsigned int *end)
{
  MEntity_ptr cell;
  int *cell_nodeids;
  int nn, lid;


  cell = cell_id_to_handle[cellid];
      
  if (celldim == 3) {
    List_ptr rverts = MR_Vertices(cell);
 
    nn = List_Num_Entries(rverts);
    assert ((unsigned int) (end - begin) >= nn);
    
    cell_nodeids = new int[nn];
    
    for (int i = 0; i < nn; i++) {
      lid = MEnt_ID(List_Entry(rverts,i));
      cell_nodeids[i] = lid-1;
    }
    
    List_Delete(rverts);
  }
  else {
    List_ptr fverts = MF_Vertices(cell,1,0);

    nn = List_Num_Entries(fverts);
    assert ((unsigned int) (end - begin) >= nn);

    cell_nodeids = new int[nn];
    
    for (int i = 0; i < nn; i++) {
      lid = MEnt_ID(List_Entry(fverts,i));
      cell_nodeids[i] = lid-1;
    }
    
    List_Delete(fverts);
  }

  std::copy (cell_nodeids, cell_nodeids+nn, begin);

  delete [] cell_nodeids;
}




void Mesh_maps_mstk::face_to_nodes (unsigned int faceid, 
				    std::vector<unsigned int>::iterator begin, 
				    std::vector<unsigned int>::iterator end)
{
  face_to_nodes(faceid, &(*begin), &(*end));
}

  
void Mesh_maps_mstk::face_to_nodes (unsigned int faceid, unsigned int *begin, 
				    unsigned int *end) {
  MEntity_ptr face;
  int *face_nodeids;
  int nn, lid;

  face = face_id_to_handle[faceid];

  if (facedim == 2) {
    
    List_ptr fverts = MF_Vertices(face,1,0);
    assert(fverts != NULL);

    nn = List_Num_Entries(fverts);
    assert ((unsigned int) (end - begin) >= nn);
    
    face_nodeids = new int[nn];
    
    if (faceflip[faceid]) {
      for (int i = nn-1; i >= 0; i--) {
	lid = MEnt_ID(List_Entry(fverts,i));
	face_nodeids[nn-i-1] = lid-1;
      }
    }  
    else {
      for (int i = 0; i < nn; i++) {
	lid = MEnt_ID(List_Entry(fverts,i));
	face_nodeids[i] = lid-1;
      }
    }
  }
  else {
    std::cerr << "Not implemented for 2D" << std::endl;
    assert(facedim == 2);
  }

  std::copy (face_nodeids, face_nodeids+nn, begin);

  delete [] face_nodeids;
}
  


void Mesh_maps_mstk::node_to_coordinates (unsigned int node_id, 
					  std::vector<double>::iterator begin, 
					  std::vector<double>::iterator end)
{
  node_to_coordinates(node_id, &(*begin), &(*end));
}

void Mesh_maps_mstk::node_to_coordinates (unsigned int node_id, 
					  double *begin, 
					  double *end)
{
  MEntity_ptr vtx;

  assert ((unsigned int) (end-begin) >= spacedim);
    
  double coords[3];

  vtx = vtx_id_to_handle[node_id];

  MV_Coords(vtx,coords);

  std::copy (coords, coords+spacedim, begin);

}



void Mesh_maps_mstk::cell_to_coordinates (unsigned int cellid, 
					  std::vector<double>::iterator begin, 
					  std::vector<double>::iterator end)
{
  cell_to_coordinates(cellid, &(*begin), &(*end));
}


void Mesh_maps_mstk::cell_to_coordinates (unsigned int cellid, 
					  double *begin, double *end)
{
  MEntity_ptr cell;
  double *coords;
  int nn, result;


  if (celldim == 3) {
    cell = cell_id_to_handle[cellid];
      
    List_ptr rverts = MR_Vertices(cell);
    
    nn = List_Num_Entries(rverts);
    assert ((unsigned int) (end - begin) >= nn);

    coords = new double[spacedim*nn];

    for (int i = 0; i < nn; i++)
      MV_Coords(List_Entry(rverts,i),&(coords[i*spacedim]));
    
    List_Delete(rverts);

    std::copy (coords, coords+spacedim*nn, begin);

    delete [] coords;
  }
  else {
    std::cerr << "Not implemented for 2D" << std::endl;
    assert (celldim == 3);
  }
}




void Mesh_maps_mstk::face_to_coordinates (unsigned int faceid, 
					  std::vector<double>::iterator begin, 
					  std::vector<double>::iterator end)
{
  face_to_coordinates(faceid, &(*begin), &(*end));
}
  
void Mesh_maps_mstk::face_to_coordinates (unsigned int faceid, 
					  double *begin, 
					  double *end)
{
    MEntity_ptr face;
    double *coords;
    int nn;

    face = face_id_to_handle[faceid];

    if (facedim == 2) {
      List_ptr fverts = MF_Vertices(face,1,0);

      nn = List_Num_Entries(fverts);
      assert ((unsigned int) (end - begin) >= nn);

      coords = new double[spacedim*nn];
    
      if (faceflip[faceid]) {
	for (int i = nn-1; i >=0; i--)
	  MV_Coords(List_Entry(fverts,i),coords+spacedim*(nn-i-1));
      }
      else {
	for (int i = 0; i < nn; i++)
	  MV_Coords(List_Entry(fverts,i),coords+spacedim*i);
      }
    }
    else {
      std::cerr << "Not implemented in 2D" << std::endl;
      assert (facedim == 2);
    }

    std::copy (coords, coords+spacedim*nn, begin);

    delete [] coords;
}
  


void Mesh_maps_mstk::get_set_ids (Mesh_data::Entity_kind kind, 
			      std::vector<unsigned int>::iterator begin, 
			      std::vector<unsigned int>::iterator end) const {
  get_set_ids(kind, &(*begin), &(*end));
}


void Mesh_maps_mstk::get_set_ids (Mesh_data::Entity_kind kind, 
			      unsigned int *begin, unsigned int *end) const {
  int i, ns=0;

  switch (kind) {
  case Mesh_data::CELL: {
    int sids[nmatsets];

    for (i = 0; i < nmatsets; i++)
      sids[ns++] = matset_ids[i];

    assert ((unsigned int) (end - begin) >= ns);

    std::copy (sids, sids + ns, begin);

    return;
    break;
  }

  case Mesh_data::FACE: {
    int sids[nsidesets];

    for (i = 0; i < nsidesets; i++)
      sids[ns++] = sideset_ids[i];

    assert ((unsigned int) (end - begin) >= ns);

    std::copy (sids, sids + ns, begin);

    return;
    break;
  }

  case Mesh_data::NODE: {
    int sids[nnodesets];

    for (i = 0; i < nnodesets; i++)
	sids[ns++] = nodeset_ids[i];

    assert ((unsigned int) (end - begin) >= ns);

    std::copy (sids, sids + ns, begin);

    return;
    break;
  }

  default:
    return;
  }
    
}



void Mesh_maps_mstk::get_set (unsigned int set_id, 
			      Mesh_data::Entity_kind kind, 
			      Element_Category category,
			      std::vector<unsigned int>::iterator begin, 
			      std::vector<unsigned int>::iterator end) const {
  get_set(set_id, kind, category, &(*begin), &(*end));
}

void Mesh_maps_mstk::get_set (unsigned int set_id, 
			      Mesh_data::Entity_kind kind, 
			      Element_Category category,
			      unsigned int *begin, unsigned int *end) const {

  int idx, i, lid;
  MSet_ptr mset, mset1;
  MEntity_ptr ment;
  stringstream setname;

  switch (kind) {
  case Mesh_data::CELL: {
    int ncells;
    int *cells;

    setname << "matset_" << set_id;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    ncells = MSet_Num_Entries(mset1);
    if (ncells == 0) return;

    if (category == OWNED)
      mset = MSets_Subtract(mset1,GhostCells);
    else if (category == GHOST)
      mset = MSets_Subtract(mset1,OwnedCells);
    else
      mset = mset1;

    if (!mset) return;

    ncells = MSet_Num_Entries(mset);

    assert ((unsigned int) (end - begin) >= ncells);


    cells = new int[ncells];
    i = 0; idx = 0;
    while ((ment = MSet_Next_Entry(mset,&idx))) {
      lid = MEnt_ID(ment);
      cells[i++] = lid-1;
    }

    std::copy(cells, cells + ncells, begin);

    delete [] cells;
    if (category != USED)
      MSet_Delete(mset);

    return;
    break;
  }
  case Mesh_data::FACE: {
    int nfaces;
    int *faces;

    setname << "sideset_" << set_id;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    nfaces = MSet_Num_Entries(mset1);
    if (nfaces == 0) return;

    if (category == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedFaces);
    else if (category == GHOST)
      mset = MSets_Subtract(mset1,OwnedFaces);
    else
      mset = mset1;

    if (!mset) return;

    nfaces = MSet_Num_Entries(mset);

    assert ((unsigned int) (end - begin) >= nfaces);

    faces = new int[nfaces];
    i = 0; idx = 0;
    while ((ment = MSet_Next_Entry(mset,&idx))) {
      lid = MEnt_ID(ment);
      faces[i++] = lid-1;
    }

    std::copy(faces, faces + nfaces, begin);
      
    delete [] faces;
    if (category != USED)
      MSet_Delete(mset);
    return;
    break;
  }
  case Mesh_data::NODE: {
    int nnodes;
    int *nodes;

    setname << "nodeset_" << set_id;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    nnodes = MSet_Num_Entries(mset1);
    if (nnodes == 0) return;
    
    if (category == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedVerts);
    else if (category == GHOST)
      mset = MSets_Subtract(mset1,OwnedVerts);
    else
      mset = mset1;

    if (!mset) return;

    nnodes = MSet_Num_Entries(mset);

    assert ((unsigned int) (end - begin) >= nnodes);

    nodes = new int[nnodes];
    i = 0; idx = 0;
    while ((ment = MSet_Next_Entry(mset,&idx))) {
      lid = MEnt_ID(ment);
      nodes[i++] = lid-1;
    }

    std::copy(nodes, nodes + nnodes, begin);
      
    delete [] nodes;

    if (category != USED)
      MSet_Delete(mset);

    return;
    break;
  }
  default:
    return;
  }
  
}


unsigned int Mesh_maps_mstk::num_sets(Mesh_data::Entity_kind kind) const {
  int n = 0;
    
  switch (kind) {
  case Mesh_data::CELL:
    return nmatsets;
  case Mesh_data::FACE:
    return nsidesets;
  case Mesh_data::NODE:
    return nnodesets;
  default:
    return 0;
  }
}

bool Mesh_maps_mstk::valid_set_id (unsigned int id, Mesh_data::Entity_kind kind) const {
  int n = 0;

  switch (kind) {
  case Mesh_data::CELL:
    for (int i = 0; i < nmatsets; i++)
      if (matset_ids[i] == id) return true;
    break;
  case Mesh_data::FACE:
    for (int i = 0; i < nsidesets; i++)
      if (sideset_ids[i] == id) return true;
    break;
  case Mesh_data::NODE:
    for (int i = 0; i < nnodesets; i++)
      if (nodeset_ids[i] == id) return true;
    break;
  default:
    return false;
  }    

  return false;
}

unsigned int Mesh_maps_mstk::get_set_size(unsigned int set_id, 
				      Mesh_data::Entity_kind kind, 
				      Element_Category category) const {   

  stringstream setname;
  MSet_ptr mset1, mset;
  int nent;

  switch (kind) {
  case Mesh_data::CELL: {
    int ncells;

    setname << "matset_" << set_id;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    ncells = MSet_Num_Entries(mset1);
    if (ncells == 0) return 0;

    if (category == OWNED)
      mset = MSets_Subtract(mset1,GhostCells);
    else if (category == GHOST)
      mset = MSets_Subtract(mset1,OwnedCells);
    else
      mset = mset1;

    if (!mset) return 0;

    nent = MSet_Num_Entries(mset);
    if (category != USED)
      MSet_Delete(mset);

    return nent;
    break;
  }

  case Mesh_data::FACE: {
    int nfaces;

    setname << "sideset_" << set_id;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    nfaces = MSet_Num_Entries(mset1);
    if (nfaces == 0) return 0;

    if (category == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedFaces);
    else if (category == GHOST)
      mset = MSets_Subtract(mset1,OwnedFaces);
    else
      mset = mset1;
    
    if (!mset) return 0;

    nent = MSet_Num_Entries(mset);
    if (category != USED)
      MSet_Delete(mset);

    return nent;
    break;
  }

  case Mesh_data::NODE: {
    int nnodes;

    setname << "nodeset_" << set_id;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    nnodes = MSet_Num_Entries(mset1);
    if (nnodes == 0) return 0;

    if (category == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedVerts);
    else if (category == GHOST)
      mset = MSets_Subtract(mset1,OwnedVerts);
    else
      mset = mset1;

    if (!mset) return 0;

    nent = MSet_Num_Entries(mset);
    if (category != USED)
      MSet_Delete(mset);

    return nent;
    break;
  }

  default:
    return 0;
  }
    
}


// Epetra map for cells - basically a structure specifying the
// global IDs of cells owned or used by this processor

// Amanzi/Epetra want global IDs to start at 0

void Mesh_maps_mstk::init_cell_map ()
{
  int *cell_gids;
  int ncell, idx, i;
  MEntity_ptr ment;

  if (!serial_run) {

    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells

    int nowned = MSet_Num_Entries(OwnedCells);
    int nnotowned = MSet_Num_Entries(GhostCells);

    cell_gids = new int[nowned+nnotowned];
    
    ncell = nowned;
   
    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(OwnedCells,&idx)))
      cell_gids[i++] = MEnt_GlobalID(ment)-1;

    cell_map_wo_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);
    


    ncell += nnotowned;

    idx = 0; 
    while ((ment = MSet_Next_Entry(GhostCells,&idx)))
      cell_gids[i++] = MEnt_GlobalID(ment)-1;
    
    cell_map_w_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);

  }
  else {
    ncell = MSet_Num_Entries(AllCells);
    cell_gids = new int[ncell];

    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(AllCells,&idx)))      
      cell_gids[i++] = MEnt_ID(ment)-1;

    cell_map_wo_ghosts_ = new Epetra_Map(-1,ncell,cell_gids,0,*epcomm);
  }

  delete [] cell_gids;

}




// Epetra map for faces - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_maps_mstk::init_face_map ()
{
  int *face_gids;
  int nface, idx, i;
  MEntity_ptr ment;

  if (!serial_run) {

    // For parallel runs create map without and with ghost cells included
    // Also, put in owned cells before the ghost cells

    int nowned = MSet_Num_Entries(OwnedFaces);
    int nnotowned = MSet_Num_Entries(NotOwnedFaces);

    face_gids = new int[nowned+nnotowned];
    
    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(OwnedFaces,&idx)))
      face_gids[i++] = MEnt_GlobalID(ment)-1;

    nface = nowned;
    
    face_map_wo_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*epcomm);


    idx = 0;
    while ((ment = MSet_Next_Entry(NotOwnedFaces,&idx))) 
      face_gids[i++] = MEnt_GlobalID(ment)-1;

    nface += nnotowned;

    face_map_w_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*epcomm);

  }
  else {
    nface = MSet_Num_Entries(AllFaces);
    face_gids = new int[nface];

    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(AllFaces,&idx)))
      face_gids[i++] = MEnt_ID(ment)-1;

    face_map_wo_ghosts_ = new Epetra_Map(-1,nface,face_gids,0,*epcomm);
  }

  delete [] face_gids;

}




// Epetra map for nodes - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_maps_mstk::init_node_map ()
{
  int *vert_gids;
  int nvert, idx, i;
  MEntity_ptr ment;

  if (!serial_run) {

    // For parallel runs create map without and with ghost verts included
    // Also, put in owned cells before the ghost verts

    int nowned = MSet_Num_Entries(OwnedVerts);
    int nnotowned = MSet_Num_Entries(NotOwnedVerts);

    vert_gids = new int[nowned+nnotowned];
    
    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(OwnedVerts,&idx)))
      vert_gids[i++] = MEnt_GlobalID(ment)-1;

    nvert = nowned;
    
    node_map_wo_ghosts_ = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);
    


    idx = 0;
    while ((ment = MSet_Next_Entry(NotOwnedVerts,&idx)))
      vert_gids[i++] = MEnt_GlobalID(ment)-1;

    nvert += nnotowned;

    node_map_w_ghosts_ = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);

  }
  else {
    nvert = MSet_Num_Entries(AllVerts);

    vert_gids = new int[nvert];

    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(AllVerts,&idx)))
      vert_gids[i++] = MEnt_ID(ment)-1;

    node_map_wo_ghosts_ = new Epetra_Map(-1,nvert,vert_gids,0,*epcomm);
  }

  delete [] vert_gids;

}


unsigned int Mesh_maps_mstk::GID(unsigned int lid, Mesh_data::Entity_kind kind) {
  MEntity_ptr ent;
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

  if (serial_run)
    return MEnt_ID(ent)-1;
  else
    return MEnt_GlobalID(ent)-1;
}



