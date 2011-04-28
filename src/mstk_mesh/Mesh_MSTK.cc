#include "Mesh_MSTK.hh"
//#include <Teuchos_RCP.hpp>

using namespace std;

using namespace Amanzi;
using namespace AmanziMesh;
using namespace AmanziGeometry;

// Constructor - load up mesh from file

Mesh_MSTK::Mesh_MSTK (const char *filename, MPI_Comm incomm, int space_dimension)
{
  int ok;
  int ring = 1; // One layer of ghost cells in parallel meshes
  int with_attr = 1;  // update of attributes in parallel meshes

  clear_internals_();

  MSTK_Init();

  set_space_dimension(space_dimension);

  comm = incomm;
  epcomm = new Epetra_MpiComm(comm);
  MPI_Comm_rank(comm,&myprocid);
  MPI_Comm_size(comm,&numprocs);

  serial_run =  (!comm || numprocs == 1) ? true : false;

  if (serial_run) {

    // Load serial mesh

    mesh = MESH_New(F1);
    ok = MESH_ImportFromExodusII(mesh,filename);

    set_cell_dimension((MESH_Num_Regions(mesh) ? 3 : 2));

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

    set_cell_dimension(dim);
  }

  if (!ok) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << "Failed to load " << filename << " on processor " << myprocid << std::endl;
    assert(ok);
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




Mesh_MSTK::~Mesh_MSTK() {
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

void Mesh_MSTK::clear_internals_ () 
{ 

  faceflip = NULL;

  cell_map_w_ghosts_ = cell_map_wo_ghosts_ = NULL;
  face_map_w_ghosts_ = face_map_wo_ghosts_ = NULL;
  node_map_w_ghosts_ = node_map_wo_ghosts_ = NULL;

  nmatsets = 0;
  nsidesets = 0;
  nnodesets = 0;
  matset_ids = sideset_ids = nodeset_ids = NULL;

}




void Mesh_MSTK::init_id_handle_maps() {
  int i, lid, nv, nf, nc, idx;
  MVertex_ptr vtx;
  MEntity_ptr genface;  // Mesh face in 3D, edge in 2D
  MEntity_ptr gencell;  // Mesh region in 3D, face in 2D

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
  while ((genface = MSet_Next_Entry(OwnedFaces,&idx))) {
    MEnt_Set_ID(genface,lid);
    face_id_to_handle[lid-1] = genface;
    lid++;
  }
  
  idx = 0;
  while ((genface = MSet_Next_Entry(NotOwnedFaces,&idx))) {
    MEnt_Set_ID(genface,lid);
    face_id_to_handle[lid-1] = genface;
    lid++;
  }
    


  nc = MSet_Num_Entries(AllCells);

  cell_id_to_handle.reserve(nc);

  idx = 0; lid = 1;
  while ((gencell = MSet_Next_Entry(OwnedCells,&idx))) {
    MEnt_Set_ID(gencell,lid);
    cell_id_to_handle[lid-1] = gencell;
    lid++;
  }
    
  idx = 0;
  while ((gencell = MSet_Next_Entry(GhostCells,&idx))) {
    MEnt_Set_ID(gencell,lid);
    cell_id_to_handle[lid-1] = gencell;
    lid++;
  }
    
}



void Mesh_MSTK::init_global_ids() {
    
}


void Mesh_MSTK::init_pvert_lists() {
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


void Mesh_MSTK::init_pface_lists() {
  int idx = 0;

  // Get all faces on this processor 

  if (cell_dimension() == 3) {

    MFace_ptr face;

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
  else if (cell_dimension() == 2) {

    MEdge_ptr edge;

    AllFaces = MSet_New(mesh,"AllFaces",MFACE);
    NotOwnedFaces = MSet_New(mesh,"NotOwnedFaces",MFACE);
    OwnedFaces = MSet_New(mesh,"OwnedFaces",MFACE);

    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {
      MSet_Add(AllFaces,edge);
      if (ME_PType(edge) == PGHOST)
	MSet_Add(NotOwnedFaces,edge);
      else
	MSet_Add(OwnedFaces,edge);
    }
  }
  else {
    std::cerr << "Not implemented for face dimension" << std::endl;
  }

  return;
}


void Mesh_MSTK::init_pface_dirs() {
  MRegion_ptr region0, region1;
  MFace_ptr face, face0, face1;
  MEdge_ptr edge;
  MAttrib_ptr attfc0, attfc1;
  int idx;
  int local_regid0, local_regid1;
  int remote_regid0, remote_regid1;
  int local_faceid0, local_faceid1;
  int remote_faceid0, remote_faceid1;


  if (serial_run) {
    faceflip = new bool[MSet_Num_Entries(AllFaces)];
    for (int i = 0; i < MSet_Num_Entries(AllFaces); i++) faceflip[i] = false;
  }
  else {
    // Do some additional processing to see if ghost faces and their masters
    // are oriented the same way; if not, turn on flag to flip the directions
    // when returning to the application code

    attfc0 = MAttrib_New(mesh,"TMP_FC0_ATT",INT,MFACE);
    attfc1 = MAttrib_New(mesh,"TMP_FC1_ATT",INT,MFACE);

    if (cell_dimension() == 3) {
    
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
    else if (cell_dimension() == 2) {


      idx = 0;
      while ((edge = MESH_Next_Edge(mesh,&idx))) {
	if (ME_PType(edge) != PINTERIOR) {
	  List_ptr efaces = ME_Faces(edge);
	    
	  face0 = List_Entry(efaces,0);
	  if (MF_EdgeDir(face0,edge) != 1) {
	    face1 = face0;
	    MEnt_Set_AttVal(edge,attfc1,MEnt_GlobalID(face1),0.0,NULL);

	    face0 = List_Entry(efaces,1);
	    if (face0) {
	      if (MF_EdgeDir(face0,edge) == -1)     // Sanity check
		MEnt_Set_AttVal(edge,attfc0,MEnt_GlobalID(face0),0.0,NULL);
	      else
		std::cerr << "Two faces using edge in same direction in 2D mesh" << std::endl;
	    }
	  }
	  else
	    MEnt_Set_AttVal(edge,attfc0,MEnt_GlobalID(face0),0.0,NULL);
	}
      }
    }



    MSTK_UpdateAttr(mesh, myprocid, numprocs, comm);



    faceflip = new bool[MSet_Num_Entries(AllFaces)];
    for (int i = 0; i < MSet_Num_Entries(AllFaces); i++) faceflip[i] = false;
    
    if (cell_dimension() == 3) {
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
    else if (cell_dimension() == 2) {
      double rval;
      void *pval;

      idx = 0;
      while ((edge = MSet_Next_Entry(NotOwnedFaces,&idx))) {
      
	MEnt_Get_AttVal(edge,attfc0,&remote_regid0,&rval,&pval);
	MEnt_Get_AttVal(edge,attfc1,&remote_regid1,&rval,&pval);
      
	List_ptr efaces = ME_Faces(edge);
	face0 = List_Entry(efaces,0);
	face1 = List_Entry(efaces,1);
	if (MF_EdgeDir(face0,edge) != 1) {
	  face0 = List_Entry(efaces,1);
	  face1 = List_Entry(efaces,0);
	}
	local_faceid0 = face0 ? MEnt_GlobalID(face0) : 0;
	local_faceid1 = face1 ? MEnt_GlobalID(face1) : 0;
      
	if (remote_faceid1 == local_faceid0 || 
	    remote_faceid0 == local_faceid1) {
	  int lid = MEnt_ID(edge);
	  faceflip[lid-1] = true;
	}
	else { // Sanity Check
	
	  if (remote_faceid1 != local_faceid1 &&
	      remote_faceid0 != local_faceid0) {
	  
	    cout << "Face cells mismatch between master and ghost (processor " << myprocid << ")" << std::endl;
	    cout << " Face " << MEnt_GlobalID(edge) << std::endl;
	    cout << "Remote cells " << remote_faceid0 << " " << remote_faceid1 << std::endl;
	    cout << "Local cells " << local_faceid0 << " " << local_faceid1 << std::endl;
	  }
	}
	List_Delete(efaces);
      }

    }
  }    
}


void Mesh_MSTK::init_pcell_lists() {
  int idx = 0;

  if (cell_dimension() == 3) {
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
  else if (cell_dimension() == 2) {
    MFace_ptr face;

    AllCells = MSet_New(mesh,"AllCells",MFACE);
    OwnedCells = MSet_New(mesh,"OwnedCells",MFACE);
    GhostCells = MSet_New(mesh,"GhostCells",MFACE);

    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      MSet_Add(AllCells,face);
      if (MF_PType(face) == PGHOST)
	MSet_Add(GhostCells,face);
      else
	MSet_Add(OwnedCells,face);
    }
  }
  else {
    std::cerr << "Implemented only for 2D and 3D" << std::endl;
  }

  return;
}



void Mesh_MSTK::init_set_info() {
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
      if (cell_dimension() == 3) { // create an empty set
	mset = MSet_New(mesh,setname.str().c_str(),MREGION);
      }
      else if (cell_dimension() == 2) {
	mset = MSet_New(mesh,setname.str().c_str(),MFACE);
      }
    }
  }
  

  for (i = 0; i < nsidesets; i++) {
    stringstream setname;
  
    setname << "sideset_" << sideset_ids[i];
    mset = MESH_MSetByName(mesh,setname.str().c_str());
    if (!mset) {
      if (cell_dimension() == 3) {
	mset = MSet_New(mesh,setname.str().c_str(),MFACE);
      }
      else if (cell_dimension() == 2) {
	mset = MSet_New(mesh,setname.str().c_str(),MFACE);
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

    
// Number of entities of any kind (cell, face, node) and in a
// particular category (OWNED, GHOST, USED)
    
unsigned int Mesh_MSTK::num_entities (const Entity_kind kind, 
				      const Parallel_type ptype) const
{
  const int rank = (int) kind;
  const int index = ((int) ptype) - 1;

  switch (kind) {
  case NODE:

    switch (ptype) {
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


  case FACE:
    switch (ptype) {
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


  case CELL:

    switch (ptype) {
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



// Get faces of a cell.

// On a distributed mesh, this will return all the faces of the
// cell OWNED or GHOST. The faces will be returned in a standard
// order according to Exodus II convention.
    

void Mesh_MSTK::cell_get_faces (const Entity_ID cellid, 
				std::vector<Entity_ID> *faceids)
{
  MEntity_ptr cell;

  int nfstd[4] = {4,5,5,6}; // expected number of faces for standard 3D cells
  int ncfstd[4][6] = { // Expected number of nodes for each face of 3D cell
    {3,3,3,3,0,0},         // Tet
    {4,4,4,3,3,0},         // Prism
    {3,3,3,3,4,0},         // Pyramid
    {4,4,4,4,4,4}          // Hex
  };
  int cfstd[4][6][4] = { // Expected cell-face-node pattern in 3D
    // Tet
    {{0,1,3,-1},{1,2,3,-1},{2,0,3,-1},{0,2,1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}},
    // Triangular prism
    {{0,1,4,3},{1,2,5,4},{2,0,3,5},{0,2,1,-1},{3,4,5,-1}},
    // Pyramid
    {{0,1,4,-1},{1,2,4,-1},{2,3,4,-1},{3,0,4,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}},
    // Hex
    {{0,1,5,4},{1,2,6,5},{2,3,7,6},{0,4,7,3},{0,3,2,1},{4,5,6,7}}
  };

  assert(faceids != NULL);

  faceids->clear();

  cell = cell_id_to_handle[cellid];
      
  if (cell_dimension() == 3) {
    int it;
    
    List_ptr rverts = MR_Vertices((MRegion_ptr)cell);
    List_ptr rfaces = MR_Faces((MRegion_ptr)cell);
    
    int nv = List_Num_Entries(rverts);
    int nf = List_Num_Entries(rfaces);
    
    if (nf == 4)                 // Tet
      it = 0;
    else if (nf == 5 && nv == 6) 
      it = 1;                    // prism
    else if (nf == 5 && nv == 5) 
      it = 2;                    // pyramid
    else if (nf == 6 && nv == 8) 
      it = 3;                    // Hex
    else
      it = -1;

    if (it >= 0) {

      // Have to re-sort the faces according a specific template for 
      // standard element types
      
      for (int i = 0; i < nfstd[it]; i++) {
	
	// Search for a face that has all the expected nodes
	
	bool found = false;
	
	for (int j = 0; j < nfstd[it]; j++) {
	  
	  bool all_present = true;
	  MFace_ptr face = List_Entry(rfaces,j);	  
	  List_ptr fverts = MF_Vertices(face,1,0);
	  
	  // Check if this face has all the expected nodes
	  
	  for (int k = 0; k < ncfstd[it][i]; k++) {
	    MVertex_ptr vtx = List_Entry(rverts,cfstd[it][i][k]);
	    
	    if (!List_Contains(fverts,vtx)) {
	      all_present = false;
	      break;
	    }
	  }
	  
	  if (all_present) {
	    int lid = MEnt_ID(face);
	    faceids->push_back(lid-1);
	    found = true;
	    break;
	  }
	  
	  List_Delete(fverts);
	  
	} // for (int j = 0; j < nfstd[it]; j++) 
	
	if (!found) {
	  std::cerr << "Could not find face in element" << std::endl;
	  assert(found);
	}	   
	
      } // for (int i = 0; i < nfstd[it]; i++) 
      
    }
    else {
      for (int j = 0; j < nf; j++) {
	MFace_ptr face = List_Entry(rfaces,j);
	int lid = MEnt_ID(face);
	faceids->push_back(lid-1);
      }
    }

    List_Delete(rfaces);
    List_Delete(rverts);
  }
  else {  // cell_dimension() = 2; surface or 2D mesh

    List_ptr fedges = MF_Edges((MFace_ptr)cellid,1,0);
    MEdge_ptr edge;
    int idx = 0;
    while ((edge = List_Next_Entry(fedges,&idx))) {
      int lid = MEnt_ID(edge);
      faceids->push_back(lid-1);
    }

  }
  
} // cell_get_faces



// Get directions in which a cell uses face
// In 3D, direction is 1 if face normal points out of cell
// and -1 if face normal points into cell
// In 2D, direction is 1 if face/edge is defined in the same
// direction as the cell polygon, and -1 otherwise
    
void Mesh_MSTK::cell_get_face_dirs (const Entity_ID cellid, 
				    std::vector<int> *face_dirs) 
{
  MEntity_ptr cell;
  std::vector<Entity_ID> faceids;
  int j,nf, result;

  assert(face_dirs != NULL);

  face_dirs->clear();
  
  cell = cell_id_to_handle[cellid];



  if (cell_dimension() == 3) {

    cell_get_faces(cellid,&faceids);
    nf = faceids.size();
   
    for (int i = 0; i < nf; i++) {
      MFace_ptr face = face_id_to_handle[faceids[i]];

      int cell_facedir = MR_FaceDir(cell,face) == 1 ? 1 : -1;

      // If this is a ghost face and the master has the opposite direction
      // we are supposed to flip it

      if (faceflip[faceids[i]]) cell_facedir *= -1;
      face_dirs->push_back(cell_facedir);
    }
    
  }
  else {
    cell_get_faces(cellid,&faceids);
    int ne = faceids.size();
   
    for (int i = 0; i < ne; i++) {
      MEdge_ptr edge = face_id_to_handle[faceids[i]];

      int cell_facedir = MF_EdgeDir(cell,edge) == 1 ? 1 : -1;

      // If this is a ghost face and the master has the opposite direction
      // we are supposed to flip it

      if (faceflip[faceids[i]]) cell_facedir *= -1;
      face_dirs->push_back(cell_facedir);
    }
    

  }
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

void Mesh_MSTK::cell_get_nodes (const Entity_ID cellid, 
				std::vector<Entity_ID> *nodeids)
{
  MEntity_ptr cell;
  int nn, lid;

  assert(nodeids != NULL);

  nodeids->clear();

  cell = cell_id_to_handle[cellid];
      
  if (cell_dimension() == 3) {                    // Volume mesh
    List_ptr rverts = MR_Vertices(cell);
 
    nn = List_Num_Entries(rverts);
    
    for (int i = 0; i < nn; i++) {
      lid = MEnt_ID(List_Entry(rverts,i));
      nodeids->push_back(lid-1);
    }
    
    List_Delete(rverts);
  }
  else {                                 // Surface mesh
    List_ptr fverts = MF_Vertices(cell,1,0);

    nn = List_Num_Entries(fverts);
    
    for (int i = 0; i < nn; i++) {
      lid = MEnt_ID(List_Entry(fverts,i));
      nodeids->push_back(lid-1);
    }
    
    List_Delete(fverts);
  }

}




    
// Get nodes of face 
// On a distributed mesh, all nodes (OWNED or GHOST) of the face 
// are returned
// In 3D, the nodes of the face are returned in ccw order consistent
// with the face normal
// In 2D, nfnodes is 2

void Mesh_MSTK::face_get_nodes (const Entity_ID faceid, 
				std::vector<Entity_ID> *nodeids)
{
  MEntity_ptr genface;
  int nn, lid;

  assert(nodeids != NULL);

  nodeids->clear();

  genface = face_id_to_handle[faceid];

  if (cell_dimension() == 3) {   // Volume mesh
    
    List_ptr fverts = MF_Vertices(genface,1,0);
    assert(fverts != NULL);

    nn = List_Num_Entries(fverts);
    
    if (faceflip[faceid]) {
      for (int i = nn-1; i >= 0; i--) {
	lid = MEnt_ID(List_Entry(fverts,i));
	nodeids->push_back(lid-1);
      }
    }  
    else {
      for (int i = 0; i < nn; i++) {
	lid = MEnt_ID(List_Entry(fverts,i));
	nodeids->push_back(lid-1);
      }
    }

    List_Delete(fverts);
  }
  else {                // Surface mesh or 2D mesh

    if (faceflip[faceid]) {
      nodeids->push_back(MEnt_ID(ME_Vertex(genface,1)));
      nodeids->push_back(MEnt_ID(ME_Vertex(genface,0)));
    }
    else {
      nodeids->push_back(MEnt_ID(ME_Vertex(genface,0)));
      nodeids->push_back(MEnt_ID(ME_Vertex(genface,1)));
    }
  }

}
  




// Cells of type 'ptype' connected to a node
    
void Mesh_MSTK::node_get_cells (const Entity_ID nodeid, 
				const Parallel_type ptype,
				std::vector<Entity_ID> *cellids)
{
  int idx, lid;
  List_ptr cell_list;
  MEntity_ptr ment;

  assert (cellids != NULL);

  cellids->clear();


  MVertex_ptr mv = (MVertex_ptr) vtx_id_to_handle[nodeid];

  if (cell_dimension() == 3)
    cell_list = MV_Regions(mv);
  else
    cell_list = MV_Faces(mv);

  idx = 0;
  while ((ment = List_Next_Entry(cell_list,&idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == GHOST || ptype == USED) {
	lid = MEnt_ID(ment);
	cellids->push_back(lid-1);
      }
    }
    else {
      if (ptype == OWNED || ptype == USED) {
	lid = MEnt_ID(ment);
	cellids->push_back(lid-1);
      }
    }      
  }

  List_Delete(cell_list);
}
    
  // Faces of type 'ptype' connected to a node
    
void Mesh_MSTK::node_get_faces (const Entity_ID nodeid, 
				const Parallel_type ptype,
				std::vector<Entity_ID> *faceids)
{
  int idx, lid;
  List_ptr face_list;
  MEntity_ptr ment;

  assert(faceids != NULL);

  faceids->clear();

  MVertex_ptr mv = (MVertex_ptr) vtx_id_to_handle[nodeid];

  if (cell_dimension() == 3)
    face_list = MV_Faces(mv);
  else
    face_list = MV_Edges(mv);

  idx = 0;
  while ((ment = List_Next_Entry(face_list,&idx))) {
    if (MEnt_PType(ment) == PGHOST) {
      if (ptype == GHOST || ptype == USED) {
	lid = MEnt_ID(ment);
	faceids->push_back(lid-1);
      }
    }
    else {
      if (ptype == OWNED || ptype == USED) {
	lid = MEnt_ID(ment);
	faceids->push_back(lid-1);
      }
    }      
  }

  List_Delete(face_list);
}



    
// Get faces of ptype of a particular cell that are connected to the
// given node
    
void Mesh_MSTK::node_get_cell_faces (const Entity_ID nodeid, 
				     const Entity_ID cellid,
				     const Parallel_type ptype,
				     std::vector<Entity_ID> *faceids)
{
  int idx, lid;
  List_ptr cell_list;
  MEntity_ptr ment;
  MRegion_ptr mr;
  MFace_ptr mf;
  MEdge_ptr me;

  assert(faceids != NULL);

  faceids->clear();

  MVertex_ptr mv = (MVertex_ptr) vtx_id_to_handle[nodeid];

  if (cell_dimension() == 3) {
    mr = (MFace_ptr) cell_id_to_handle[cellid];
    List_ptr rfaces = MR_Faces(mr);
    idx = 0;
    while ((mf = List_Next_Entry(rfaces,&idx))) {
      if (!MF_UsesEntity(mf,mv,MVERTEX)) continue;

      if (MEnt_PType(mf) == PGHOST) {
	if (ptype == GHOST || ptype == USED) {
	  lid = MEnt_ID(ment);
	  faceids->push_back(lid-1);
	}
      }
      else {
	if (ptype == OWNED || ptype == USED) {
	  lid = MEnt_ID(ment);
	  faceids->push_back(lid-1);
	}
      }            
    }
    List_Delete(rfaces);
  }
  else {
    mf = (MEdge_ptr) cell_id_to_handle[cellid];
    List_ptr fedges = MF_Edges(mr,1,0);
    idx = 0;
    while ((me = List_Next_Entry(fedges,&idx))) {
      if (!ME_UsesEntity(me,mv,MVERTEX)) continue;

      if (MEnt_PType(me) == PGHOST) {
	if (ptype == GHOST || ptype == USED) {
	  lid = MEnt_ID(ment);
	  faceids->push_back(lid-1);
	}
      }
      else {
	if (ptype == OWNED || ptype == USED) {
	  lid = MEnt_ID(ment);
	  faceids->push_back(lid-1);
	}
      }            
    }
    List_Delete(fedges);
  }

}    


    
// Cells connected to a face
    
void Mesh_MSTK::face_get_cells (const Entity_ID faceid, 
				const Parallel_type ptype,
				std::vector<Entity_ID> *cellids)
{
  int lid;

  assert(cellids != NULL);

  cellids->clear();

  if (cell_dimension() == 3) {
    MFace_ptr mf = (MFace_ptr) face_id_to_handle[faceid];
   
    List_ptr fregs = MF_Regions(mf);
    MRegion_ptr mr;
    int idx = 0;
    while ((mr = List_Next_Entry(fregs,&idx))) {
      if (MEnt_PType(mr) == PGHOST) {
	if (ptype == GHOST || ptype == USED) {
	  lid = MR_ID(mr);
	  cellids->push_back(lid-1);
	}
      }
      else {
	if (ptype == OWNED || ptype == USED) {
	  lid = MR_ID(mr);
	  cellids->push_back(lid-1);
	}
      }
    }
    List_Delete(fregs);

  }
  else {
    MEdge_ptr me = (MEdge_ptr) face_id_to_handle[faceid];

    List_ptr efaces = ME_Faces(me);
    MFace_ptr mf;
    int idx = 0;
    while ((mf = List_Next_Entry(efaces,&idx))) {
      if (MEnt_PType(mf) == PGHOST) {
	if (ptype == GHOST || ptype == USED) {
	  lid = MF_ID(mf);
	  cellids->push_back(lid-1);
	}
      }
      else {
	if (ptype == OWNED || ptype == USED) {
	  lid = MF_ID(mf);
	  cellids->push_back(lid-1);
	}
      }
    }
    List_Delete(efaces);

  }

}
    


// Same level adjacencies
//-----------------------

// Face connected neighboring cells of given cell

void Mesh_MSTK::cell_get_face_adj_cells(const Entity_ID cellid,
					const Parallel_type ptype,
					std::vector<Entity_ID> *fadj_cellids)
{
  int lid;

  assert(fadj_cellids != NULL);

  fadj_cellids->clear();

  if (cell_dimension() == 3) {

    MRegion_ptr mr = (MRegion_ptr) cell_id_to_handle[cellid];

    List_ptr rfaces = MR_Faces(mr);
    int idx = 0;
    MFace_ptr mf;
    while ((mf = List_Next_Entry(rfaces,&idx))) {
      List_ptr fregs = MF_Regions(mf);
      int idx2 = 0;
      MRegion_ptr mr2;
      while ((mr2 = List_Next_Entry(fregs,&idx2))) {
	if (mr2 != mr) {
	  if (MEnt_PType(mr2) == PGHOST) {
	    if (ptype == GHOST || ptype == USED) {
	      lid = MEnt_ID(mr2);
	      fadj_cellids->push_back(lid-1);
	    }
	  }
	  else {
	    if (ptype == GHOST || ptype == USED) {
	      lid = MEnt_ID(mr2);
	      fadj_cellids->push_back(lid-1);
	    }
	  }
	}
      }
      List_Delete(fregs);
    }

    List_Delete(rfaces);

  }
  else if (cell_dimension() == 2) {

    MFace_ptr mf = (MFace_ptr) cell_id_to_handle[cellid];

    List_ptr fedges = MF_Edges(mf,1,0);
    int idx = 0;
    MEdge_ptr me;
    while ((me = List_Next_Entry(fedges,&idx))) {
      List_ptr efaces = ME_Faces(me);
      int idx2 = 0;
      MFace_ptr mf2;
      while ((mf2 = List_Next_Entry(fedges,&idx2))) {
	if (mf2 != mf) {
	  if (MEnt_PType(mf2) == PGHOST) {
	    if (ptype == GHOST || ptype == USED) {
	      lid = MEnt_ID(mf2);
	      fadj_cellids->push_back(lid-1);
	    }
	  }
	  else {
	    if (ptype == GHOST || ptype == USED) {
	      lid = MEnt_ID(mf2);
	      fadj_cellids->push_back(lid-1);
	    }
	  }
	}
      }
      List_Delete(efaces);
    }

    List_Delete(fedges);

  }

}




// Node connected neighboring cells of given cell 

void Mesh_MSTK::cell_get_node_adj_cells(const Entity_ID cellid,
					const Parallel_type ptype,
					std::vector<Entity_ID> *nadj_cellids)
{
  int lid, mkid;
  List_ptr cell_list;

  assert(nadj_cellids != NULL);

  nadj_cellids->clear();

  mkid = MSTK_GetMarker();

  cell_list = List_New(0);
  if (cell_dimension() == 3) {

    MRegion_ptr mr = (MRegion_ptr) cell_id_to_handle[cellid];

    List_ptr rvertices = MR_Vertices(mr);
    int idx = 0;
    MVertex_ptr mv;
    while ((mv = List_Next_Entry(rvertices,&idx))) {
      List_ptr vregs = MV_Regions(mv);
      int idx2 = 0;
      MRegion_ptr mr2;
      while ((mr2 = List_Next_Entry(vregs,&idx2))) {
	if (mr2 != mr && !MEnt_IsMarked(mr2,mkid)) {
	  MEnt_Mark(mr2,mkid);
	  List_Add(cell_list,mr2);
	  if (MEnt_PType(mr2) == PGHOST) {
	    if (ptype == GHOST || ptype == USED) {
	      lid = MEnt_ID(mr2);
	      nadj_cellids->push_back(lid-1);
	    }
	  }
	  else {
	    if (ptype == GHOST || ptype == USED) {
	      lid = MEnt_ID(mr2);
	      nadj_cellids->push_back(lid-1);
	    }
	  }
	}
      }
      List_Delete(vregs);
    }

    List_Delete(rvertices);
    
  }
  else if (cell_dimension() == 2) {

    MFace_ptr mf = (MFace_ptr) cell_id_to_handle[cellid];

    List_ptr fverts = MF_Vertices(mf,1,0);
    int idx = 0;
    MVertex_ptr mv;
    while ((mv = List_Next_Entry(fverts,&idx))) {
      List_ptr vfaces = MV_Faces(mv);
      int idx2 = 0;
      MFace_ptr mf2;
      while ((mf2 = List_Next_Entry(vfaces,&idx2))) {
	if (mf2 != mf && !MEnt_IsMarked(mf2,mkid)) {
	  MEnt_Mark(mf2,mkid);
	  List_Add(cell_list,mf2);
	  if (MEnt_PType(mf2) == PGHOST) {
	    if (ptype == GHOST || ptype == USED) {
	      lid = MEnt_ID(mf2);
	      nadj_cellids->push_back(lid-1);
	    }
	  }
	  else {
	    if (ptype == GHOST || ptype == USED) {
	      lid = MEnt_ID(mf2);
	      nadj_cellids->push_back(lid-1);
	    }
	  }
	}
      }
      List_Delete(vfaces);
    }

    List_Delete(fverts);

  }

  List_Unmark(cell_list,mkid);
  List_Delete(cell_list);
  MSTK_FreeMarker(mkid);

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
    
Cell_type Mesh_MSTK::cell_get_type_4viz(const Entity_ID cellid) const
{

}
    
    
// See cell_get_nodes for details on node ordering
    
void Mesh_MSTK::cell_get_nodes_4viz (const Entity_ID cellid, 
					 std::vector<Entity_ID> *nodeids)
{

}
    
    



// Node coordinates - 3 in 3D and 2 in 2D
    
void Mesh_MSTK::node_get_coordinates (const Entity_ID nodeid, Point *ncoords)
{    
  MEntity_ptr vtx;
  double coords[3];

  assert(ncoords != NULL);

  vtx = vtx_id_to_handle[nodeid];

  MV_Coords(vtx,coords);

  if (space_dimension() == 3) 
    ncoords->set(coords[0],coords[1],coords[2]);
  else if (space_dimension() == 2)
    ncoords->set(coords[0],coords[1]);
}



// Coordinates of cells in standard order (Exodus II convention)
// STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
// For a general polyhedron this will return the node coordinates in
// arbitrary order
// Number of nodes is vector size divided by number of spatial dimensions

void Mesh_MSTK::cell_get_coordinates (const Entity_ID cellid, std::vector<Point> *ccoords)
{    
  MEntity_ptr cell;
  double coords[3];
  int nn, result;

  assert(ccoords != NULL);

  ccoords->clear();

  cell = cell_id_to_handle[cellid];
      
  if (cell_dimension() == 3) {
    List_ptr rverts = MR_Vertices(cell);
    
    nn = List_Num_Entries(rverts);

    for (int i = 0; i < nn; i++) {
      MV_Coords(List_Entry(rverts,i),coords);

      ((*ccoords)[i]).init(3);
      ((*ccoords)[i]).set(coords[0],coords[1],coords[2]);
    }    

    List_Delete(rverts);
  }
  else {
    List_ptr fverts = MF_Vertices(cell,1,0);
    nn = List_Num_Entries(fverts);

    for (int i = 0; i < nn; i++) {
      MV_Coords(List_Entry(fverts,i),coords);

      ((*ccoords)[i]).init(2);
      ((*ccoords)[i]).set(coords[0],coords[1],coords[2]);
    }
  }
}




// Face coordinates - conventions same as face_to_nodes call 
// Number of nodes is the vector size divided by number of spatial dimensions
    
void Mesh_MSTK::face_get_coordinates (const Entity_ID faceid, std::vector<Point> *fcoords)
{
  double coords[3];
  int nn;

  assert(fcoords != NULL);

  fcoords->clear();

  if (space_dimension() == 3) {
    if (cell_dimension() == 3) {
      MFace_ptr face = face_id_to_handle[faceid];
      List_ptr fverts = MF_Vertices(face,1,0);

      nn = List_Num_Entries(fverts);
	
      if (faceflip[faceid]) {
	for (int i = nn-1; i >=0; i--) {
	  MV_Coords(List_Entry(fverts,i),coords);
	  ((*fcoords)[nn-i-1]).init(3);
	  ((*fcoords)[nn-i-i]).set(coords[0],coords[1],coords[2]);
	}
      }
      else {
	for (int i = 0; i < nn; i++) {
	  MV_Coords(List_Entry(fverts,i),coords);
	  ((*fcoords)[i]).init(3);
	  ((*fcoords)[i]).set(coords[0],coords[1],coords[2]);
	}	  
      }
    }
    else { // Surface mesh embedded in 3D
      MEdge_ptr edge = face_id_to_handle[faceid];
      MVertex_ptr ev[2];
      if (!faceflip[faceid]) {
	ev[0] = ME_Vertex(edge,0);
	ev[1] = ME_Vertex(edge,1);
      }
      else {
	ev[1] = ME_Vertex(edge,0);
	ev[0] = ME_Vertex(edge,1);
      }

      MV_Coords(ev[0],coords);
      ((*fcoords)[0]).init(3);
      ((*fcoords)[0]).set(coords[0],coords[1],coords[2]);
	
      MV_Coords(ev[1],coords);
      ((*fcoords)[1]).init(3);
      ((*fcoords)[1]).set(coords[0],coords[1],coords[2]);
    }
  }
  else {  // 2D Mesh
    MEdge_ptr edge;
    MVertex_ptr ev[2];

    edge = face_id_to_handle[faceid];
    if (!faceflip[faceid]) {
      ev[0] = ME_Vertex(edge,0);
      ev[1] = ME_Vertex(edge,1);
    }
    else {
      ev[1] = ME_Vertex(edge,0);
      ev[0] = ME_Vertex(edge,1);
    }
      
    MV_Coords(ev[0],coords);
    ((*fcoords)[0]).init(2);
    ((*fcoords)[0]).set(coords[0],coords[1]);
      
    MV_Coords(ev[1],coords);
    ((*fcoords)[1]).init(2);
    ((*fcoords)[1]).set(coords[0],coords[1]);
  }	  
}
  


  // Ids of sets containing entities of 'kind'

void Mesh_MSTK::get_set_ids (const Entity_kind kind, std::vector<Set_ID> *setids) 
{
  int i, ns=0;

  assert(setids != NULL);

  setids->clear();

  switch (kind) {
  case CELL: {
    int sids[nmatsets];

    for (i = 0; i < nmatsets; i++)
      setids->push_back(matset_ids[i]);

    return;
    break;
  }

  case FACE: {
    int sids[nsidesets];

    for (i = 0; i < nsidesets; i++)
      setids->push_back(sideset_ids[i]);

    return;
    break;
  }

  case NODE: {
    int sids[nnodesets];

    for (i = 0; i < nnodesets; i++)
      setids->push_back(nodeset_ids[i]);

    return;
    break;
  }

  default:
    return;
  }
    
}



// Get list of entities of type 'category' in set

void Mesh_MSTK::get_set_entities (const Set_ID setid, 
				  const Entity_kind kind, 
				  const Parallel_type ptype, 
				  std::vector<Entity_ID> *setents) 
{
  int idx, i, lid;
  MSet_ptr mset, mset1;
  MEntity_ptr ment;
  stringstream setname;


  assert(setents != NULL);
  
  setents->clear();

  switch (kind) {
  case CELL: {
    int ncells;
    int *cells;

    setname << "matset_" << setid;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    ncells = MSet_Num_Entries(mset1);
    if (ncells == 0) return;

    if (ptype == OWNED)
      mset = MSets_Subtract(mset1,GhostCells);
    else if (ptype == GHOST)
      mset = MSets_Subtract(mset1,OwnedCells);
    else
      mset = mset1;

    if (!mset) return;

    ncells = MSet_Num_Entries(mset);

    idx = 0;
    while ((ment = MSet_Next_Entry(mset,&idx))) {
      lid = MEnt_ID(ment);
      setents->push_back(lid-1);
    }

    if (ptype != USED)
      MSet_Delete(mset);

    return;
    break;
  }
  case FACE: {

    setname << "sideset_" << setid;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    if (!MSet_Num_Entries(mset1)) return;

    if (ptype == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedFaces);
    else if (ptype == GHOST)
      mset = MSets_Subtract(mset1,OwnedFaces);
    else
      mset = mset1;

    if (!mset) return;

    idx = 0;
    while ((ment = MSet_Next_Entry(mset,&idx))) {
      lid = MEnt_ID(ment);
      setents->push_back(lid-1);
    }

    if (ptype != USED)
      MSet_Delete(mset);

    return;
    break;
  }
  case NODE: {

    setname << "nodeset_" << setid;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    if (!MSet_Num_Entries(mset1)) return;
    
    if (ptype == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedVerts);
    else if (ptype == GHOST)
      mset = MSets_Subtract(mset1,OwnedVerts);
    else
      mset = mset1;

    if (!mset) return;

    idx = 0;
    while ((ment = MSet_Next_Entry(mset,&idx))) {
      lid = MEnt_ID(ment);
      setents->push_back(lid-1);
    }

    if (ptype != USED)
      MSet_Delete(mset);

    return;
    break;
  }
  default:
    return;
  }
  
}


// Number of sets containing entities of type 'kind' in mesh
    
unsigned int Mesh_MSTK::num_sets(const Entity_kind kind) const
{    
  int n = 0;
    
  switch (kind) {
  case CELL:
    return nmatsets;
  case FACE:
    return nsidesets;
  case NODE:
    return nnodesets;
  default:
    return 0;
  }
}


  // Is this is a valid ID of a set containing entities of 'kind'

bool Mesh_MSTK::valid_set_id (const Set_ID setid, const Entity_kind kind) const 
{
  int n = 0;

  switch (kind) {
  case CELL:
    for (int i = 0; i < nmatsets; i++)
      if (matset_ids[i] == setid) return true;
    break;
  case FACE:
    for (int i = 0; i < nsidesets; i++)
      if (sideset_ids[i] == setid) return true;
    break;
  case NODE:
    for (int i = 0; i < nnodesets; i++)
      if (nodeset_ids[i] == setid) return true;
    break;
  default:
    return false;
  }    

  return false;
}




  // Get number of entities of type 'ptype' in set

unsigned int Mesh_MSTK::get_set_size (const Set_ID setid, 
				      const Entity_kind kind,
				      const Parallel_type ptype)
{
  stringstream setname;
  MSet_ptr mset1, mset;
  int nent;

  switch (kind) {
  case CELL: {

    setname << "matset_" << setid;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    if (!MSet_Num_Entries(mset1)) return 0;

    if (ptype == OWNED)
      mset = MSets_Subtract(mset1,GhostCells);
    else if (ptype == GHOST)
      mset = MSets_Subtract(mset1,OwnedCells);
    else
      mset = mset1;

    if (!mset) return 0;

    nent = MSet_Num_Entries(mset);
    if (ptype != USED)
      MSet_Delete(mset);

    return nent;
    break;
  }

  case FACE: {

    setname << "sideset_" << setid;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    if (!MSet_Num_Entries(mset1)) return 0;

    if (ptype == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedFaces);
    else if (ptype == GHOST)
      mset = MSets_Subtract(mset1,OwnedFaces);
    else
      mset = mset1;
    
    if (!mset) return 0;

    nent = MSet_Num_Entries(mset);
    if (ptype != USED)
      MSet_Delete(mset);

    return nent;
    break;
  }

  case NODE: {

    setname << "nodeset_" << setid;
    mset1 = MESH_MSetByName(mesh,setname.str().c_str());
    assert (mset1 != NULL);

    if (!MSet_Num_Entries(mset1)) return 0;

    if (ptype == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedVerts);
    else if (ptype == GHOST)
      mset = MSets_Subtract(mset1,OwnedVerts);
    else
      mset = mset1;

    if (!mset) return 0;

    nent = MSet_Num_Entries(mset);
    if (ptype != USED)
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

void Mesh_MSTK::init_cell_map ()
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

void Mesh_MSTK::init_face_map ()
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

void Mesh_MSTK::init_node_map ()
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


// Global ID of any entity

unsigned int Mesh_MSTK::GID(const Entity_ID lid, const Entity_kind kind) const
{
  MEntity_ptr ent;
  unsigned int gid;

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

  if (serial_run)
    return MEnt_ID(ent)-1;
  else
    return MEnt_GlobalID(ent)-1;
}

