#include "Mesh_MSTK.hh"
//#include <Teuchos_RCP.hpp>

using namespace std;


namespace Amanzi
{

namespace AmanziMesh
{


Mesh_MSTK::Mesh_MSTK (const char *filename, MPI_Comm incomm,
		      const AmanziGeometry::GeometricModelPtr& gm)
{
  // Assume three dimensional problem if constructor called without 
  // the space_dimension parameter

  int space_dimension = 3;

  Mesh_MSTK(filename,incomm,space_dimension,gm);

}

// Constructor - load up mesh from file

Mesh_MSTK::Mesh_MSTK (const char *filename, MPI_Comm incomm, int space_dimension,
		      const AmanziGeometry::GeometricModelPtr& gm)
{
  int ok;
  int ring = 1; // One layer of ghost cells in parallel meshes
  int with_attr = 1;  // update of attributes in parallel meshes

  clear_internals_();

  MSTK_Init();

  Mesh::set_geometric_model(gm);

  set_space_dimension(space_dimension);

  mpicomm = incomm;
  MPI_Comm_rank(mpicomm,&myprocid);
  MPI_Comm_size(mpicomm,&numprocs);

  serial_run =  (!mpicomm || numprocs == 1) ? true : false;

  Mesh::set_comm(mpicomm);

    if (myprocid == 0) {
    int DebugWait=0;
    while (DebugWait);
    }

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

    ok = ok & MSTK_Mesh_Distribute(&mesh,dim,ring,with_attr,myprocid,
					numprocs,mpicomm);

    if (myprocid == 0)
      MESH_Delete(globalmesh);

    set_cell_dimension(dim);
  }

  if (!ok) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << "Failed to load " << filename << " on processor " << myprocid << std::endl;
    assert(ok);
  }


  // Pre-process the mesh to remove degenerate edges
  
  collapse_degen_edges();

  label_celltype();




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
  
  if (gm != NULL)
    init_set_info();

}




Mesh_MSTK::~Mesh_MSTK() {
  delete cell_map_wo_ghosts_;
  delete cell_map_w_ghosts_;
  delete face_map_wo_ghosts_;
  delete face_map_w_ghosts_;
  delete node_map_wo_ghosts_;
  delete node_map_w_ghosts_;
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

  std::vector<Entity_ID> *orig_vids;
  void *pval;
  MRegion_ptr region;
  int idx = 0;
  while ((region = MESH_Next_Region(mesh,&idx))) {
    MEnt_Get_AttVal(region,orig_celltopo_att,NULL,NULL,&pval);
    orig_vids = (std::vector<Entity_ID> *) pval;
    delete orig_vids;
  }

  MFace_ptr face;
  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx))) {
    MEnt_Get_AttVal(face,orig_celltopo_att,NULL,NULL,&pval);
    orig_vids = (std::vector<Entity_ID> *) pval;
    delete orig_vids;
  }
  MAttrib_Delete(orig_celltopo_att);

  MAttrib_Delete(orig_celltype_att);
    
  MAttrib_Delete(celltype_att);

  // MESH_Delete(mesh);
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


// Get cell type

Cell_type Mesh_MSTK::cell_get_type(const Entity_ID cellid) const {
  MEntity_ptr cell;
  int ival;
  Cell_type celltype;
  
  cell = cell_id_to_handle[cellid];
  
  MEnt_Get_AttVal(cell,celltype_att,&ival,NULL,NULL);
  celltype = (Cell_type) ival;

  return celltype;

} // Mesh_MSTK::cell_get_type



// Get faces of a cell.

// On a distributed mesh, this will return all the faces of the
// cell OWNED or GHOST. The faces will be returned in a standard
// order according to Exodus II convention.
    

void Mesh_MSTK::cell_get_faces (const Entity_ID cellid, 
				std::vector<Entity_ID> *faceids) const
{
  MEntity_ptr cell;

  assert(faceids != NULL);

  faceids->clear();

  cell = cell_id_to_handle[cellid];
      
  if (cell_dimension() == 3) {
    int celltype;

    celltype = cell_get_type(cellid);

    List_ptr rverts = MR_Vertices((MRegion_ptr)cell);
    List_ptr rfaces = MR_Faces((MRegion_ptr)cell);
    
    int nv = List_Num_Entries(rverts);
    int nf = List_Num_Entries(rfaces);
    
    if (celltype >= TET && celltype <= HEX) {

      // Have to re-sort the faces according a specific template for 
      // standard element types
      
      int nfstd = nface_std[celltype];

      for (int i = 0; i < nfstd; i++) {
	
	// Search for a face that has all the expected nodes
	
	bool found = false;
	
	for (int j = 0; j < nf; j++) {
	  
	  bool all_present = true;
	  MFace_ptr face = List_Entry(rfaces,j);	  
	  List_ptr fverts = MF_Vertices(face,1,0);
	  
	  // Check if this face has all the expected nodes
	  
	  int nfn = nfnodes_std[celltype][i];

	  for (int k = 0; k < nfn; k++) {
	    int nodenum = fnodes_std[celltype][i][k];
	    MVertex_ptr vtx = List_Entry(rverts,nodenum);
	    
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
	  
	} // for (int j = 0; j < nf; j++) 
	
	if (!found) {
	  std::cerr << "Could not find face in element" << std::endl;
	  assert(found);
	}	   
	
      } // for (int i = 0; i < nf; i++) 
      
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

    List_ptr fedges = MF_Edges((MFace_ptr)cell,1,0);
    MEdge_ptr edge;
    int idx = 0;
    while ((edge = List_Next_Entry(fedges,&idx))) {
      int lid = MEnt_ID(edge);
      faceids->push_back(lid-1);
    }

  }
  
} // Mesh_MSTK::cell_get_faces



// Get directions in which a cell uses face
// In 3D, direction is 1 if face normal points out of cell
// and -1 if face normal points into cell
// In 2D, direction is 1 if face/edge is defined in the same
// direction as the cell polygon, and -1 otherwise
    
void Mesh_MSTK::cell_get_face_dirs (const Entity_ID cellid, 
				    std::vector<int> *face_dirs) const 
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
} // Mesh_MSTK::cell_get_face_dirs



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
				std::vector<Entity_ID> *nodeids) const
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

}  // Mesh_MSTK::cell_get_nodes




    
// Get nodes of face 
// On a distributed mesh, all nodes (OWNED or GHOST) of the face 
// are returned
// In 3D, the nodes of the face are returned in ccw order consistent
// with the face normal
// In 2D, nfnodes is 2

void Mesh_MSTK::face_get_nodes (const Entity_ID faceid, 
				std::vector<Entity_ID> *nodeids) const
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

} // Mesh_MSTK::face_get_nodes
  




// Cells of type 'ptype' connected to a node
    
void Mesh_MSTK::node_get_cells (const Entity_ID nodeid, 
				const Parallel_type ptype,
				std::vector<Entity_ID> *cellids) const
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

} // Mesh_MSTK::node_get_cells


    
// Faces of type 'ptype' connected to a node
    
void Mesh_MSTK::node_get_faces (const Entity_ID nodeid, 
				const Parallel_type ptype,
				std::vector<Entity_ID> *faceids) const
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

} // Mesh_MSTK::node_get_faces



    
// Get faces of ptype of a particular cell that are connected to the
// given node
    
void Mesh_MSTK::node_get_cell_faces (const Entity_ID nodeid, 
				     const Entity_ID cellid,
				     const Parallel_type ptype,
				     std::vector<Entity_ID> *faceids) const
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

} // Mesh_MSTK::node_get_cell_faces


    
// Cells connected to a face
    
void Mesh_MSTK::face_get_cells (const Entity_ID faceid, 
				const Parallel_type ptype,
				std::vector<Entity_ID> *cellids) const
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

} // Mesh_MSTK::face_get_cells
    


// Same level adjacencies
//-----------------------

// Face connected neighboring cells of given cell

void Mesh_MSTK::cell_get_face_adj_cells(const Entity_ID cellid,
					const Parallel_type ptype,
					std::vector<Entity_ID> *fadj_cellids) const
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

} // Mesh_MSTK::cell_get_face_adj_cells




// Node connected neighboring cells of given cell 

void Mesh_MSTK::cell_get_node_adj_cells(const Entity_ID cellid,
					const Parallel_type ptype,
					std::vector<Entity_ID> *nadj_cellids) const
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

} // Mesh_MSTK::cell_get_node_adj_cells


    
//
// Mesh Topology for viz  
//----------------------
//
// We need a special function because certain types of degenerate
// hexes will not be recognized as any standard element type (hex,
// pyramid, prism or tet). The original topology of this element 
// without any collapsed nodes will be returned by this call.


// Get cell type for viz (will return original type if cell has been
// modified)

Cell_type Mesh_MSTK::cell_get_type_4viz(const Entity_ID cellid) const {
  MEntity_ptr cell;
  Cell_type celltype;
  int ival;
  
  cell = cell_id_to_handle[cellid];
  
  MEnt_Get_AttVal(cell,orig_celltype_att,&ival,NULL,NULL);
  celltype = (Cell_type) ival;

  if (celltype == 0)
    return cell_get_type(cellid);
  else
    return celltype;

} // Mesh_MSTK::cell_get_type_4viz


void Mesh_MSTK::cell_get_nodes_4viz (const Entity_ID cellid, 
				     std::vector<Entity_ID> *nodeids) const
{
  MEntity_ptr cell;
  std::vector<Entity_ID> *orig_vids;
  void *pval;

  assert(nodeids != NULL);

  nodeids->clear();

  cell = cell_id_to_handle[cellid];
      
  MEnt_Get_AttVal(cell,orig_celltopo_att,NULL,NULL,&pval);
  if (pval) {
    orig_vids = (std::vector<Entity_ID> *) pval;
    *nodeids = *orig_vids; // copy content from *orig_vids to *nodeids
  }
  else
    cell_get_nodes(cellid,nodeids);

} // Mesh_MSTK::cell_get_nodes_4viz





// Node coordinates - 3 in 3D and 2 in 2D
    
void Mesh_MSTK::node_get_coordinates (const Entity_ID nodeid, AmanziGeometry::Point *ncoords) const
{    
  MEntity_ptr vtx;
  double coords[3];

  assert(ncoords != NULL);

  vtx = vtx_id_to_handle[nodeid];

  MV_Coords(vtx,coords);

  if (space_dimension() == 3) {
    ncoords->init(3);
    ncoords->set(coords[0],coords[1],coords[2]);
  }
  else if (space_dimension() == 2) {
    ncoords->init(2);
    ncoords->set(coords[0],coords[1]);
  }

} // Mesh_MSTK::node_get_coordinates


// Coordinates of cells in standard order (Exodus II convention)
// STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
// For a general polyhedron this will return the node coordinates in
// arbitrary order
// Number of nodes is vector size divided by number of spatial dimensions

void Mesh_MSTK::cell_get_coordinates (const Entity_ID cellid, std::vector<AmanziGeometry::Point> *ccoords) const
{    
  MEntity_ptr cell;
  double coords[3];
  int nn, result;
  int spdim = space_dimension(), celldim = cell_dimension();
  AmanziGeometry::Point xyz(spdim);

  assert(ccoords != NULL);

  ccoords->clear();

  cell = cell_id_to_handle[cellid];
      
  if (celldim == 3) {
    List_ptr rverts = MR_Vertices(cell);
    
    nn = List_Num_Entries(rverts);

    for (int i = 0; i < nn; i++) {
      MV_Coords(List_Entry(rverts,i),coords);

      xyz.set(coords[0],coords[1],coords[2]);
      ccoords->push_back(xyz);
    }    

    List_Delete(rverts);
  }
  else if (celldim == 2) {
    List_ptr fverts = MF_Vertices(cell,1,0);
    nn = List_Num_Entries(fverts);

    for (int i = 0; i < nn; i++) {
      MV_Coords(List_Entry(fverts,i),coords);

      if (spdim == 2)
	xyz.set(coords[0],coords[1]);
      else
	xyz.set(coords[0],coords[1],coords[2]);

      ccoords->push_back(xyz);
    }
  }

} // Mesh_MSTK::cell_get_coordinates_4viz




// Face coordinates - conventions same as face_to_nodes call 
// Number of nodes is the vector size divided by number of spatial dimensions
    
void Mesh_MSTK::face_get_coordinates (const Entity_ID faceid, std::vector<AmanziGeometry::Point> *fcoords) const
{
  double coords[3];
  int nn;
  int spdim = space_dimension(), celldim = cell_dimension();

  assert(fcoords != NULL);

  fcoords->clear();

  if (spdim == 3) {
    AmanziGeometry::Point xyz(3);
    if (celldim == 3) {
      MFace_ptr face = face_id_to_handle[faceid];
      List_ptr fverts = MF_Vertices(face,1,0);

      nn = List_Num_Entries(fverts);
	
      if (faceflip[faceid]) {
	for (int i = nn-1; i >=0; i--) {
	  MV_Coords(List_Entry(fverts,i),coords);

	  xyz.set(coords[0],coords[1],coords[2]);
	  fcoords->push_back(xyz);
	}
      }
      else {
	for (int i = 0; i < nn; i++) {
	  MV_Coords(List_Entry(fverts,i),coords);

	  xyz.set(coords[0],coords[1],coords[2]);
	  fcoords->push_back(xyz);
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
      xyz.set(coords[0],coords[1],coords[2]);
      fcoords->push_back(xyz);
	
      MV_Coords(ev[1],coords);
      xyz.set(coords[0],coords[1],coords[2]);
      fcoords->push_back(xyz);
    }
  }
  else {  // 2D Mesh
    AmanziGeometry::Point xyz(2);
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
    xyz.set(coords[0],coords[1]);
    fcoords->push_back(xyz);
      
    MV_Coords(ev[1],coords);
    xyz.set(coords[0],coords[1]);
    fcoords->push_back(xyz);
  }	  
} // Mesh_MSTK::face_get_coordinates
  





// Get list of entities of type 'category' in set specified by setname

void Mesh_MSTK::get_set_entities (const std::string setname, 
				  const Entity_kind kind, 
				  const Parallel_type ptype, 
				  std::vector<Entity_ID> *setents) const 
{
  int idx, i, lid;
  MSet_ptr mset, mset1;
  MEntity_ptr ment;
  bool found(false);
  int celldim = Mesh::cell_dimension();
  int spacedim = Mesh::space_dimension();

  assert(setents != NULL);
  
  setents->clear();

  AmanziGeometry::GeometricModelPtr gm = Mesh::geometric_model();

  // Is there an appropriate region by this name?

  AmanziGeometry::RegionPtr rgn = gm->FindRegion(setname);

  // Did not find the region
  
  if (rgn == NULL) 
    {
      std::cerr << "Geometric model has no region named " << setname << std::endl;
      std::cerr << "Cannot construct set by this name" << std::endl;
      throw std::exception();
    }



  // Region is of type labeled set and a mesh set should have been
  // initialized from the input file
  
  if (rgn->type() == AmanziGeometry::LABELEDSET)
    {
      AmanziGeometry::LabeledSetRegionPtr lsrgn = dynamic_cast<AmanziGeometry::LabeledSetRegionPtr> (rgn);
      std::string label = lsrgn->label();
      std::string entity_type = lsrgn->entity_str();

      if ((kind == CELL && entity_type != "CELL") ||
          (kind == FACE && entity_type != "FACE") ||
          (kind == NODE && entity_type != "NODE"))
        {
          std::cerr << "Found labeled set region named " << setname << " but it" << std::endl;
          std::cerr << "contains entities of type " << entity_type << ", not the requested type" << std::endl;

          throw std::exception();
        } 

      char internal_name[256];

      if (entity_type == "CELL")
        sprintf(internal_name,"matset_%s",label.c_str());
      else if (entity_type == "FACE")
        sprintf(internal_name,"sideset_%s",label.c_str());
      else if (entity_type == "NODE")
        sprintf(internal_name,"nodeset_%s",label.c_str());

      mset1 = MESH_MSetByName(mesh,internal_name);
      
      if (!mset1)
        {
          std::cerr << "Mesh set " << setname << " should have been read in" << std::endl;
          std::cerr << "as set " << label << " from the mesh file. Its absence" << std::endl;
          std::cerr << "indicates an error in the input file or in the mesh file" << std::endl;
          
          throw std::exception();
        }
    }
  else
    {
      mset1 = MESH_MSetByName(mesh,setname.c_str());
    }

  if (mset1 == NULL) 
    {

      // This mesh entity set did not exist in the input mesh file
      // Create it based on the region defintion

      switch (kind) 
        {

          // cellsets

        case CELL:

          if (rgn->dimension() != celldim) 
            {
              std::cerr << "No region of dimension " << celldim << " defined in geometric model" << std::endl;
              std::cerr << "Cannot construct cell set from region " << setname << std::endl;
            }


          if (rgn->type() == AmanziGeometry::BOX) 
            {
              if (celldim == 3)
                mset1 = MSet_New(mesh,setname.c_str(),MREGION);
              else
                mset1 = MSet_New(mesh,setname.c_str(),MFACE);

              int ncell = num_entities(CELL, USED);
              
              for (int icell = 0; icell < ncell; icell++)
                {
                  if (rgn->inside(cell_centroid(icell)))
                    {
                      MSet_Add(mset1,cell_id_to_handle[icell]);
                    }
                }
            }
          else
            {
              std::cerr << "Region type not applicable/supported for cell sets" << std::endl;
              throw std::exception();
            }
      
          break;
      

          // sidesets

        case FACE:

          if (rgn->dimension() != celldim-1) 
            {
              std::cerr << "No region of dimension " << celldim-1 << " defined in geometric model" << std::endl;
              std::cerr << "Cannot construct cell set from region " << setname << std::endl;
            }

          if (rgn->type() == AmanziGeometry::BOX) 
            {
              if (celldim == 3)   // 3D meshes, 2D sidesets
                mset1 = MSet_New(mesh,setname.c_str(),MFACE);
              else
                mset1 = MSet_New(mesh,setname.c_str(),MEDGE);

              int nface = num_entities(FACE, USED);
              
              for (int iface = 0; iface < nface; iface++)
                {
                  if (rgn->inside(face_centroid(iface)))
                    {
                      MSet_Add(mset1,face_id_to_handle[iface]);
                    }
                }
            }
          else if (rgn->type() == AmanziGeometry::PLANE) 
            {
              if (celldim == 3)
                mset1 = MSet_New(mesh,setname.c_str(),MFACE);
              else
                mset1 = MSet_New(mesh,setname.c_str(),MEDGE);

              
              int nface = num_entities(FACE, USED);
              
              for (int iface = 0; iface < nface; iface++)
                {
                  std::vector<AmanziGeometry::Point> fcoords(spacedim);

                  face_get_coordinates(iface, &fcoords);

                  bool on_plane = true;
                  for (int j = 0; j < fcoords.size(); j++)
                    {
                      if (!rgn->inside(fcoords[j]))
                          {
                            on_plane = false;
                            break;
                          }
                    }
                  
                  if (on_plane)
                    {
                      MSet_Add(mset1,face_id_to_handle[iface]);
                    }
                }
            }
          else 
            {
              std::cerr << "Region type not applicable/supported for face sets" << std::endl;
              throw std::exception();
            }

          break;


          // Nodesets

        case NODE:

          if (rgn->type() == AmanziGeometry::BOX ||
              rgn->type() == AmanziGeometry::PLANE) 
            {
              mset1 = MSet_New(mesh,setname.c_str(),MVERTEX);

              int nnode = num_entities(NODE, USED);

              for (int inode = 0; inode < nnode; inode++)
                {
                  AmanziGeometry::Point vpnt(spacedim);

                  node_get_coordinates(inode, &vpnt);
                  
                  if (rgn->inside(vpnt)) 
                    {
                      MSet_Add(mset1,vtx_id_to_handle[inode]);
                    }
                }
            }
          else 
            {
              std::cerr << "Region type not applicable/supported for node sets" << std::endl;
              throw std::exception();
            }

          break;
        }
    }
  

  switch (kind) {
  case CELL: {
    int ncells;
    int *cells;

    if (!MSet_Num_Entries(mset1)) return;

    if (ptype == OWNED)
      mset = MSets_Subtract(mset1,GhostCells);
    else if (ptype == GHOST)
      mset = MSets_Subtract(mset1,OwnedCells);
    else
      mset = mset1;

    if (!mset || !MSet_Num_Entries(mset)) return;

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

    if (!MSet_Num_Entries(mset1)) return;

    if (ptype == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedFaces);
    else if (ptype == GHOST)
      mset = MSets_Subtract(mset1,OwnedFaces);
    else
      mset = mset1;

    if (!mset || !MSet_Num_Entries(mset)) return;

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

    if (!MSet_Num_Entries(mset1)) return;
    
    if (ptype == OWNED)
      mset = MSets_Subtract(mset1,NotOwnedVerts);
    else if (ptype == GHOST)
      mset = MSets_Subtract(mset1,OwnedVerts);
    else
      mset = mset1;

    if (!mset || !MSet_Num_Entries(mset)) return;

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
  
} // Mesh_MSTK::get_set_entities (by set name) 


void Mesh_MSTK::get_set_entities (const char *setname,
				  const Entity_kind kind, 
				  const Parallel_type ptype, 
				  std::vector<Entity_ID> *setents) const 
{
  std::string setname1(setname);

  get_set_entities(setname1,kind,ptype,setents);
  
} // Mesh_MSTK::get_set_entities (by set name) 


void Mesh_MSTK::get_set_entities (const Set_ID setid, 
				  const Entity_kind kind, 
				  const Parallel_type ptype, 
				  std::vector<Entity_ID> *setents) const 
{
  AmanziGeometry::GeometricModelPtr gm = Mesh::geometric_model();
  AmanziGeometry::RegionPtr rgn = gm->FindRegion(setid);

  std::cerr << "DEPRECATED METHOD!" << std::endl;
  std::cerr << "Call get_set_entities with setname instead of setid" << std::endl;

  if (!rgn)
    {
      std::cerr << "No region with id" << setid << std::endl;
    }

  get_set_entities(rgn->name(),kind,ptype,setents);
  
} // Mesh_MSTK::get_set_entities (by set id) 



// Get number of entities of type 'ptype' in set

unsigned int Mesh_MSTK::get_set_size (const char *setname,
                                      const Entity_kind kind, 
                                      const Parallel_type ptype) const 
{
  Entity_ID_List setents;
  std::string setname1(setname);

  get_set_entities(setname1,kind,ptype,&setents);
  
  return setents.size();
  
} // Mesh_MSTK::get_set_size (by set name)

// Get number of entities of type 'ptype' in set

unsigned int Mesh_MSTK::get_set_size (const std::string setname, 
                                      const Entity_kind kind, 
                                      const Parallel_type ptype) const 
{
  Entity_ID_List setents;

  get_set_entities(setname,kind,ptype,&setents);
  
  return setents.size();
  
} // Mesh_MSTK::get_set_size (by set name)



// Get number of entities of type 'ptype' in set

unsigned int Mesh_MSTK::get_set_size (const Set_ID setid, 
                                      const Entity_kind kind, 
                                      const Parallel_type ptype) const 
{
  
  Entity_ID_List setents;

  get_set_entities(setid,kind,ptype,&setents);
  
  return setents.size();
  
} // Mesh_MSTK::get_set_size (by set id)




// Epetra map for cells - basically a structure specifying the
// global IDs of cells owned or used by this processor

// Amanzi/Epetra want global IDs to start at 0

void Mesh_MSTK::init_cell_map ()
{
  int *cell_gids;
  int ncell, idx, i;
  MEntity_ptr ment;
  const Epetra_Comm *epcomm = get_comm();

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

} // Mesh_MSTK::init_cell_map




// Epetra map for faces - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_MSTK::init_face_map ()
{
  int *face_gids;
  int nface, idx, i;
  MEntity_ptr ment;
  const Epetra_Comm *epcomm = get_comm();

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

} // Mesh_MSTK::init_face_map




// Epetra map for nodes - basically a structure specifying the
// global IDs of cells owned or used by this processor

void Mesh_MSTK::init_node_map ()
{
  int *vert_gids;
  int nvert, idx, i;
  MEntity_ptr ment;
  const Epetra_Comm *epcomm = get_comm();

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

} // Mesh_MSTK::init_node_map


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

} // Mesh_MSTK::GID





// Some initializations

void Mesh_MSTK::clear_internals_ () 
{ 

  faceflip = NULL;

  cell_map_w_ghosts_ = cell_map_wo_ghosts_ = NULL;
  face_map_w_ghosts_ = face_map_wo_ghosts_ = NULL;
  node_map_w_ghosts_ = node_map_wo_ghosts_ = NULL;

} // Mesh_MSTK::clear_internals




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
    
} // Mesh_MSTK::init_id_handle_maps



void Mesh_MSTK::init_global_ids() {
    
} // Mesh_MSTK::init_global_ids


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

} // Mesh_MSTK::init_pvert_lists


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
} // Mesh_MSTK::init_pface_lists


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



    MSTK_UpdateAttr(mesh, myprocid, numprocs, mpicomm);



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
} // Mesh_MSTK::init_pface_dirs


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
} // Mesh_MSTK::init_pcell_lists





void Mesh_MSTK::init_set_info() {
  int idx;
  MSet_ptr mset;
  char setname[256];
  
  AmanziGeometry::GeometricModelPtr gm = Mesh::geometric_model();

  if (gm == NULL) { 
    std::cerr << "Need region definitions to initialize sets" << std::endl;
    return;
  }
    

  unsigned int ngr = gm->Num_Regions();


  for (int i = 0; i < ngr; i++) {

    AmanziGeometry::RegionPtr rgn = gm->Region_i(i);

    strcpy(setname,rgn->name().c_str());

    mset = MESH_MSetByName(mesh,setname);

    if (mset) 
      {

	// The only time a mesh set by this name should already exist is
	// if the region is of type labeled set. In that case, verify
	// that the entity types in the set are the same as the one
	// requested in the region

	if (rgn->type() == AmanziGeometry::LABELEDSET) 
	  {
            AmanziGeometry::LabeledSetRegionPtr lsrgn =
              dynamic_cast<AmanziGeometry::LabeledSetRegionPtr> (rgn);

	    MType entdim = MSet_EntDim(mset);
	    if (Mesh::space_dimension() == 3) 
	      {
		if ((lsrgn->entity_str() == "CELL" && entdim != MREGION) ||
		    (lsrgn->entity_str() == "FACE" && entdim != MFACE) ||
		    (lsrgn->entity_str() == "NODE" && entdim != MVERTEX))
		  {
		    std::cerr << "Mismatch of entity type in labeled set region and mesh set" << std::endl;
		    throw std::exception();
		  }
	      }
	    else if (Mesh::space_dimension() == 2)
	      {
		if ((lsrgn->entity_str() == "CELL" && entdim != MFACE) ||
		    (lsrgn->entity_str() == "FACE" && entdim != MEDGE) ||
		    (lsrgn->entity_str() == "NODE" && entdim != MVERTEX))
		  {
		    std::cerr << "Mismatch of entity type in labeled set region and mesh set" << std::endl;
		    throw std::exception();
		  }
	      }
	  }
      }
    else 
      { 
	// Create it on demand later
      }
  }
  

}  /* Mesh_MSTK::init_set_info */



void Mesh_MSTK::collapse_degen_edges() {
  const int topoflag=1;
  std::vector<unsigned int> *orig_vids;
  int idx, idx2, evgid0, evgid1;
  MVertex_ptr ev0, ev1, vkeep, vdel;
  MEdge_ptr edge;
  MFace_ptr face;
  MRegion_ptr region;
  List_ptr eregs, efaces, vregs, vfaces;
  double len2;
  int ival;
  void *pval;
  Cell_type celltype;

  if (cell_dimension() == 2) {
    orig_celltype_att = MAttrib_New(mesh,"Orig_Cell_type",INT,MFACE);
    orig_celltopo_att = MAttrib_New(mesh,"Orig_Cell_Topo",POINTER,MFACE);
  }
  else {
    orig_celltype_att = MAttrib_New(mesh,"Orig_Cell_type",INT,MREGION);
    orig_celltopo_att = MAttrib_New(mesh,"Orig_Cell_Topo",POINTER,MREGION);
  }

  idx = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx))) {

    len2 = ME_LenSqr(edge);

    if (len2 <= 1.0e-32) {

      /* Degenerate edge  - must collapse */

      /* First record the cell configuration and store it for use by viz */

      eregs = ME_Regions(edge);
      if (eregs) {
	       
	idx2 = 0;
	while ((region = List_Next_Entry(eregs,&idx2))) {

	  orig_vids = new std::vector<unsigned int>;

	  cell_get_nodes(MR_ID(region),orig_vids);

	  MEnt_Get_AttVal(region,celltype_att,&ival,NULL,NULL);
	  celltype = (Cell_type) ival;

	  MEnt_Set_AttVal(region,orig_celltype_att,celltype,0.0,NULL);
	  MEnt_Set_AttVal(region,orig_celltopo_att,0,0.0,orig_vids);
	}

	List_Delete(eregs);
      }
      else {

	efaces = ME_Faces(edge);

	idx2 = 0;
	while ((face = List_Next_Entry(efaces,&idx2))) {

	  orig_vids = new std::vector<unsigned int>;

	  cell_get_nodes(MF_ID(face),orig_vids);

	  MEnt_Get_AttVal(face,celltype_att,&ival,NULL,NULL);
	  celltype = (Cell_type) ival;

	  MEnt_Set_AttVal(face,orig_celltype_att,celltype,0.0,NULL);
	  MEnt_Set_AttVal(face,orig_celltopo_att,0,0.0,orig_vids);
	}

	List_Delete(efaces);
      }
      

      /* Now collapse choosing the vertex to be deleted and vertex to
	 be kept consistently. If topological constraints permit,
	 collapse the vertex with the higher global ID to the vertex
	 with the lower global ID. If they do not, reverse the
	 order. Since global IDs and topological constraints are the
	 same for master and slave edges and their nodes, we will not
	 have conflict between processors */

      ev0 = ME_Vertex(edge,0); evgid0 = MEnt_GlobalID(ev0);
      ev1 = ME_Vertex(edge,1); evgid1 = MEnt_GlobalID(ev1);

      if (evgid0 < evgid1) {
	vkeep = ev0;
	vdel = ev1;
      }
      else {
	vkeep = ev1;
	vdel = ev0;
      }

      vkeep = ME_Collapse(edge, vkeep, topoflag);

      if (!vkeep) {
	vkeep = vdel;
	vdel = (vkeep == ev0) ? ev1 : ev1;

	vkeep = ME_Collapse(edge, vkeep, topoflag);
      }

      if (!vkeep) {
	cout << "Could not collapse degenerate edge. Expect computational issues with connected elements" << std::endl;
	return;
      }


      /* In the stored configuration of the original cell, replace the
	 ID of the deleted vertex with the retained vertex */

      vregs = MV_Regions(vkeep);
      if (vregs) {

	idx = 0;
	while ((region = List_Next_Entry(vregs,&idx))) {
	  MEnt_Get_AttVal(region,orig_celltopo_att,NULL,NULL,&pval);
	  orig_vids = (std::vector<Entity_ID> *)pval;

          // Will this work? Don't we have to use push_back?

	  for (int i = 0; i < orig_vids->size(); i++) 
	    if ((*orig_vids)[i] == MV_ID(vdel)) {
	        throw std::exception();
	      (*orig_vids)[i] = (Entity_ID) MV_ID(vkeep);
            }
	}

	List_Delete(vregs);
      }
      else {
	
	vfaces = MV_Faces(vkeep);
	if (vfaces) {

	  idx = 0;
	  while ((face = List_Next_Entry(vfaces,&idx))) {
	    MEnt_Get_AttVal(region,orig_celltopo_att,NULL,NULL,&pval);
	    orig_vids = (std::vector<Entity_ID> *)pval;

            // will this work? Don't we have to use push_back
	    for (int i = 0; i < orig_vids->size(); i++)
	      if ((*orig_vids)[i] = MV_ID(vdel)) {
	        throw std::exception();
		(*orig_vids)[i] = (Entity_ID) MV_ID(vkeep);
              }
	  }
	}
      }

    }

  }

} /* end Mesh_MSTK::collapse_degen_edges */


void Mesh_MSTK::label_celltype() {
  Cell_type ctype;
  int idx, idx2, nrv, nrf, nfv, nquads;
  MFace_ptr face;
  MRegion_ptr region;
  List_ptr rfaces, rverts;

  if (cell_dimension() == 2) 
    celltype_att = MAttrib_New(mesh,"Cell_type",INT,MFACE);
  else
    celltype_att = MAttrib_New(mesh,"Cell_type",INT,MREGION);

  if (cell_dimension() == 2) {

    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      nfv = MF_Num_Vertices(face);

      switch (nfv) {
      case 3:
	ctype = TRI;
	break;
      case 4:
	ctype = QUAD;
	break;
      default:
	ctype = POLYGON;
      }
      MEnt_Set_AttVal(face,celltype_att,ctype,0.0,NULL);
    }
      
  }
  else if (cell_dimension() == 3) {

    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {

      rverts = MR_Vertices(region);
      nrv = List_Num_Entries(rverts);
      List_Delete(rverts);

      nrf = MR_Num_Faces(region);

      switch (nrf) {
      case 4:
	ctype = TET;
	break;
      case 5:

	  nquads = 0;
	  rfaces = MR_Faces(region);	  
	  idx2 = 0;
	  while ((face = List_Next_Entry(rfaces,&idx2)))
	    if (MF_Num_Vertices(face) == 4)
	      nquads++;
	  List_Delete(rfaces);

	  switch (nquads) {
	  case 1:
	    ctype = PYRAMID;
	    break;
	  case 3:
	    ctype = PRISM;
	    break;
	  default:
	    ctype = POLYHED;
	  }

	break;

      case 6:

	  nquads = 0;
	  rfaces = MR_Faces(region);	  
	  idx2 = 0;
	  while ((face = List_Next_Entry(rfaces,&idx2)))
	    if (MF_Num_Vertices(face) == 4)
	      nquads++;
	  List_Delete(rfaces);

	  ctype = (nquads == 6) ? HEX : POLYHED;

	break;
      default:

	ctype = POLYHED;

      }

      MEnt_Set_AttVal(region,celltype_att,ctype,0.0,NULL);

    }

  }

} /* Mesh_MSTK::label_celltypes */


//
// Epetra maps
//------------
    
    
inline 
const Epetra_Map& Mesh_MSTK::cell_epetra_map (const bool include_ghost) const
{
  if (serial_run)
    return *cell_map_wo_ghosts_;
  else
    return (include_ghost ? *cell_map_w_ghosts_ : *cell_map_wo_ghosts_);
}
    
inline 
const Epetra_Map& Mesh_MSTK::face_epetra_map (const bool include_ghost) const
{
  if (serial_run)
    return *face_map_wo_ghosts_;
  else
    return (include_ghost ? *face_map_w_ghosts_ : *face_map_wo_ghosts_);
}
    
inline 
const Epetra_Map& Mesh_MSTK::node_epetra_map (const bool include_ghost) const
{
  if (serial_run)
    return *node_map_wo_ghosts_;
  else
    return (include_ghost ? *node_map_w_ghosts_ : *node_map_wo_ghosts_);
}


} // close namespace AmanziMesh
} // close namespace Amanzi
