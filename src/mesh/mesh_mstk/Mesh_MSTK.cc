//#include <Teuchos_RCP.hpp>

#include "dbc.hh"

#include "Mesh_MSTK.hh"


using namespace std;


namespace Amanzi
{

namespace AmanziMesh
{


//--------------------------------------
// Constructor - load up mesh from file
//--------------------------------------

Mesh_MSTK::Mesh_MSTK (const char *filename, const Epetra_MpiComm *incomm,
		      const AmanziGeometry::GeometricModelPtr& gm) :
  mpicomm(incomm->GetMpiComm())
{  

  // Assume three dimensional problem if constructor called without 
  // the space_dimension parameter

  int ok;


  // Pre-processing (init, MPI queries etc)

  int space_dim = 3;
  pre_create_steps_(space_dim, incomm, gm);




  if (myprocid == 0) {
    int DebugWait=0;
    while (DebugWait);
  }


  if (serial_run) {

    // Load serial mesh

    mesh = MESH_New(F1);
    ok = MESH_ImportFromExodusII(mesh,filename);

    set_cell_dimension((MESH_Num_Regions(mesh) ? 3 : 2));

    if (cell_dimension() == 2 && space_dim == 3) {
      
      // Check if this is a completely planar mesh 
      // in which case one can label the space dimension as 2
      
      MVertex_ptr mv, mv0 = MESH_Vertex(mesh,0);
      double vxyz[3], z0;

      MV_Coords(mv0,vxyz);
      z0 = vxyz[2];

      bool planar = true;
      int idx = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx))) {
        MV_Coords(mv,vxyz);
        if (z0 != vxyz[2]) {
          planar = false;
          break;
        }
      }
       
      if (planar) {
        space_dim = 2;
        set_space_dimension(space_dim);
      }
    }

    myprocid = 0;
  }
  else {
    Mesh_ptr globalmesh;
    int topo_dim=3; // What is the topological dimension of the mesh
    int ring = 1; // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes
    
    if (myprocid == 0) {
      globalmesh = MESH_New(F1);
      ok = MESH_ImportFromExodusII(globalmesh,filename);
      
      topo_dim = MESH_Num_Regions(globalmesh) ? 3 : 2;
      
      mesh = globalmesh;

      if (cell_dimension() == 2 && space_dim == 3) {
        
        // Check if this is a completely planar mesh 
        // in which case one can label the space dimension as 2
        
        MVertex_ptr mv, mv0 = MESH_Vertex(mesh,0);
        double vxyz[3], z0;
        
        MV_Coords(mv0,vxyz);
        z0 = vxyz[2];
        
        bool planar = true;
        int idx = 0;
        while ((mv = MESH_Next_Vertex(mesh,&idx))) {
          MV_Coords(mv,vxyz);
          if (z0 != vxyz[2]) {
            planar = false;
            break;
          }
        }
        
        if (planar) {
          space_dim = 2;
          MPI_Scatter(&space_dim,1,MPI_INT,&space_dim,1,MPI_INT,0,mpicomm);
        }
      }

    }
    else {
      mesh = MESH_New(UNKNOWN_REP);
      ok = 1;
    }

    ok = ok & MSTK_Mesh_Distribute(&mesh,&topo_dim,ring,with_attr,myprocid,
                                   numprocs,mpicomm);

    if (myprocid == 0)
      MESH_Delete(globalmesh);

    set_cell_dimension(topo_dim);
    set_space_dimension(space_dim);
  }

  if (!ok) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << "Failed to load " << filename << " on processor " << myprocid << std::endl;
    assert(ok);
  }



  // Do all the processing required for setting up the mesh for Amanzi 
  
  post_create_steps_();

}


//--------------------------------------
// Constructor - load up mesh from file
//--------------------------------------

Mesh_MSTK::Mesh_MSTK (const char *filename, const Epetra_MpiComm *incomm, 
		      int space_dimension,
		      const AmanziGeometry::GeometricModelPtr& gm) :
  mpicomm(incomm->GetMpiComm())
{
  int ok;

  pre_create_steps_(space_dimension, incomm, gm);


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
    int topo_dim=3; // What is the topological dimension of the mesh
    int ring = 1; // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes
    
    if (myprocid == 0) {
      globalmesh = MESH_New(F1);
      ok = MESH_ImportFromExodusII(globalmesh,filename);
      
      topo_dim = MESH_Num_Regions(globalmesh) ? 3 : 2;
      
      mesh = globalmesh;
    }
    else {
      mesh = MESH_New(UNKNOWN_REP);
      ok = 1;
    }

    ok = ok & MSTK_Mesh_Distribute(&mesh,&topo_dim,ring,with_attr,myprocid,
                                   numprocs,mpicomm);

    if (myprocid == 0)
      MESH_Delete(globalmesh);
    
    set_cell_dimension(topo_dim);
  }

  if (!ok) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << "Failed to load " << filename << " on processor " << myprocid << std::endl;
    assert(ok);
  }





  // Do all the processing required for setting up the mesh for Amanzi 
  
  post_create_steps_();

}


//--------------------------------------
// Construct a 3D regular hexahedral mesh internally
//--------------------------------------


Mesh_MSTK::Mesh_MSTK(const double x0, const double y0, const double z0,
		     const double x1, const double y1, const double z1,
		     const unsigned int nx, const unsigned int ny, 
		     const unsigned int nz, 
		     const Epetra_MpiComm *incomm,
		     const AmanziGeometry::GeometricModelPtr& gm) :
  mpicomm(incomm->GetMpiComm())
{
  int ok;

  int space_dimension = 3;
  pre_create_steps_(space_dimension, incomm, gm);




  if (serial_run) {

    // Load serial mesh

    mesh = MESH_New(F1);
    ok = generate_regular_mesh(mesh,x0,y0,z0,x1,y1,z1,nx,ny,nz);

    set_cell_dimension(3);

    myprocid = 0;
  }
  else {
    Mesh_ptr globalmesh;
    int topo_dim=3; // What is the topological dimension of the mesh
    int ring = 1; // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes

    
    if (myprocid == 0) {
      globalmesh = MESH_New(F1);

      ok = generate_regular_mesh(globalmesh,x0,y0,z0,x1,y1,z1,nx,ny,nz);
      
      mesh = globalmesh;
    }
    else {
      mesh = MESH_New(UNKNOWN_REP);
      ok = 1;
    }

    ok = ok & MSTK_Mesh_Distribute(&mesh,&topo_dim,ring,with_attr,myprocid,
					numprocs,mpicomm);

    if (myprocid == 0)
      MESH_Delete(globalmesh);

    set_cell_dimension(topo_dim);
  }

  if (!ok) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << "Failed to generate mesh on processor " << myprocid << std::endl;
    assert(ok);
  }




  // Do all the processing required for setting up the mesh for Amanzi 
  
  post_create_steps_();

}



//--------------------------------------
// Construct a 2D regular quadrilateral mesh internally
//--------------------------------------


Mesh_MSTK::Mesh_MSTK(const double x0, const double y0,
		     const double x1, const double y1,
		     const int nx, const int ny, 
		     const Epetra_MpiComm *incomm,
		     const AmanziGeometry::GeometricModelPtr& gm) :
  mpicomm(incomm->GetMpiComm())
{
  int ok;

  int space_dim = 2;
  pre_create_steps_(space_dim, incomm, gm);




  if (myprocid == 0) {
    int DebugWait=0;
    while (DebugWait);
    }


  int topo_dim=space_dim; // What is the topological dimension of the mesh
  set_cell_dimension(topo_dim);

  if (serial_run) {

    // Load serial mesh

    mesh = MESH_New(F1);
    ok = generate_regular_mesh(mesh,x0,y0,x1,y1,nx,ny);


    myprocid = 0;
  }
  else {
    Mesh_ptr globalmesh;
    int ring = 1; // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes

    if (myprocid == 0) {
      globalmesh = MESH_New(F1);

      ok = generate_regular_mesh(globalmesh,x0,y0,x1,y1,nx,ny);
      
      mesh = globalmesh;
    }
    else {
      mesh = MESH_New(UNKNOWN_REP);
      ok = 1;
    }

    ok = ok & MSTK_Mesh_Distribute(&mesh,&topo_dim,ring,with_attr,myprocid,
					numprocs,mpicomm);

    if (myprocid == 0)
      MESH_Delete(globalmesh);
  }

  if (!ok) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << "Failed to generate mesh on processor " << myprocid << std::endl;
    assert(ok);
  }




  // Do all the processing required for setting up the mesh for Amanzi 
  
  post_create_steps_();

}



//-------------------------------------- 
// Construct a 2D or 3D regular mesh using input from the
// GenerationSpec class 
//--------------------------------------

Mesh_MSTK::Mesh_MSTK(const GenerationSpec& gspec,
		     const Epetra_MpiComm *incomm,
		     const AmanziGeometry::GeometricModelPtr& gm) :
  mpicomm(incomm->GetMpiComm())
{
  int ok;


  // Get info about the domain from the generation specification class

  AmanziGeometry::Point p0(gspec.domain().point0());
  AmanziGeometry::Point p1(gspec.domain().point1());

  int space_dim = p0.dim();
  pre_create_steps_(space_dim, incomm, gm);




  if (myprocid == 0) {
    int DebugWait=0;
    while (DebugWait);
  }

  int topo_dim=space_dim;
  set_cell_dimension(topo_dim);

  if (serial_run) {

    // Load serial mesh

    mesh = MESH_New(F1);

    if (topo_dim == 2) {
      ok = generate_regular_mesh(mesh,p0.x(),p0.y(),p1.x(),p1.y(),
				 gspec.xcells(),gspec.ycells());
    }
    else if (topo_dim == 3) {
      ok = generate_regular_mesh(mesh,p0.x(),p0.y(),p0.z(),p1.x(),p1.y(),p1.z(),
				 gspec.xcells(),gspec.ycells(),gspec.zcells());
    }


    myprocid = 0;
  }
  else {
    Mesh_ptr globalmesh;
    int ring = 1; // One layer of ghost cells in parallel meshes
    int with_attr = 1;  // update of attributes in parallel meshes
    
    if (myprocid == 0) {
      globalmesh = MESH_New(F1);

      if (topo_dim == 2) {
	ok = generate_regular_mesh(globalmesh,p0.x(),p0.y(),p1.x(),p1.y(),
				   gspec.xcells(),gspec.ycells());
      }
      else if (topo_dim == 3) {
	ok = generate_regular_mesh(globalmesh,
				   p0.x(),p0.y(),p0.z(),p1.x(),p1.y(),p1.z(),
				   gspec.xcells(),gspec.ycells(),gspec.zcells());
      }
      
      mesh = globalmesh;
    }
    else {
      mesh = MESH_New(UNKNOWN_REP);
      ok = 1;
    }

    ok = ok & MSTK_Mesh_Distribute(&mesh,&topo_dim,ring,with_attr,myprocid,
					numprocs,mpicomm);

    if (myprocid == 0)
      MESH_Delete(globalmesh);
  }

  if (!ok) {
    std::cerr << "FAILED" << std::endl;
    std::cerr << "Failed to generate mesh on processor " << myprocid << std::endl;
    assert(ok);
  }




  // Do all the processing required for setting up the mesh for Amanzi 
  
  post_create_steps_();

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

  MESH_Delete(mesh);
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
    

// For the sake of efficiency, this routine is using implicit
// knowledge of how MSTK operates internally and computes MR_Vertices

// Also, for the sake of efficiency, this code is repeated in 
// cell_get_face_dirs 
// If you change this routine (cell_get_faces), check if you 
// need to make the same changes in (cell_get_face_dirs)


void Mesh_MSTK::cell_get_faces (const Entity_ID cellid, 
				std::vector<Entity_ID> *faceids) const
{
  MEntity_ptr cell;

  assert(faceids != NULL);

  faceids->clear();

  cell = cell_id_to_handle[cellid];
      
  if (cell_dimension() == 3) {

    int celltype = cell_get_type(cellid);

    List_ptr rfaces = MR_Faces((MRegion_ptr)cell);   

    /* base face */

    MFace_ptr face0 = List_Entry(rfaces,0);
    int fdir0 = MR_FaceDir_i((MRegion_ptr)cell,0);


    if (celltype >= TET && celltype <= HEX) {
      int lid;
      
      int mkid = MSTK_GetMarker();
      MEnt_Mark(face0,mkid);

      /* Add all lateral faces first (faces adjacent to the base face) */

      List_ptr fedges0 = MF_Edges(face0,!fdir0,0);
      int idx = 0;
      MEdge_ptr fe;
      while ((fe = List_Next_Entry(fedges0,&idx))) {

        /* Is there an unprocessed face in this region that is
           adjacent to this edge */
        
        int idx2 = 0;
        MFace_ptr fadj;
        while ((fadj = List_Next_Entry(rfaces,&idx2))) {

          if (!MEnt_IsMarked(fadj,mkid)) {

            if (MF_UsesEntity(fadj,fe,MEDGE)) {
              lid = MEnt_ID(fadj);
              faceids->push_back(lid-1);
              
              MEnt_Mark(fadj,mkid);
            }
          }

        }
      }
      List_Delete(fedges0);

      /* Add the base face */

      lid = MEnt_ID(face0);
      faceids->push_back(lid-1);

      /* If there is a last remaining face, it is the top face */

      MFace_ptr fopp;
      idx = 0;
      while ((fopp = List_Next_Entry(rfaces,&idx))) {

        if (!MEnt_IsMarked(fopp,mkid)) {
          lid = MEnt_ID(fopp);
          faceids->push_back(lid-1);
        }

      }


      List_Unmark(rfaces,mkid);
      MSTK_FreeMarker(mkid);
    }
    else {
      int idx = 0;
      MFace_ptr face;
      while ((face = List_Next_Entry(rfaces,&idx))) {
	int lid = MEnt_ID(face);
	faceids->push_back(lid-1);
      }
    }
    
    List_Delete(rfaces);
  }
  else {  // cell_dimension() = 2; surface or 2D mesh

    List_ptr fedges = MF_Edges((MFace_ptr)cell,1,0);
    MEdge_ptr edge;
    int idx = 0;
    while ((edge = List_Next_Entry(fedges,&idx))) {
      int lid = MEnt_ID(edge);
      faceids->push_back(lid-1);
    }
    List_Delete(fedges);
  }
  
} // Mesh_MSTK::cell_get_faces



// Get directions in which a cell uses face
// In 3D, direction is 1 if face normal points out of cell
// and -1 if face normal points into cell
// In 2D, direction is 1 if face/edge is defined in the same
// direction as the cell polygon, and -1 otherwise

// For the sake of efficiency, we are repeating the code from above
// If you change this routine (cell_get_face_dirs), check if you 
// need to make the same changes in (cell_get_faces)
    
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

    int celltype = cell_get_type(cellid);

    List_ptr rfaces = MR_Faces((MRegion_ptr)cell);   

    /* base face */

    MFace_ptr face0 = List_Entry(rfaces,0);
    int fdir0 = MR_FaceDir_i((MRegion_ptr)cell,0);


    if (celltype >= TET && celltype <= HEX) {
      int lid;

      int mkid = MSTK_GetMarker();
      MEnt_Mark(face0,mkid);

      /* Add all lateral faces first (faces adjacent to the base face) */

      List_ptr fedges0 = MF_Edges(face0,!fdir0,0);
      int idx = 0;
      MEdge_ptr fe;
      while ((fe = List_Next_Entry(fedges0,&idx))) {

        /* Is there an unprocessed face in this region that is
           adjacent to this edge */
        
        int idx2 = 0;
        MFace_ptr fadj; 
        int i = 0;
        while ((fadj = List_Next_Entry(rfaces,&idx2))) {

          if (!MEnt_IsMarked(fadj,mkid)) {

            if (MF_UsesEntity(fadj,fe,MEDGE)) {

              int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;
              
              lid = MEnt_ID(fadj);
              if (faceflip[lid-1]) fdir *= -1;
              
              face_dirs->push_back(fdir);
              
              MEnt_Mark(fadj,mkid);
            }
          }

          i++;
        }
      }
      List_Delete(fedges0);

      /* Add the base face */

      lid = MEnt_ID(face0);
      fdir0 = fdir0 ? 1 : -1;
      if (faceflip[lid-1]) fdir0 *= -1;
      face_dirs->push_back(fdir0);

      /* If there is a last remaining face, it is the top face */

      MFace_ptr fopp;
      idx = 0;
      int i = 0;
      while ((fopp = List_Next_Entry(rfaces,&idx))) {

        if (!MEnt_IsMarked(fopp,mkid)) {

          int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;

          lid = MEnt_ID(fopp);
          if (faceflip[lid-1]) fdir *= -1;

          face_dirs->push_back(fdir);
        }

        i++;
      }


      List_Unmark(rfaces,mkid);
      MSTK_FreeMarker(mkid);
    }
    else {
      int idx = 0;
      MFace_ptr face;
      int i = 0;
      while ((face = List_Next_Entry(rfaces,&idx))) {

        int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;
        
	int lid = MEnt_ID(face);
        if (faceflip[lid-1]) fdir *= -1;

	face_dirs->push_back(fdir);

        i++;
      }
    }
    
    List_Delete(rfaces);
  }
  else {  // cell_dimension() = 2; surface or 2D mesh

    List_ptr fedges = MF_Edges((MFace_ptr)cell,1,0);

    MEdge_ptr edge;
    int idx = 0;
    int i = 0;
    while ((edge = List_Next_Entry(fedges,&idx))) {

      int fdir = (MF_EdgeDir_i((MFace_ptr)cell,i) == 1) ? 1 : -1;

      int lid = MEnt_ID(edge);
      if (faceflip[lid-1]) fdir *= -1;

      face_dirs->push_back(fdir);

      i++;
    }
    List_Delete(fedges);
  }

} // Mesh_MSTK::cell_get_face_dirs



  // Get faces of a cell and directions in which the cell uses the face 

  // The Amanzi coding guidelines regarding function arguments is purposely
  // violated here to allow for a default input argument

  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order

  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise

void Mesh_MSTK::cell_get_faces_and_dirs (const Entity_ID cellid,
                                         Entity_ID_List *faceids,
                                         std::vector<int> *face_dirs,
					 const bool ordered) const 
{

  MEntity_ptr cell;
  int j,nf, result;

  ASSERT(faceids != NULL);
  ASSERT(face_dirs != NULL);

  faceids->clear();
  face_dirs->clear();
  
  
  cell = cell_id_to_handle[cellid];


  if (cell_dimension() == 3) {

    List_ptr rfaces = MR_Faces((MRegion_ptr)cell);   

    int celltype = cell_get_type(cellid);

    if (ordered && (celltype >= TET && celltype <= HEX)) {
      int lid;


      /* base face */

      MFace_ptr face0 = List_Entry(rfaces,0);
      int fdir0 = MR_FaceDir_i((MRegion_ptr)cell,0);

      /* Markers for faces to avoid searching */

      int mkid = MSTK_GetMarker();
      MEnt_Mark(face0,mkid);

      /* Add all lateral faces first (faces adjacent to the base face) */

      List_ptr fedges0 = MF_Edges(face0,!fdir0,0);
      int idx = 0;
      MEdge_ptr fe;
      while ((fe = List_Next_Entry(fedges0,&idx))) {

        /* Is there an unprocessed face in this region that is
           adjacent to this edge */
        
        int idx2 = 0;
        MFace_ptr fadj; 
        int i = 0;
        while ((fadj = List_Next_Entry(rfaces,&idx2))) {

          if (fadj != face0 && !MEnt_IsMarked(fadj,mkid)) {

            if (MF_UsesEntity(fadj,fe,MEDGE)) {

              int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;
              
              lid = MEnt_ID(fadj);
              if (faceflip[lid-1]) fdir *= -1;
              
              faceids->push_back(lid-1);
              face_dirs->push_back(fdir);
              
              MEnt_Mark(fadj,mkid);
            }
          }

          i++;
        }
      }
      List_Delete(fedges0);

      /* Add the base face */

      lid = MEnt_ID(face0);
      fdir0 = fdir0 ? 1 : -1;
      if (faceflip[lid-1]) fdir0 *= -1;

      faceids->push_back(lid-1);
      face_dirs->push_back(fdir0);

      /* If there is a last remaining face, it is the top face */

      MFace_ptr fopp;
      idx = 0;
      int i = 0;
      while ((fopp = List_Next_Entry(rfaces,&idx))) {

        if (fopp != face0 && !MEnt_IsMarked(fopp,mkid)) {

          int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;

          lid = MEnt_ID(fopp);
          if (faceflip[lid-1]) fdir *= -1;

          faceids->push_back(lid-1);
          face_dirs->push_back(fdir);

          break;
        }

        i++;
      }


      List_Unmark(rfaces,mkid);
      MSTK_FreeMarker(mkid);
    }
    else {
      int idx = 0;
      MFace_ptr face;
      int i = 0;
      while ((face = List_Next_Entry(rfaces,&idx))) {

        int fdir = (MR_FaceDir_i((MRegion_ptr)cell,i) == 1) ? 1 : -1;
        
	int lid = MEnt_ID(face);
        if (faceflip[lid-1]) fdir *= -1;

        faceids->push_back(lid-1);
	face_dirs->push_back(fdir);

        i++;
      }
    }
    
    List_Delete(rfaces);
  }
  else {  // cell_dimension() = 2; surface or 2D mesh

    List_ptr fedges = MF_Edges((MFace_ptr)cell,1,0);

    MEdge_ptr edge;
    int idx = 0;
    int i = 0;
    while ((edge = List_Next_Entry(fedges,&idx))) {

      int fdir = (MF_EdgeDir_i((MFace_ptr)cell,i) == 1) ? 1 : -1;

      int lid = MEnt_ID(edge);
      if (faceflip[lid-1]) fdir *= -1;

      faceids->push_back(lid-1);
      face_dirs->push_back(fdir);

      i++;
    }
    List_Delete(fedges);
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
      nodeids->push_back(MEnt_ID(ME_Vertex(genface,1))-1);
      nodeids->push_back(MEnt_ID(ME_Vertex(genface,0))-1);
    }
    else {
      nodeids->push_back(MEnt_ID(ME_Vertex(genface,0))-1);
      nodeids->push_back(MEnt_ID(ME_Vertex(genface,1))-1);
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
      while ((mf2 = List_Next_Entry(efaces,&idx2))) {
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

    List_Delete(fverts);
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

      List_Delete(fverts);
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
  


// Modify a node's coordinates

void Mesh_MSTK::node_set_coordinates(const AmanziMesh::Entity_ID nodeid, 
                                      const double *coords) {
  MVertex_ptr v = vtx_id_to_handle[nodeid];
  MV_Set_Coords(v,(double *)coords);
}

void Mesh_MSTK::node_set_coordinates(const AmanziMesh::Entity_ID nodeid, 
                                     const AmanziGeometry::Point coords) {
  MVertex_ptr v = vtx_id_to_handle[nodeid];

  double coordarray[3] = {0.0,0.0,0.0};
  for (int i = 0; i < Mesh::space_dimension(); i++)
    coordarray[i] = coords[i];

  MV_Set_Coords(v,(double *)coordarray);
}




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
  const Epetra_Comm *epcomm = get_comm();

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

      // Make sure we retrieved a mesh set with the right kind of entities

      MType entdim;

      switch (kind)
        {
        case CELL:
          if (celldim == 3)
            entdim = MREGION;
          else if (celldim == 2)
            entdim = MFACE;
          break;
        case FACE:
          if (celldim == 3)
            entdim = MFACE;
          else if (celldim == 2)
            entdim = MEDGE;
          break;
        case NODE:
          entdim = MVERTEX;
        }

      // If not, can we find a mesh set with the right name and right
      // kind of entities

      char setname1[256];

      if (mset1 && MSet_EntDim(mset1) != entdim) 
        {
          idx = 0;
          while ((mset1 = MESH_Next_MSet(mesh,&idx))) 
            {
              MSet_Name(mset1,setname1);
              
              if (MSet_EntDim(mset) == entdim &&
                  strcmp(setname1,setname.c_str()) == 0)
                break;
            }
        }
    }

  if (mset1 == NULL) 
    {

      // This mesh entity set did not exist in the input mesh file
      // Create it based on the region defintion

      switch (kind) 
        {

          // cellsets

        case CELL:

          if (rgn->type() == AmanziGeometry::BOX ||
              rgn->type() == AmanziGeometry::COLORFUNCTION) 
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
          else if (rgn->type() == AmanziGeometry::POINT)
            {
              if (celldim == 3)
                mset1 = MSet_New(mesh,setname.c_str(),MREGION);
              else
                mset1 = MSet_New(mesh,setname.c_str(),MFACE);

              AmanziGeometry::Point vpnt(spacedim);
              AmanziGeometry::Point rgnpnt(spacedim);

              rgnpnt = ((AmanziGeometry::PointRegionPtr)rgn)->point();


              int nnode = num_entities(NODE, USED);
              double mindist2 = 1.e+16;
              int minnode = -1;

              int inode;
              for (inode = 0; inode < nnode; inode++)
                {
                  
                  node_get_coordinates(inode, &vpnt);
                  
                  double dist2 = (vpnt-rgnpnt)*(vpnt-rgnpnt);
 
                  if (dist2 < mindist2) 
                    {
                      mindist2 = dist2;
                      minnode = inode;
                      if (mindist2 <= 1.0e-32)
                        break;
                    }
                }

              Entity_ID_List cells, cells1;

              node_get_cells(minnode,USED,&cells);

              int ncells = cells.size();
              for (int ic = 0; ic < ncells; ic++) 
                {
                  Entity_ID icell = cells[ic];
                  
                  // Check if point is contained in cell
                  
                  if (point_in_cell(rgnpnt,icell))
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
              rgn->type() == AmanziGeometry::PLANE ||
              rgn->type() == AmanziGeometry::POINT) 
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

                      // Only one node per point region

                      if (rgn->type() == AmanziGeometry::POINT)
                        break;      
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


  /* Check if no processor got any mesh entities */

  int nent_loc, nent_glob;

  nent_loc = MSet_Num_Entries(mset1);
  epcomm->SumAll(&nent_loc,&nent_glob,1);

  if (nent_glob == 0) {
    std::stringstream stream;
    stream << "Could not retrieve any mesh entities for set " << setname << std::endl;
    Errors::Message mesg(stream.str());
    Exceptions::amanzi_throw(mesg);
  }
  

  mset = NULL;
  switch (kind) {
  case CELL: {
    int ncells;
    int *cells;

    if (nent_loc) {
      if (ptype == OWNED)
        mset = MSets_Subtract(mset1,GhostCells);
      else if (ptype == GHOST)
        mset = MSets_Subtract(mset1,OwnedCells);
      else
        mset = mset1;      
    }
    

    /* Check if there were no entities left on any processor after
       extracting the appropriate category of entities */
    
    nent_loc = mset ? MSet_Num_Entries(mset) : 0;
    
    epcomm->SumAll(&nent_loc,&nent_glob,1);
    
    if (nent_glob == 0) {
      std::stringstream stream;
      stream << "Could not retrieve any mesh entities for set " << setname << std::endl;
      Errors::Message mesg(stream.str());
      Exceptions::amanzi_throw(mesg);
    }
    
    /* If this processor has no entities, nothing else to do */
    
    if (mset) {
      idx = 0;
      while ((ment = MSet_Next_Entry(mset,&idx))) {
        lid = MEnt_ID(ment);
        setents->push_back(lid-1);
      }
      
      if (ptype != USED)
        MSet_Delete(mset);
    }

    return;
    break;
  }
  case FACE: {

    if (nent_loc) {
      if (ptype == OWNED)
        mset = MSets_Subtract(mset1,NotOwnedFaces);
      else if (ptype == GHOST)
        mset = MSets_Subtract(mset1,OwnedFaces);
      else
        mset = mset1;
    }

    /* Check if there were no entities left on any processor after
       extracting the appropriate category of entities */
 
    nent_loc = mset ? MSet_Num_Entries(mset) : 0;

    epcomm->SumAll(&nent_loc,&nent_glob,1);

    if (nent_glob == 0) {
      std::stringstream stream;
      stream << "Could not retrieve any mesh entities for set " << setname << std::endl;
      Errors::Message mesg(stream.str());
      Exceptions::amanzi_throw(mesg);
    }

    /* If this processor has no entities, nothing else to do */

    if (mset) {
      idx = 0;
      while ((ment = MSet_Next_Entry(mset,&idx))) {
        lid = MEnt_ID(ment);
        setents->push_back(lid-1);
      }

      if (ptype != USED)
        MSet_Delete(mset);
    }
      
    return;
    break;
  }
  case NODE: {

    if (nent_loc) {    
      if (ptype == OWNED)
        mset = MSets_Subtract(mset1,NotOwnedVerts);
      else if (ptype == GHOST)
        mset = MSets_Subtract(mset1,OwnedVerts);
      else
        mset = mset1;
    }
      
    /* Check if there were no entities left on any processor after
       extracting the appropriate category of entities */
 
    nent_loc = mset ? MSet_Num_Entries(mset) : 0;

    epcomm->SumAll(&nent_loc,&nent_glob,1);

    if (nent_glob == 0) {
      std::stringstream stream;
      stream << "Could not retrieve any mesh entities for set " << setname << std::endl;
      Errors::Message mesg(stream.str());
      Exceptions::amanzi_throw(mesg);
    }

    /* If this processor has no entities, nothing else to do */

    if (mset) {
      idx = 0;
      while ((ment = MSet_Next_Entry(mset,&idx))) {
        lid = MEnt_ID(ment);
        setents->push_back(lid-1);
      }

      if (ptype != USED)
        MSet_Delete(mset);
    }
      
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



// Procedure to perform all the post-mesh creation steps in a constructor

void Mesh_MSTK::post_create_steps_()
{


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
  
  if (Mesh::geometric_model() != NULL)
    init_set_info();

}


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

    if (cell_dimension() == 3) {
      attfc0 = MAttrib_New(mesh,"TMP_FC0_ATT",INT,MFACE);
      attfc1 = MAttrib_New(mesh,"TMP_FC1_ATT",INT,MFACE);
    }
    else if (cell_dimension() == 2) {
      attfc0 = MAttrib_New(mesh,"TMP_FC0_ATT",INT,MEDGE);
      attfc1 = MAttrib_New(mesh,"TMP_FC1_ATT",INT,MEDGE);
    }

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
	      if (MF_EdgeDir(face0,edge) == 1)     // Sanity check
		MEnt_Set_AttVal(edge,attfc0,MEnt_GlobalID(face0),0.0,NULL);
	      else
		std::cerr << "Two faces using edge in same direction in 2D mesh" << std::endl;
	    }
	  }
	  else {
	    MEnt_Set_AttVal(edge,attfc0,MEnt_GlobalID(face0),0.0,NULL);
            face1 = List_Entry(efaces,1);
            if (face1)
              MEnt_Set_AttVal(edge,attfc1,MEnt_GlobalID(face1),0.0,NULL);
          }
          List_Delete(efaces);
	}
      }

    }  // else if (cell_dimension() == 2)



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
      
	MEnt_Get_AttVal(edge,attfc0,&remote_faceid0,&rval,&pval);
	MEnt_Get_AttVal(edge,attfc1,&remote_faceid1,&rval,&pval);
      
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

    if (mset) {

      // The only time a mesh set by this name should already exist is
      // if the region is of type labeled set. In that case, verify
      // that the entity types in the set are the same as the one
      // requested in the region
      
      if (rgn->type() == AmanziGeometry::LABELEDSET) {

        AmanziGeometry::LabeledSetRegionPtr lsrgn =
          dynamic_cast<AmanziGeometry::LabeledSetRegionPtr> (rgn);
        
        MType entdim = MSet_EntDim(mset);
        if (Mesh::space_dimension() == 3) {

          if ((lsrgn->entity_str() == "CELL" && entdim != MREGION) ||
              (lsrgn->entity_str() == "FACE" && entdim != MFACE) ||
              (lsrgn->entity_str() == "NODE" && entdim != MVERTEX)) {
            std::cerr << "Mismatch of entity type in labeled set region and mesh set" << std::endl;
            throw std::exception();
          }
        }
        else if (Mesh::space_dimension() == 2) {
          if ((lsrgn->entity_str() == "CELL" && entdim != MFACE) ||
              (lsrgn->entity_str() == "FACE" && entdim != MEDGE) ||
              (lsrgn->entity_str() == "NODE" && entdim != MVERTEX)) {
            std::cerr << "Mismatch of entity type in labeled set region and mesh set" << std::endl;
            throw std::exception();
          }
        }
      }
    }
    else { 
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

          List_Delete(vfaces);
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



int Mesh_MSTK::generate_regular_mesh(Mesh_ptr mesh, double x0, double y0, 
				     double z0, double x1, double y1, 
				     double z1, int nx, int ny, int nz)
{
/*

  Index directions for classification templates

  k   j
  |  /
  | /
  |/___ i


  Model vertex, edge and face enumeration for classification templates 


         MODEL                   MODEL                  MODEL
         VERTICES                EDGES                  FACES

     7 ____________ 8          ______7_____           ____________  
      /|          /|          /|          /|         /|      2   /| 
     / |         / |       12/ |8      11/ | 	    / |  4      / | 
   5/___________/6 |        /_____3_____/  |6	   /___________/  | 
    |  |        |  |        |  |        |  | 	   |  |        | 5| 
    |  |________|__|        |  |_____5__|__| 	   |6 |_1______|__| 
    |  /3       |  /4      4|  /        |  / 	   |  /        |  / 
    | /         | /         | /9       2| /10	   | /      3  | /  
    |/__________|/          |/__________|/   	   |/__________|/   
   1             2                1
                                                   
                                                    Front  - Face1
						    Back   - Face2
						    Bottom - Face3
						    Top    - Face4
						    Left   - Face6
						    Right  - Face5

  Classification of mesh regions onto multiple material regions is not done
  here since the "geometric model" could have overlapping regions. Instead
  mesh sets are created as necessary based on point location in regions.

*/

  int i, j, k, ii, jj, kk, gid, gdim;
  double xyz[3], dx, dy, dz;
  MVertex_ptr ***verts, mv, rverts[8], fverts[4], everts[2];
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  int vgid_tmpl[3][3][3] = {{{1,4,5},{9,6,12},{3,8,7}},{{1,1,3},{3,1,4},{5,2,7}},{{2,2,6},{10,5,11},{4,6,8}}};
  int vgdim_tmpl[3][3][3]= {{{0,1,0},{1,2,1}, {0,1,0}},{{1,2,1},{2,3,2},{1,2,1}},{{0,1,0},{1,2,1},{0,1,0}}};
  int egdim_tmpl[3][3] = {{1,2,1},{2,3,2},{1,2,1}};
  int egid_tmpl2[3][3] = {{4,6,8},{1,1,2},{2,5,6}};  /* Y direction edges (iterating over i,k) */
  int egid_tmpl1[3][3] = {{9,6,12},{3,1,4},{10,5,11}}; /* Z direction edges (iterating over i,j)*/
  int egid_tmpl0[3][3] = {{1,1,3},{3,1,4},{5,2,7}}; /* X direction edges (iterating over j,k) */
  int fgdim_tmpl[3] = {2,3,2};
  int fgid_tmpl0[3] = {6,1,5};
  int fgid_tmpl1[3] = {3,1,4};
  int fgid_tmpl2[3] = {1,1,2};

  dx = (x1-x0)/nx;
  dy = (y1-y0)/ny;
  dz = (z1-z0)/nz;

  verts = (MVertex_ptr ***) malloc((nx+1)*sizeof(MVertex_ptr **));
  for (j = 0; j < nx+1; j++) {
    verts[j] = (MVertex_ptr **) malloc((ny+1)*sizeof(MVertex_ptr *)); 
    for (k = 0; k < ny+1; k++) 
      verts[j][k] = (MVertex_ptr *) malloc((nz+1)*sizeof(MVertex_ptr));
  }

  for (k = 0; k < nz+1; k++) {
    xyz[2] = (k == nz) ? z1 : z0 + k*dz;
    kk =  (k%nz) ? 1 : (k ? 2 : 0);

    for (j = 0; j < ny+1; j++) {
      xyz[1] = (j == ny) ? y1 : y0 + j*dy;      
      jj = (j%ny) ? 1 : (j ? 2 : 0);

      for (i = 0; i < nx+1; i++) {
	xyz[0] = (i == nx) ? x1 : x0 + i*dx;
	ii = (i%nx) ? 1 : (i ? 2 : 0);
	
	mv = MV_New(mesh);
	MV_Set_Coords(mv,xyz);	
	verts[i][j][k] = mv;

	gdim  = vgdim_tmpl[ii][jj][kk];
	MV_Set_GEntDim(mv,gdim);

	gid = vgid_tmpl[ii][jj][kk];
	MV_Set_GEntID(mv,gid);
      }
    }
  }


  /* Create the edges explicitly to get the classification right */
  for (i = 0; i < nx+1; i++) {
    for (j = 0; j < ny+1; j++) {
      for (k = 0; k < nz; k++) {
	me = ME_New(mesh);

	everts[0] = verts[i][j][k];
	everts[1] = verts[i][j][k+1];
	ME_Set_Vertex(me,0,everts[0]);
	ME_Set_Vertex(me,1,everts[1]);

	ii = (i%nx) ? 1 : (i ? 2 : 0);
	jj = (j%ny) ? 1 : (j ? 2 : 0);
	gdim = egdim_tmpl[ii][jj];
	gid = egid_tmpl2[ii][jj];

	ME_Set_GEntDim(me,gdim);
	ME_Set_GEntID(me,gid);
      }
    }
  }
	
  for (i = 0; i < nx+1; i++) {
    for (k = 0; k < nz+1; k++) {
      for (j = 0; j < ny; j++) {
	me = ME_New(mesh);

	everts[0] = verts[i][j][k];
	everts[1] = verts[i][j+1][k];
	ME_Set_Vertex(me,0,everts[0]);
	ME_Set_Vertex(me,1,everts[1]);

	ii = (i%nx) ? 1 : (i ? 2 : 0);
	kk = (k%nz) ? 1 : (k ? 2 : 0);
	gdim = egdim_tmpl[ii][kk];
	gid = egid_tmpl1[ii][kk];

	ME_Set_GEntDim(me,gdim);
	ME_Set_GEntID(me,gid);
      }
    }
  }
	
  for (j = 0; j < ny+1; j++) {
    for (k = 0; k < nz+1; k++) {
      for (i = 0; i < nx; i++) {
	me = ME_New(mesh);

	everts[0] = verts[i][j][k];
	everts[1] = verts[i+1][j][k];
	ME_Set_Vertex(me,0,everts[0]);
	ME_Set_Vertex(me,1,everts[1]);

	jj = (j%ny) ? 1 : (j ? 2 : 0);
	kk = (k%nz) ? 1 : (k ? 2 : 0);
	gdim = egdim_tmpl[jj][kk];
	gid = egid_tmpl0[jj][kk];

	ME_Set_GEntDim(me,gdim);
	ME_Set_GEntID(me,gid);
      }
    }
  }
	


  /* Create the faces explicitly to get the classification right */
  for (i = 0; i < nx+1; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
	mf = MF_New(mesh);

	fverts[0] = verts[i][j][k];
	fverts[1] = verts[i][j+1][k];
	fverts[2] = verts[i][j+1][k+1];
	fverts[3] = verts[i][j][k+1];
	MF_Set_Vertices(mf,4,fverts);

	ii = (i%nx) ? 1 : (i ? 2 : 0);
	gdim = fgdim_tmpl[ii];
	gid = fgid_tmpl0[ii];

	MF_Set_GEntDim(mf,gdim);
	MF_Set_GEntID(mf,gid);
      }
    }
  }
	
  for (j = 0; j < ny+1; j++) {
    for (i = 0; i < nx; i++) {
      for (k = 0; k < nz; k++) {
	mf = MF_New(mesh);

	fverts[0] = verts[i][j][k];
	fverts[1] = verts[i+1][j][k];
	fverts[2] = verts[i+1][j][k+1];
	fverts[3] = verts[i][j][k+1];
	MF_Set_Vertices(mf,4,fverts);

	jj = (j%ny) ? 1 : (j ? 2 : 0);
	gdim = fgdim_tmpl[jj];
	gid = fgid_tmpl1[jj];

	MF_Set_GEntDim(mf,gdim);
	MF_Set_GEntID(mf,gid);
      }
    }
  }
	
  for (k = 0; k < nz+1; k++) {
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
	mf = MF_New(mesh);

	fverts[0] = verts[i][j][k];
	fverts[1] = verts[i+1][j][k];
	fverts[2] = verts[i+1][j+1][k];
	fverts[3] = verts[i][j+1][k];
	MF_Set_Vertices(mf,4,fverts);

	kk = (k%nz) ? 1 : (k ? 2 : 0);
	gdim = fgdim_tmpl[kk];
	gid = fgid_tmpl2[kk];

	MF_Set_GEntDim(mf,gdim);
	MF_Set_GEntID(mf,gid);
      }
    }
  }
	

  /* Not the most efficient way but the easiest to code */

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
	mr = MR_New(mesh);
	MR_Set_GEntID(mr,1);
	
	rverts[0] = verts[i][j][k];       rverts[1] = verts[i+1][j][k]; 
	rverts[2] = verts[i+1][j+1][k];   rverts[3] = verts[i][j+1][k];
	rverts[4] = verts[i][j][k+1];     rverts[5] = verts[i+1][j][k+1]; 
	rverts[6] = verts[i+1][j+1][k+1]; rverts[7] = verts[i][j+1][k+1];

	MR_Set_Vertices(mr, 8, rverts, 6, NULL);
      }
    }
  }
      
  for (i = 0; i < nx+1; i++) {
    for (j = 0; j < ny+1; j++)
      free(verts[i][j]);
    free(verts[i]);
  }
  free(verts);

  return 1;
}


int Mesh_MSTK::generate_regular_mesh(Mesh_ptr mesh, double x0, double y0, 
				     double x1, double y1, int nx, int ny)
{
  int i, j, dir[4];
  double xyz[3], llx, lly, urx, ury, dx, dy;
  MVertex_ptr **verts, v0, v1, mv;
  MEdge_ptr fedges[4], me;
  MFace_ptr mf;

  dx = (x1-x0)/nx;
  dy = (y1-y0)/ny;

  verts = (MVertex_ptr **) malloc((nx+1)*sizeof(MVertex_ptr *));
  for (i = 0; i < nx+1; i++)
    verts[i] = (MVertex_ptr *) malloc((ny+1)*sizeof(MVertex_ptr));
 
  xyz[2] = 0.0;
  for (j = 0; j < ny+1; j++) {
    xyz[1] = (j == ny) ? y1 : y0 + j*dy;

    for (i = 0; i < nx+1; i++) {
      xyz[0] = (i == nx) ? x1 : x0 + i*dx;

      mv = MV_New(mesh);
      MV_Set_Coords(mv,xyz);

      if (i == 0) {
	if (j == 0) {
	  MV_Set_GEntDim(mv,0);
	  MV_Set_GEntID(mv,1);
	}
	else if (j == ny) {
	  MV_Set_GEntDim(mv,0);
	  MV_Set_GEntID(mv,4);	  
	}
	else {
	  MV_Set_GEntDim(mv,1);
	  MV_Set_GEntID(mv,4);
	}
      }
      else if (i == nx) {
	if (j == 0) {
	  MV_Set_GEntDim(mv,0);
	  MV_Set_GEntID(mv,2);
	}
	else if (j == ny) {
	  MV_Set_GEntDim(mv,0);
	  MV_Set_GEntID(mv,3);	  
	}
	else {
	  MV_Set_GEntDim(mv,1);
	  MV_Set_GEntID(mv,2);
	}
      }
      else {
	if (j == 0) {
	  MV_Set_GEntDim(mv,1);
	  MV_Set_GEntID(mv,1);
	}
	else if (j == ny) {
	  MV_Set_GEntDim(mv,1);
	  MV_Set_GEntID(mv,3);
	}
	else {
	  MV_Set_GEntDim(mv,2);
	  MV_Set_GEntID(mv,1);
	}
      }

      verts[i][j] = mv;
    }
  }


  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      mf = MF_New(mesh);
      
      /* edge 0 */
      v0 = verts[i][j];
      v1 = verts[i+1][j];
      fedges[0] = MVs_CommonEdge(v0,v1);
      if (fedges[0])
	dir[0] = (ME_Vertex(fedges[0],0) == v0) ? 1 : 0;
      else {
	me = ME_New(mesh);
	
	ME_Set_Vertex(me,0,v0);
	ME_Set_Vertex(me,1,v1);
	
	if (j == 0) {
	  ME_Set_GEntDim(me,1);
	  ME_Set_GEntID(me,1);
	}
	else {
	  ME_Set_GEntDim(me,2);
	  ME_Set_GEntID(me,1);
	}
	
	fedges[0] = me;
	dir[0] = 1;
      }
      
      
      /* edge 1 */
      v0 = verts[i+1][j];
      v1 = verts[i+1][j+1];
      fedges[1] = MVs_CommonEdge(v0,v1);
      if (fedges[1])
	dir[1] = (ME_Vertex(fedges[1],0) == v0) ? 1 : 0;
      else {
	me = ME_New(mesh);
	
	ME_Set_Vertex(me,0,v0);
	ME_Set_Vertex(me,1,v1);
	
	if (i+1 == nx) {
	  ME_Set_GEntDim(me,1);
	  ME_Set_GEntID(me,2);
	}
	else {
	  ME_Set_GEntDim(me,2);
	  ME_Set_GEntID(me,1);
	}
	
	fedges[1] = me;
	dir[1] = 1;
      }
      
      
      /* edge 2 */
      v0 = verts[i+1][j+1];
      v1 = verts[i][j+1];
      fedges[2] = MVs_CommonEdge(v0,v1);
      if (fedges[2])
	dir[2] = (ME_Vertex(fedges[2],0) == v0) ? 1 : 0;
      else {
	me = ME_New(mesh);
	
	ME_Set_Vertex(me,0,v0);
	ME_Set_Vertex(me,1,v1);
	
	if (j+1 == nx) {
	  ME_Set_GEntDim(me,1);
	  ME_Set_GEntID(me,3);
	}
	else {
	  ME_Set_GEntDim(me,2);
	  ME_Set_GEntID(me,1);
	}
	
	fedges[2] = me;
	dir[2] = 1;
      }
      
      
      /* edge 3 */
      v0 = verts[i][j+1];
      v1 = verts[i][j];
      fedges[3] = MVs_CommonEdge(v0,v1);
      if (fedges[3])
	dir[3] = (ME_Vertex(fedges[3],0) == v0) ? 1 : 0;
      else {
	me = ME_New(mesh);
	
	ME_Set_Vertex(me,0,v0);
	ME_Set_Vertex(me,1,v1);
	
	if (i == 0) {
	  ME_Set_GEntDim(me,1);
	  ME_Set_GEntID(me,4);
	}
	else {
	  ME_Set_GEntDim(me,2);
	  ME_Set_GEntID(me,1);
	}
	
	fedges[3] = me;
	dir[3] = 1;
      }


      MF_Set_Edges(mf,4,fedges,dir);

      MF_Set_GEntDim(mf,2);
      MF_Set_GEntID(mf,1);
    }
  }
   
  for (i = 0; i < nx+1; i++)
    free(verts[i]);
  free(verts);

  return 1;
}

void Mesh_MSTK::pre_create_steps_(const int space_dimension, 
                                  const Epetra_MpiComm *comm, 
                                  const AmanziGeometry::GeometricModelPtr& gm) 
{
  clear_internals_();

  MSTK_Init();

  Mesh::set_comm(comm);
  Mesh::set_geometric_model(gm);

  set_space_dimension(space_dimension);

  MPI_Comm_rank(mpicomm,&myprocid);
  MPI_Comm_size(mpicomm,&numprocs);

  serial_run =  (!mpicomm || numprocs == 1) ? true : false;

}

} // close namespace AmanziMesh
} // close namespace Amanzi
