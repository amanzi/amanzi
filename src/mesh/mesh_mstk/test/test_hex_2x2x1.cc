#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"

#include "MeshAudit.hh"

#include "Epetra_Map.h"
#include "AmanziComm.hh"


TEST(MSTK_HEX_2x2x1)
{

  int i, j, k, err, nc, nf, nv;
  std::vector<Amanzi::AmanziMesh::Entity_ID> faces(6), cnodes(8), fnodes(6), expfacenodes(4);
  std::vector<int> facedirs(6);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  int NV = 18;
  int NF = 20;
  int NC = 4;
  double xyz[18][3] = {{-0.5,-0.5, 0.25},
		       {-0.5,-0.5,-0.25},
		       {-0.5, 0,  -0.25},
		       {-0.5, 0,   0.25},
		       { 0,  -0.5, 0.25},
		       { 0,  -0.5,-0.25},
		       { 0,   0,  -0.25},
		       { 0,   0,   0.25},
		       {-0.5, 0.5,-0.25},
		       {-0.5, 0.5, 0.25},
		       { 0,   0.5,-0.25},
		       { 0,   0.5, 0.25},
		       { 0.5,-0.5, 0.25},
		       { 0.5,-0.5,-0.25},
		       { 0.5, 0,  -0.25},
		       { 0.5, 0,   0.25},
		       { 0.5, 0.5,-0.25},
		       { 0.5, 0.5, 0.25}};

  Amanzi::AmanziMesh::Entity_ID cnstd[8] = {0,1,2,3,4,5,6,7};
  Amanzi::AmanziMesh::Entity_ID cfstd[6][4] = {{0,1,5,4},    // Expected cell-face-node pattern
			      {1,2,6,5},
			      {2,3,7,6},
			      {3,0,4,7},
			      {0,3,2,1},
			      {4,5,6,7}};

  auto comm = Amanzi::getDefaultComm(); 


  // Load a simple 4 element hex mesh 

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_MSTK("test/hex_2x2x1_ss.exo",comm));
  
  
  nv = mesh->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NV,nv);
  
  for (i = 0; i < nv; i++) {
    Amanzi::AmanziGeometry::Point coords(mesh->space_dimension());
    
    // coords.init(mesh->space_dimension()); 
    
    mesh->node_get_coordinates(i,&coords);
    CHECK_EQUAL(xyz[i][0],coords[0]);
    CHECK_EQUAL(xyz[i][1],coords[1]);
    CHECK_EQUAL(xyz[i][2],coords[2]);
  }
  
  
  nf = mesh->num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NF,nf);
  
  nc = mesh->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC,nc);
  
  for (i = 0; i < nc; i++) {
    mesh->cell_get_nodes(i,&cnodes);
    mesh->cell_get_faces_and_dirs(i,&faces,&facedirs,true);
    
    for (j = 0; j < 6; j++) {
      
      mesh->face_get_nodes(faces[j],&fnodes);
      mesh->face_get_coordinates(faces[j],&fcoords);
      

      for (k = 0; k < 4; k++)
	expfacenodes[k] = cnodes[cfstd[j][k]];

      // The order of nodes returned may be different from what we expected
      // So make sure we have a matching node to start with
      
      int k0 = -1;
      int found = 0;
      for (k = 0; k < 4; k++) {
	if (expfacenodes[k] == fnodes[0]) {
	  k0 = k;
	  found = 1;
	  break;
	}
      }
      
      CHECK_EQUAL(found,1); 
      
      if (facedirs[j] == 1) {
	for (k = 0; k < 4; k++) {
	  CHECK_EQUAL(expfacenodes[(k0+k)%4],fnodes[k]);
	  CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+k)%4]],fcoords[k],3);
	}
      }
      else {
	for (k = 0; k < 4; k++) {
	  CHECK_EQUAL(expfacenodes[(k0+4-k)%4],fnodes[k]);
	  CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+4-k)%4]],fcoords[k],3);
	}
      }
    }

  
    mesh->cell_get_coordinates(i,&ccoords);
    
    for (j = 0; j < 8; j++)
      CHECK_ARRAY_EQUAL(xyz[cnodes[cnstd[j]]],ccoords[j],3);
  }


  std::vector<Amanzi::AmanziMesh::Entity_ID>  c2f(6);
  Epetra_Map cell_map(mesh->cell_map(true));
  Epetra_Map face_map(mesh->face_map(false));

  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++) {
    CHECK_EQUAL(cell_map.GID(c),mesh->GID(c,Amanzi::AmanziMesh::CELL));
    mesh->cell_get_faces(c, &c2f);

    for (j=0; j<6; j++) {
      int f = face_map.LID(mesh->GID(c2f[j],Amanzi::AmanziMesh::FACE));
      CHECK( f == c2f[j] );
    }
  }

  std::stringstream fname;
  fname << "test/mstk_hex_2x2x1.out";
  std::ofstream fout(fname.str().c_str());
  Amanzi::MeshAudit auditor(mesh,fout);
  auditor.Verify();
}

