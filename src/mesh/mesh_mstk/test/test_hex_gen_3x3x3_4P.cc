#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"

#include "Epetra_Map.h"
#include "AmanziComm.hh"

#include "MeshAudit.hh"


// Test for generation of hex mesh distributed over 4 processors

TEST(MSTK_HEX_GEN_3x3x3_4P)
{
  std::vector<Amanzi::AmanziMesh::Entity_ID> faces(6), nodes(8);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  auto comm = Amanzi::getDefaultComm();
  int rank = comm->MyPID();
  int size = comm->NumProc();
  CHECK_EQUAL(4,size);
  if (size != 4) {
    std::cerr << "Test must be run with 4 processors" << std::endl;
    //    return;
  }

  // if (rank == 0) {
  int DebugWait = 0;
  while (DebugWait);
  // }

  // Create a 3x3x3 cell hex mesh

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,1.0,1.0,1.0,3,3,3,comm));


  std::vector<Amanzi::AmanziMesh::Entity_ID>  c2f(6);
  Epetra_Map cell_map(mesh->cell_map(false));
  Epetra_Map face_map(mesh->face_map(true));

  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      CHECK_EQUAL(cell_map.GID(c),mesh->GID(c,Amanzi::AmanziMesh::CELL));
      mesh->cell_get_faces(c, &c2f, true);

      for (int j=0; j<6; j++)
	{
	  int f = face_map.LID(mesh->GID(c2f[j],Amanzi::AmanziMesh::FACE));
	  CHECK_EQUAL( f,c2f[j] );
	  if (f != c2f[j]) {
	    std::cout << std::endl;
	    std::cout << "Processor ID " << rank << std::endl;
	    std::cout << "Cell ID " << cell_map.GID(c) << std::endl;
	    std::cout << "Problem face c2f[j] = " << c2f[j] << " GID = " << mesh->GID(c2f[j],Amanzi::AmanziMesh::FACE) << " f = " << f << std::endl;
	    std::cout << std::endl;
	  }
	}

    }

  std::stringstream fname;
  fname << "test/mstk_hex_gen_3x3x3_4P." << rank << ".out";
  std::ofstream fout(fname.str().c_str());
  Amanzi::MeshAudit auditor(mesh,fname);
  auditor.Verify();  

}

