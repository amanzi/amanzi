#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"

#include "MeshAudit.hh"

#include "AmanziMap.hh"
#include "AmanziComm.hh"


TEST(MSTK_HEX_3x3x3)
{

  int i, j, k, err, nc, nf, nv;
  Amanzi::AmanziMesh::Entity_ID faces[6], nodes[8];

  int NV = 64;
  int NF = 108;
  int NC = 27;

  auto comm = Amanzi::getDefaultComm();

  // Load a mesh consisting of 3x3x3 elements

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_MSTK("test/hex_3x3x3_sets.exo",comm));

  nf = mesh->num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NF,nf);

  nc = mesh->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC,nc);


  std::vector<Amanzi::AmanziMesh::Entity_ID>  c2f(6);
  auto cell_map = mesh->cell_map(false);
  auto face_map = mesh->face_map(false);
  for (int c=cell_map->getMinLocalIndex(); c<=cell_map->getMaxLocalIndex(); c++)
    {
      CHECK_EQUAL(cell_map->getGlobalElement(c),mesh->getGlobalElement(c,Amanzi::AmanziMesh::CELL));
      mesh->cell_get_faces(c, &c2f, true);
      for (int j=0; j<6; j++)
	{
	  int f = face_map->getLocalElement(c2f[j]);
	  CHECK( f == c2f[j] );
	}

    }

  std::stringstream fname;
  fname << "test/mstk_hex_3x3x3.out";
  std::ofstream fout(fname.str().c_str());
  Amanzi::MeshAudit auditor(mesh,fout);
  auditor.Verify();

}
