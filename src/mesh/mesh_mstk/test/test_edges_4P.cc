#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"
#include "MeshAudit.hh"

#include "Epetra_Map.h"
#include "AmanziComm.hh"

// Test edge functions in 2D

TEST(MSTK_EDGES_2D)
{
  auto comm = Amanzi::getDefaultComm();
  int size = comm->NumProc();
  CHECK_EQUAL(4,size);
  
  int DebugWait = 0;
  while (DebugWait);

  // Generate a 4x4 quad mesh distributed over four processors
  
  bool request_faces = true, request_edges = true;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
    mesh(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,2.0,1.0,4,4,comm,Teuchos::null,
					   Teuchos::null,request_faces,request_edges));

  // Check that we get the expected number of edges

  int ne_owned = mesh->num_entities(Amanzi::AmanziMesh::EDGE,
				    Amanzi::AmanziMesh::Parallel_type::OWNED);
  int ne_all = mesh->num_entities(Amanzi::AmanziMesh::EDGE,
				  Amanzi::AmanziMesh::Parallel_type::ALL);
  CHECK(ne_owned <= ne_all);

  // This assumes a symmetric partitioning - not always the case with
  // ZOLTAN graph partitioning

  //  CHECK_EQUAL(24,ne_all);

  // In 2D, faces and edges are the same - so face global IDs and edge
  // global IDs for a cell must match

  int nc_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL,
				    Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < nc_owned; ++c) {
    Amanzi::AmanziMesh::Entity_ID_List cedges, cfaces, fedges;
    std::vector<int> cfdirs, fedirs, cedirs;    

    mesh->cell_2D_get_edges_and_dirs(c,&cedges,&cedirs);
    mesh->cell_get_faces_and_dirs(c,&cfaces,&cfdirs);

    for (int e = 0; e < cedges.size(); ++e) {
      CHECK_EQUAL(mesh->GID(cedges[e],Amanzi::AmanziMesh::EDGE), 
		  mesh->GID(cfaces[e],Amanzi::AmanziMesh::FACE));

      // Also, see if the direction and vector we got for edges of 2D
      // cell is consistent with the direction and normal vector we
      // got for the faces of the cell
      
      CHECK_EQUAL(cedirs[e],cfdirs[e]);

      Amanzi::AmanziGeometry::Point evec(2), fnormal(2), ftangent(2);

      evec = mesh->edge_vector(cedges[e])*cedirs[e];

      fnormal = mesh->face_normal(cfaces[e])*cfdirs[e];
      ftangent.set(-fnormal[1],fnormal[0]);

      CHECK_EQUAL(evec[0],ftangent[0]);
      CHECK_EQUAL(evec[1],ftangent[1]);
    }


    for (int f = 0; f < cfaces.size(); ++f) {
      mesh->face_get_edges_and_dirs(cfaces[f],&fedges,&fedirs);

      CHECK_EQUAL(1,fedges.size()); // face is same as edge in 2D
      CHECK_EQUAL(1,fedirs[0]); // direction is always 1
      
      // check the face-edges to cell-edges map

      std::vector<int> map;

      mesh->face_to_cell_edge_map(cfaces[f],c,&map);

      for (int e = 0; e < fedges.size(); ++e)
	CHECK_EQUAL(fedges[e],cedges[map[e]]);
    }
  }

  // owing to how we constructed the mesh, the length of horizontal edges 
  // should be 0.5 and vertical edges 0.25

  for (int e = 0; e < ne_owned; ++e) {
    Amanzi::AmanziGeometry::Point evec(2);
    double elen;

    evec = mesh->edge_vector(e);
    elen = mesh->edge_length(e);
    if (evec[1] == 0.0) {
      CHECK_EQUAL(0.5,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
    else if (evec[0] == 0.0) {
      CHECK_EQUAL(0.25,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
  }
}



// Test edge functions in 3D

TEST(MSTK_EDGES_3D)
{
  auto comm = Amanzi::getDefaultComm();
  int size = comm->NumProc();
  CHECK_EQUAL(4,size);

  //  if (rank == 0) {
  int DebugWait = 0;
  while (DebugWait);
  //  }

  // Generate a 4x4x4 quad mesh distributed over four processors
  
  bool request_faces = true, request_edges = true;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
    mesh(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,2.0,1.0,4.0,4,4,4,
					   comm,Teuchos::null,Teuchos::null,request_faces,
					   request_edges));

  // How many owned and used edges are there?

  int ne_owned = mesh->num_entities(Amanzi::AmanziMesh::EDGE,
				    Amanzi::AmanziMesh::Parallel_type::OWNED);
  int ne_all = mesh->num_entities(Amanzi::AmanziMesh::EDGE,
				  Amanzi::AmanziMesh::Parallel_type::ALL);

  // Check that we got a non-zero number

  CHECK(ne_owned != 0);
  CHECK(ne_all != 0);  


  // Go through the cells and retrieve their edges to make sure it
  // works correctly. Also, get the faces of the cells and the edges
  // of these faces and do additional checks

  int nc_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL,
				    Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < nc_owned; ++c) {
    Amanzi::AmanziMesh::Entity_ID_List cedges, cfaces, fedges;
    std::vector<int> cfdirs, fedirs;    

    mesh->cell_get_edges(c,&cedges);
    mesh->cell_get_faces_and_dirs(c,&cfaces,&cfdirs);

    for (int f = 0; f < cfaces.size(); ++f) {
      mesh->face_get_edges_and_dirs(cfaces[f],&fedges,&fedirs);

      // check the face-edges to cell-edges map

      std::vector<int> map;
      mesh->face_to_cell_edge_map(cfaces[f],c,&map);

      for (int e = 0; e < fedges.size(); ++e)
	CHECK_EQUAL(fedges[e],cedges[map[e]]);
    }
  }

  // owing to how we constructed the mesh, the length of x-direction
  // should be 0.5, y-direction edges should
  // 0.25 and z-direction edges should be 1.0

  for (int e = 0; e < ne_owned; ++e) {
    Amanzi::AmanziGeometry::Point evec(2);
    double elen;

    evec = mesh->edge_vector(e);
    elen = mesh->edge_length(e);
    if (evec[0] != 0.0 && evec[1] == 0.0) {
      CHECK_EQUAL(0.5,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
    else if (evec[0] == 0.0 && evec[1] != 0.0) {
      CHECK_EQUAL(0.25,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
    else  {
      CHECK_EQUAL(1.0,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
  }

}

