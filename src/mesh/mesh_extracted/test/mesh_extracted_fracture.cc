/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include <UnitTest++.h>

// Amanzi::Mesh
#include "MeshAudit.hh"
#include "MeshExtractedManifold.hh"
#ifdef HAVE_MESH_MOAB
#include "Mesh_MOAB.hh"
#endif
#ifdef HAVE_MESH_MSTK
#include "Mesh_MSTK.hh"
#endif
#include "Mesh_simple.hh"

/* **************************************************************** */
void RunTest(const std::string regname, int* cells, int* edges) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read parameter list
  std::string xmlFileName = "test/mesh_extracted_fracture.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  std::string setname(regname);

  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  auto mesh_list = Teuchos::sublist(plist, "mesh", false);

  for (int i = 0; i < 3; ++i) {
    RCP<const Mesh> mesh3D;
    if (i == 0) {
#ifdef HAVE_MESH_MSTK
      std::cout << "\nMesh framework: MSTK (" << regname << ")\n";
      // mesh3D = Teuchos::rcp(new Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, comm, gm, mesh_list, true, true));
      mesh3D = Teuchos::rcp(new Mesh_MSTK("test/mesh_extracted_fracture.exo", comm, gm, mesh_list, true, true));
#endif
    } else if (i == 1 && comm->NumProc() == 1) {
      std::cout << "\nMesh framework: simple\n";
      mesh3D = Teuchos::rcp(new Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, comm, gm, mesh_list, true, true));
    } else if (i == 2) {
#ifdef HAVE_MESH_MOAB
      if (comm->NumProc() > 1) continue;
      std::cout << "\nMesh framework: MOAB\n";
      mesh3D = Teuchos::rcp(new Mesh_MOAB("test/mesh_extracted_fracture.exo", comm, gm, mesh_list, true, true));
#endif
    }
    if (mesh3D == Teuchos::null) continue;

    // extract fractures mesh
    try {
      RCP<Mesh> mesh = Teuchos::rcp(new MeshExtractedManifold(mesh3D, setname, AmanziMesh::FACE,
                                                              comm, gm, plist, true, false));

      int ncells = mesh->cell_map(false).NumGlobalElements();
      int nfaces = mesh->face_map(false).NumGlobalElements();
      int mfaces = (i == 0) ? mesh->exterior_face_map(false).NumGlobalElements() : 0;
      std::cout << i << " pid=" << comm->MyPID() << " cells: " << ncells 
                     << " faces: " << nfaces << " bnd faces: " << mfaces << std::endl;
      CHECK(cells[i] == ncells);
      CHECK(edges[i] == mfaces);

      // verify mesh 
      MeshAudit audit(mesh);
      int ok = audit.Verify();
      CHECK(ok == 0);

    } catch (...) {
      std::cout << "Framework failed.\n";
    }
  }
}


TEST(MESH_EXTRACTED_FRACTURE_NETWORK2) {
  int cells[3] = {108, 200, 0};
  int edges[3] = {48, 0, 0};
  RunTest("fractures-two", cells, edges);
}

TEST(MESH_EXTRACTED_SURFACE) {
  int cells[3] = {9, 25, 0};
  int edges[3] = {12, 0, 0};
  RunTest("Left side", cells, edges);
}

TEST(MESH_EXTRACTED_FRACTURE_NETWORK3) {
  int cells[3] = {108, 200, 0};
  int edges[3] = {72, 0, 0};
  RunTest("fractures-three", cells, edges);
}

