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
#include "MeshExtractedManifold.hh"
#ifdef HAVE_MOAB_MESH
#include "Mesh_MOAB.hh"
#endif
#ifdef HAVE_MSTK_MESH
#include "Mesh_MSTK.hh"
#endif
#include "Mesh_simple.hh"

/* **************************************************************** */
TEST(MESH_EXTRACTED_FRACTURES) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read parameter list
  std::string xmlFileName = "test/mesh_extracted_fracture.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  std::string exoname("test/mesh_extracted_fracture.exo");
  std::string setname("fractures");

  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  auto mesh_list = Teuchos::sublist(plist, "mesh", false);

  for (int i = 0; i < 2; ++i) {
    RCP<const Mesh> mesh3D;
    if (i == 0) {
#ifdef HAVE_MSTK_MESH
      std::cout << "\nMesh framework: MSTK\n";
      mesh3D = Teuchos::rcp(new Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, comm, gm, mesh_list, true, true));
      mesh3D->write_to_exodus_file(exoname);
#endif
    } else if (i == 1) {
      std::cout << "\nMesh framework: simple\n";
      mesh3D = Teuchos::rcp(new Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, comm, gm, mesh_list, true, true));
    } else if (i == 2) {
#ifdef HAVE_MOAB_MESH
      std::cout << "\nMesh framework: MOAB\n";
      mesh3D = Teuchos::rcp(new Mesh_MOAB(exoname, comm, gm, mesh_list, true, true));
#endif
    }
    if (mesh3D == Teuchos::null) continue;

    // extract fractures mesh
    try {
      RCP<const Mesh> mesh = Teuchos::rcp(new MeshExtractedManifold(mesh3D, setname, AmanziMesh::FACE,
                                                                    comm, gm, plist, true, false));

      int ncells = mesh->cell_map(false).NumGlobalElements();
      int nfaces = mesh->face_map(false).NumGlobalElements();
      std::cout << "pid=" << comm->MyPID() << " cells: " << ncells 
                                           << " faces: " << nfaces << std::endl;
    } catch (...) {
      std::cout << "Framework failed.\n";
    }
  }
}


