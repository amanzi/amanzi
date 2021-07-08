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
#include "../Mesh_MSTK.hh"

/* **************************************************************** */
TEST(EXTRACT_FRACTURE_MESH) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read parameter list
  std::string xmlFileName = "test/test_fracture_mesh.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  RCP<const Mesh> mesh3D = Teuchos::rcp(new Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10,
                                                      comm, gm, plist, true, true));

  // extract fractures mesh
  std::vector<std::string> setnames;
  setnames.push_back("fracture 1");
  setnames.push_back("fracture 2");

  Entity_ID_List ids;
  for (auto name : setnames) {
    Entity_ID_List ids_l;
    mesh3D->get_set_entities(name, AmanziMesh::FACE, Parallel_type::OWNED, &ids_l);
    ids.insert(ids.end(), ids_l.begin(), ids_l.end());
  }

  RCP<const Mesh> mesh = Teuchos::rcp(new Mesh_MSTK(mesh3D, ids, AmanziMesh::FACE,
                                                    false, comm, gm, plist, true, false));

  int ncells = mesh->cell_map(false).NumGlobalElements();
  int nfaces = mesh->face_map(false).NumGlobalElements();
  std::cout << "pid=" << comm->MyPID() << " cells: " << ncells 
                                       << " faces: " << nfaces << std::endl;
}





