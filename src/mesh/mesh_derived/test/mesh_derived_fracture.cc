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
#include "MeshDerived.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"

/* **************************************************************** */
TEST(DERIVED_MESH) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  std::vector<Framework> frameworks = {Framework::MSTK, Framework::MOAB, Framework::SIMPLE};

  // read parameter list
  std::string xmlFileName = "test/mesh_derived_fracture.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  std::string exoname("test/mesh_derived_fracture.exo");
  std::string setname("fractures");

  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  for (auto frm = frameworks.begin(); frm != frameworks.end(); ++frm) {
    // create a mesh framework
    std::cout << "\nMesh framework: " << framework_names.at(*frm) << std::endl;
    MeshFactory meshfactory(comm, gm, plist);
    meshfactory.set_preference(Preference({*frm}));
    RCP<const Mesh> mesh3D;

    if (framework_names.at(*frm) != "MOAB")
      mesh3D = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, true, true);
    else 
      mesh3D = meshfactory.create(exoname, true, true);

    if (framework_names.at(*frm) == "MSTK") mesh3D->write_to_exodus_file(exoname);

    // extract fractures mesh
    try {
      RCP<const Mesh> mesh = Teuchos::rcp(new MeshDerived(mesh3D, setname, AmanziMesh::FACE,
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


