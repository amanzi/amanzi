/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

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
#include "Mesh_MSTK.hh"

/* **************************************************************** */
TEST(MESH_SURFACE_FLATTENED) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read parameter list
  std::string xmlFileName = "test/mesh_extracted_fracture.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  auto mesh_list = Teuchos::sublist(plist, "mesh", false);
  mesh_list->set<bool>("request edges", true);
  RCP<MeshFramework> mesh3D = Teuchos::rcp(new Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, comm, gm, mesh_list));

  bool flatten = true;
  auto mesh3D_cache = Teuchos::rcp(new Mesh(mesh3D, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()))); 
  RCP<MeshFramework> mesh = Teuchos::rcp(new MeshExtractedManifold(
      mesh3D_cache, "Top side", AmanziMesh::Entity_kind::FACE, comm, gm, plist, flatten));
  
  RCP<Mesh> mesh_cached = Teuchos::rcp(new Mesh(mesh, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()))); 

  int ncells = mesh_cached->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces = mesh_cached->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  int nnodes = mesh_cached->getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
  int mfaces = mesh_cached->getNumEntities(Entity_kind::BOUNDARY_FACE, Parallel_kind::OWNED);
  std::cout << " pid=" << comm->MyPID() << " cells: " << ncells 
            << " faces: " << nfaces << " bnd faces: " << mfaces << std::endl;

  CHECK(ncells == 100);
  CHECK(nfaces == 220);
  CHECK(nnodes == 121);

  AmanziGeometry::Point p;
  for (int n = 0; n < nnodes; ++n) {
    p = mesh_cached->getNodeCoordinate(n);
    CHECK(p.dim() == 2);
    CHECK(p[0] >= 0.0 && p[1] <= 1.0);
  }

  for (int n = 0; n < nfaces; ++n) {
    p = mesh_cached->getFaceCentroid(n);
    CHECK(p.dim() == 2);
    CHECK(p[0] >= 0.0 && p[1] <= 1.0);
  }

  for (int n = 0; n < ncells; ++n) {
    p = mesh_cached->getCellCentroid(n);
    CHECK(p.dim() == 2);
    CHECK(p[0] > 0.0 && p[1] < 1.0);
  }
}
