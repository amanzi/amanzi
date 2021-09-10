/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include <UnitTest++.h>

#include <iostream>

#include "AmanziComm.hh"
#include "Geometry.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshException.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"

using namespace Amanzi;

TEST(MESH_DEFORM2D)
{
  auto comm = getDefaultComm();

  // We are not including MOAB or SIMPLE since they cannot generate in 2D
  std::vector<AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

  if (AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(AmanziMesh::Framework::MSTK);
    framework_names.push_back("MSTK");
  }

  for (const auto& frm : frameworks) {
    // Set the framework
    std::cout << "Testing deformation with " << AmanziMesh::framework_names.at(frm) << std::endl;

    // Create the mesh
    Teuchos::RCP<AmanziMesh::Mesh> mesh = createFrameworkStructuredUnitQuad({frm}, 10, 10, comm,
            Teuchos::null, Teuchos::null, false, 10.0, 10.0);

    // Deform the mesh
    AmanziMesh::Entity_ID_List nodeids;
    AmanziGeometry::Point_List newpos, finpos;

    int nnodes = mesh->num_entities(AmanziMesh::NODE,
                                AmanziMesh::Parallel_type::OWNED);
    for (int j = 0; j < nnodes; j++) {
      nodeids.push_back(j);
      AmanziGeometry::Point oldcoord(2),newcoord(2);
      mesh->node_get_coordinates(j,&oldcoord);
      newcoord.set(oldcoord[0],0.5*oldcoord[1]);
      newpos.push_back(newcoord);
    }

    int status = mesh->deform(nodeids,newpos,false,&finpos);
    CHECK_EQUAL(status,1);

    // If the deformation was successful, the cell volumes should be half
    // of what they were
    int ncells = mesh->num_entities(AmanziMesh::CELL,
                                    AmanziMesh::Parallel_type::ALL);

    for (int j = 0; j < ncells; j++) {
      double volume = mesh->cell_volume(j);
      CHECK_CLOSE(0.5, volume, 1.e-10);
    }
  } // for each framework i
}


TEST(MESH_GENERATED_DEFORM3D)
{
  auto comm = getDefaultComm();
  const int nproc(comm->NumProc());
  if (nproc != 1) {
    std::cout << "Parallel deformation not implemented" << std::endl;
    return;
  }

  std::vector<AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

  frameworks.push_back(AmanziMesh::Framework::SIMPLE);
  framework_names.push_back("simple");

  if (AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(AmanziMesh::Framework::MSTK);
    framework_names.push_back("MSTK");
  }


  for (const auto& frm : frameworks) {
    std::cout << "Testing deformation with " << AmanziMesh::framework_names.at(frm) << std::endl;

    // start with a mesh that will be deformed into the known mesh coordinates
    auto mesh = createFrameworkStructuredUnitHex({frm}, 3, 3, 3,
            comm, Teuchos::null, Teuchos::null, false, 1.0, 1.0, 2.0);

    // Deform the mesh
    AmanziMesh::Entity_ID_List nodeids;
    AmanziGeometry::Point_List newpos, finpos;

    int status = 0;
    int nnodes = mesh->num_entities(AmanziMesh::NODE,
            AmanziMesh::Parallel_type::OWNED);
    for (int j = 0; j < nnodes; j++) {
      nodeids.push_back(j);

      AmanziGeometry::Point oldcoord(3),newcoord(3);
      mesh->node_get_coordinates(j,&oldcoord);

      newcoord.set(oldcoord[0], oldcoord[1], oldcoord[2]/2.0);
      newpos.push_back(newcoord);
    }
    status = mesh->deform(nodeids, newpos, false, &finpos);

    // check geometry
    testMeshAudit<MeshAudit, Mesh>(mesh);
    testGeometryCube<Mesh>(mesh, 3,3,3);
  }
}


