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
#include "MeshCurved.hh"
#include "MeshFactory.hh"
#include "MeshException.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"

using namespace Amanzi;

TEST(MESH_HIGH_ORDER2D)
{
  auto comm = getDefaultComm();

  // We are not including MOAB or SIMPLE since support of high-order
  // meshes is under development
  std::vector<AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

  if (AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(AmanziMesh::Framework::MSTK);
    framework_names.push_back("MSTK");
  }

  for (const auto& frm : frameworks) {
    // Set the framework
    std::cout << "Testing high-order mesh with " << AmanziMesh::framework_names.at(frm) << std::endl;

    // Create the mesh
    Teuchos::RCP<Teuchos::ParameterList> plist;
    auto mesh = Teuchos::rcp(new AmanziMesh::MeshCurved(0.0, 0.0, 1.0, 1.0, 7, 8, comm, Teuchos::null, plist));

    int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

    // Deform the mesh
    double phi(2 * M_PI), t(0.07);
    AmanziGeometry::Point oldcoord(2), newcoord(2);

    AmanziMesh::Entity_ID_List nodeids;
    AmanziGeometry::Point_List newpos, finpos;

    for (int j = 0; j < nnodes; j++) {
      nodeids.push_back(j);
      mesh->node_get_coordinates(j, &oldcoord);
      double tmp = t * std::sin(oldcoord[0] * phi) * std::sin(oldcoord[1] * phi);
      newcoord.set(oldcoord[0] + tmp, oldcoord[1] + tmp);
      newpos.push_back(newcoord);
    }

    int status = mesh->deform(nodeids,newpos,false,&finpos);
    CHECK_EQUAL(status,1);

    // Add new nodes on mesh edges/faces
    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    auto ho_nodes_f = std::make_shared<std::vector<AmanziGeometry::Point_List> >(nfaces);

    for (int f = 0; f < nfaces; ++f) {
      oldcoord = mesh->face_centroid(f);
      double tmp = t * std::sin(oldcoord[0] * phi) * std::sin(oldcoord[1] * phi);
      newcoord.set(oldcoord[0] + tmp, oldcoord[1] + tmp);
      (*ho_nodes_f)[f].push_back(newcoord);
    }
    mesh->set_face_ho_nodes(ho_nodes_f);

    // 
    mesh->BuildCache();

    // If the deformation was successful, the total volume should be one
    double sum1(0.0), sum2(0.0), deviation(0.0);

    for (int j = 0; j < ncells; j++) {
      double volume1 = mesh->cell_volume(j);
      sum1 += volume1;

      double volume2 = mesh->cell_volume_linear(j);
      sum2 += volume2;
      deviation += std::fabs(volume1 - volume2);
    }
    CHECK_CLOSE(sum1, 1.0, 1.0e-14);
    CHECK_CLOSE(sum2, 1.0, 1.0e-14);
    CHECK(deviation > 0.1);
  }
}


