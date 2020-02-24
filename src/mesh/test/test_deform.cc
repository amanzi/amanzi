/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Rao Garimella, others
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <UnitTest++.h>

#include <iostream>

#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "Geometry.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshException.hh"

TEST(MESH_DEFORM2D)
{
  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh
  std::vector<Amanzi::AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

  if (Amanzi::AmanziMesh::framework_enabled(
        Amanzi::AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(Amanzi::AmanziMesh::Framework::MSTK);
    framework_names.push_back("MSTK");
  }

  for (int i = 0; i < frameworks.size(); i++) {
    // Set the framework
    std::cerr << "Testing deformation with " << framework_names[i] << std::endl;

    // Create the mesh
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::Preference prefs(meshfactory.preference());
      prefs.clear();
      prefs.push_back(frameworks[i]);

      meshfactory.set_preference(prefs);

      mesh = meshfactory.create(0.0, 0.0, 10.0, 10.0, 10, 10);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &aerr);

    CHECK_EQUAL(aerr, 0);


    // Deform the mesh

    Amanzi::AmanziMesh::Entity_ID_List nodeids;
    std::vector<Amanzi::AmanziGeometry::Point> finpos;
    Amanzi::AmanziGeometry::Point_List newpos;

    int status, nnodes;

    nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                                Amanzi::AmanziMesh::Parallel_type::OWNED);
    nodeids.resize(nnodes); 
    for (int j = 0; j < nnodes; j++) {
      nodeids[j] = j;

      Amanzi::AmanziGeometry::Point oldcoord(2), newcoord(2);

      mesh->node_get_coordinates(j, &oldcoord);

      newcoord.set(oldcoord[0], 0.5 * oldcoord[1]);

      newpos.push_back(newcoord);
    }

    status = mesh->deform(nodeids, newpos, false, finpos);


    CHECK_EQUAL(status, 1);


    // If the deformation was successful, the cell volumes should be half
    // of what they were

    int ncells = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                                    Amanzi::AmanziMesh::Parallel_type::ALL);

    for (int j = 0; j < ncells; j++) {
      double volume = mesh->cell_volume(j, false);
      CHECK_EQUAL(volume, 0.5);
    }

  } // for each framework i
}


TEST(MESH_DEFORM3D)
{
  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh
  std::vector<Amanzi::AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

  frameworks.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);
  framework_names.push_back("simple");

  if (Amanzi::AmanziMesh::framework_enabled(
        Amanzi::AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(Amanzi::AmanziMesh::Framework::MSTK);
    framework_names.push_back("MSTK");
  }

  for (int i = 0; i < frameworks.size(); i++) {
    // Set the framework
    std::cerr << "Testing deformation with " << framework_names[i] << std::endl;

    // Create the mesh
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::Preference prefs(meshfactory.preference());
      prefs.clear();
      prefs.push_back(frameworks[i]);

      meshfactory.set_preference(prefs);

      mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &aerr);

    CHECK_EQUAL(aerr, 0);


    // Deform the mesh

    Amanzi::AmanziMesh::Entity_ID_List nodeids;
    Amanzi::AmanziGeometry::Point_List newpos;
    std::vector<Amanzi::AmanziGeometry::Point> finpos;

    int status, nnodes;
    if (nproc == 1) {
      nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                                  Amanzi::AmanziMesh::Parallel_type::OWNED);
      nodeids.resize(nnodes);
      for (int j = 0; j < nnodes; j++) {
        double pi = 3.1415926;
        nodeids[j] = j;

        Amanzi::AmanziGeometry::Point oldcoord(3), newcoord(3);

        mesh->node_get_coordinates(j, &oldcoord);

        newcoord.set(oldcoord[0],
                     oldcoord[1],
                     oldcoord[2] - sin((oldcoord[0] + 1) * pi) *
                                     sin((oldcoord[1] + 1) * pi) * oldcoord[2] *
                                     0.2);

        newpos.push_back(newcoord);
      }

      status = mesh->deform(nodeids, newpos, false, finpos);

    } else {
      std::cerr << "Parallel deformation not implemented" << std::endl;
      status = 0;
    }

    CHECK_EQUAL(status, 1);


    // Check the deformations

    for (int j = 0; j < nnodes; j++) {
      Amanzi::AmanziGeometry::Point diff = finpos[j] - newpos[j];
      CHECK_EQUAL(diff[0], 0.0);
      CHECK_EQUAL(diff[1], 0.0);
      CHECK_EQUAL(diff[2], 0.0);
    }

  } // for each framework i
}
