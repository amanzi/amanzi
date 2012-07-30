/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   test_deform.cc
 * @author Rao V. Garimella
 * @date   Thu Apr 12, 2012
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include <Epetra_MpiComm.h>

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FrameworkTraits.hh"
#include "Geometry.hh"

TEST(MESH_DEFORM2D)
{

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh

  const Amanzi::AmanziMesh::Framework frameworks[] = {  
    Amanzi::AmanziMesh::MSTK
  };
  const char *framework_names[] = {
    "MSTK"
  };

  // const Amanzi::AmanziMesh::Framework frameworks[] = {  
  //   Amanzi::AmanziMesh::STKMESH, 
  //   Amanzi::AmanziMesh::MSTK, 
  //   Amanzi::AmanziMesh::Simple
  // };
  // const char *framework_names[] = {
  //   "stk::mesh", "MSTK", "Simple"
  // };

  const int numframeworks = sizeof(frameworks)/sizeof(Amanzi::AmanziMesh::Framework);


  Amanzi::AmanziMesh::Framework the_framework;
  for (int i = 0; i < numframeworks; i++) {


    // Set the framework

    the_framework = frameworks[i];

    if (!Amanzi::AmanziMesh::framework_available(the_framework)) continue;

    std::cerr << "Testing deformation with " << framework_names[i] << std::endl;


    // Create the mesh

    Amanzi::AmanziMesh::MeshFactory factory(&comm);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
      prefs.clear(); 
      prefs.push_back(the_framework);

      factory.preference(prefs);

      mesh = factory(0.0,0.0,1.0,1.0,10,10);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm.SumAll(&ierr, &aerr, 1);

    CHECK_EQUAL(aerr,0);


    // Deform the mesh

    Amanzi::AmanziMesh::Entity_ID_List nodeids;
    Amanzi::AmanziGeometry::Point_List newpos, finpos;

    int status, nnodes;
    if (nproc == 1) {

      nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                                 Amanzi::AmanziMesh::OWNED);

      for (int j = 0; j < nnodes; j++) {
        double pi = 3.1415926;
        nodeids.push_back(j);

        Amanzi::AmanziGeometry::Point oldcoord(2),newcoord(2);

        mesh->node_get_coordinates(j,&oldcoord);

        newcoord.set(oldcoord[0],oldcoord[1]+sin((oldcoord[0]+1)*pi)*oldcoord[1]*0.2);

        newpos.push_back(newcoord);
      }

      status = mesh->deform(nodeids,newpos,false,&finpos);
      
    } else {
      
      std::cerr << "Parallel deformation not implemented" << std::endl;
      status = 0;

    }

    CHECK_EQUAL(status,1);


    // Check the deformations

    for (int j = 0; j < nnodes; j++) {
      Amanzi::AmanziGeometry::Point diff = finpos[j]-newpos[j];
      CHECK_EQUAL(diff[0],0.0);
      CHECK_EQUAL(diff[1],0.0);
    }

  } // for each framework i


}



TEST(MESH_DEFORM3D)
{

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh

  const Amanzi::AmanziMesh::Framework frameworks[] = {  
    Amanzi::AmanziMesh::STKMESH,
    Amanzi::AmanziMesh::MSTK, 
    Amanzi::AmanziMesh::Simple
  };
  const char *framework_names[] = {
    "stk::mesh", "MSTK", "Simple"
  };
  
  const int numframeworks = sizeof(frameworks)/sizeof(Amanzi::AmanziMesh::Framework);


  Amanzi::AmanziMesh::Framework the_framework;
  for (int i = 0; i < numframeworks; i++) {


    // Set the framework

    the_framework = frameworks[i];

    if (!Amanzi::AmanziMesh::framework_available(the_framework)) continue;

    std::cerr << "Testing deformation with " << framework_names[i] << std::endl;


    // Create the mesh

    Amanzi::AmanziMesh::MeshFactory factory(&comm);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    int ierr = 0;
    int aerr = 0;
    try {
      Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
      prefs.clear(); 
      prefs.push_back(the_framework);

      factory.preference(prefs);

      mesh = factory(0.0,0.0,0.0,1.0,1.0,1.0,10,10,10);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm.SumAll(&ierr, &aerr, 1);

    CHECK_EQUAL(aerr,0);


    // Deform the mesh

    Amanzi::AmanziMesh::Entity_ID_List nodeids;
    Amanzi::AmanziGeometry::Point_List newpos, finpos;

    int status, nnodes;
    if (nproc == 1) {

      nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                                 Amanzi::AmanziMesh::OWNED);

      for (int j = 0; j < nnodes; j++) {
        double pi = 3.1415926;
        nodeids.push_back(j);

        Amanzi::AmanziGeometry::Point oldcoord(3),newcoord(3);

        mesh->node_get_coordinates(j,&oldcoord);

        newcoord.set(oldcoord[0],oldcoord[1],
                     oldcoord[2]-sin((oldcoord[0]+1)*pi)*sin((oldcoord[1]+1)*pi)*oldcoord[2]*0.2);

        newpos.push_back(newcoord);
      }

      status = mesh->deform(nodeids,newpos,false,&finpos);
      
    } else {
      
      std::cerr << "Parallel deformation not implemented" << std::endl;
      status = 0;

    }

    CHECK_EQUAL(status,1);


    // Check the deformations

    for (int j = 0; j < nnodes; j++) {
      Amanzi::AmanziGeometry::Point diff = finpos[j]-newpos[j];
      CHECK_EQUAL(diff[0],0.0);
      CHECK_EQUAL(diff[1],0.0);
      CHECK_EQUAL(diff[2],0.0);
    }

  } // for each framework i


}
