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


TEST(MESH_COLUMNS)
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
  
  const int numframeworks = sizeof(frameworks)/sizeof(Amanzi::AmanziMesh::Framework);


  Amanzi::AmanziMesh::Framework the_framework;
  for (int i = 0; i < numframeworks; i++) {

    // Set the framework

    the_framework = frameworks[i];

    if (!Amanzi::AmanziMesh::framework_available(the_framework)) continue;

    std::cerr << "Testing columns with " << framework_names[i] << std::endl;


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

      mesh = factory(0.0,0.0,0.0,1.0,1.0,1.0,4,4,4);

    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm.SumAll(&ierr, &aerr, 1);

    CHECK_EQUAL(aerr,0);


    int status, nnodes;
    int ncells = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                                    Amanzi::AmanziMesh::OWNED);

    for (int j = 0; j < ncells; j++) {
      int expcellabove = ((j+1)%4 == 0) ? -1 : j+1;
      CHECK_EQUAL(expcellabove,mesh->cell_get_cell_above(j));
      int expcellbelow = (j%4 == 0) ? -1 : j-1;
      CHECK_EQUAL(expcellbelow,mesh->cell_get_cell_below(j));
    }
      
  } // for each framework i

}
