/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   test_Hex.cc
 * @author William A. Perkins
 * @date Fri Jul 29 09:37:58 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 18, 2010 by William A. Perkins
// Last Change: Fri Jul 29 09:37:58 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <UnitTest++.h>
#include <Epetra_MpiComm.h>

#include "../Mesh_STK_Impl.hh"
#include "../Mesh_STK_factory.hh"
#include "../Data_structures.hh"
#include "HexMeshGenerator.hh"
#include "Auditor.hh"


SUITE (HexMesh)
{
  TEST (HexMesh)
  {
    const unsigned int isize(3), jsize(3), ksize(3);

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    const int nproc(comm.NumProc());
    const int me(comm.MyPID());

    stk::ParallelMachine pm(comm.Comm());

    Amanzi::AmanziMesh::Data::HexMeshGenerator g(&comm, isize, jsize, ksize);

    Teuchos::RCP<Amanzi::AmanziMesh::Data::Data> meshdata(g.generate());

    for (int p = 0; p < nproc; p++) {
      if (me == p) {
        std::cerr << std::endl;
        std::cerr << ">>>>>> Process " << p << " Begin <<<<<<" << std::endl;
        meshdata->to_stream(std::cerr, true);
        std::cerr << ">>>>>> Process " << p << " End <<<<<<" << std::endl;
        std::cerr << std::endl;
      }
      comm.Barrier();
    }

    // need to have 1-based global indexes for stk::mesh
    Teuchos::RCP<Epetra_Map> cmap(g.cellmap(true));

    CHECK_EQUAL(cmap->NumGlobalElements(), isize*jsize*ksize);
    CHECK_EQUAL(cmap->MinAllGID(), 1);
    CHECK_EQUAL(cmap->MaxAllGID(), isize*jsize*ksize);

    // need to have 1-based global indexes for stk::mesh
    Teuchos::RCP<Epetra_Map> vmap(g.vertexmap(true));
    CHECK_EQUAL(vmap->MinAllGID(), 1);
    CHECK_EQUAL(vmap->MaxAllGID(), (isize+1)*(jsize+1)*(ksize+1));


    Amanzi::AmanziMesh::STK::Mesh_STK_factory mf(pm, 1000);
    Amanzi::AmanziMesh::Data::Fields nofields;
    Amanzi::AmanziMesh::STK::Mesh_STK_Impl_p 
        mesh(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));

    CHECK_EQUAL (mesh->rank_id (), me);
        
    int lcount, gcount;
    lcount = mesh->count_entities(stk::mesh::Element,Amanzi::AmanziMesh::OWNED);
    comm.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, isize*jsize*ksize);

    lcount = mesh->count_entities(stk::mesh::Face, Amanzi::AmanziMesh::OWNED);
    comm.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, 
                 (isize  )*(jsize  )*(ksize+1) + 
                 (isize  )*(jsize+1)*(ksize ) + 
                 (isize+1)*(jsize  )*(ksize ));

    lcount = mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED);
    comm.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, (isize+1)*(jsize+1)*(ksize+1));

    stk::mesh::Part *side;

    side = mesh->get_set("West", stk::mesh::Face);
    lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::OWNED);
    comm.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, isize*jsize);

    side = mesh->get_set("East", stk::mesh::Face);
    lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::OWNED);
    comm.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, isize*jsize);

    side = mesh->get_set("South", stk::mesh::Face);
    lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::OWNED);
    comm.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, isize*ksize);

    side = mesh->get_set("North", stk::mesh::Face);
    lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::OWNED);
    comm.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, isize*ksize);

    side = mesh->get_set("Bottom", stk::mesh::Face);
    lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::OWNED);
    comm.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, isize*jsize);

    side = mesh->get_set("East", stk::mesh::Face);
    lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::OWNED);
    comm.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, jsize*ksize);

    mesh->summary(std::cerr);

    Auditor audit("stk_mesh_hextest1_", mesh);
    audit();

  }

  TEST (HexGhosting)
  {
    const unsigned int isize(1), jsize(1), ksize(8);

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    const int nproc(comm.NumProc());
    const int me(comm.MyPID());

    stk::ParallelMachine pm(comm.Comm());

    Amanzi::AmanziMesh::Data::HexMeshGenerator g(&comm, isize, jsize, ksize);

    Teuchos::RCP<Amanzi::AmanziMesh::Data::Data> meshdata(g.generate());
    Teuchos::RCP<Epetra_Map> cmap(g.cellmap(true));
    Teuchos::RCP<Epetra_Map> vmap(g.vertexmap(true));

    Amanzi::AmanziMesh::STK::Mesh_STK_factory mf(pm, 1000);
    Amanzi::AmanziMesh::Data::Fields nofields;
    Amanzi::AmanziMesh::STK::Mesh_STK_Impl_p 
        mesh(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));

    Amanzi::AmanziMesh::STK::Entity_vector e;

    int ncell(mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED));

    if (nproc > 1) {

      mesh->summary(std::cerr);

      // all processes should have at least 1 but at most 8 shared nodes

      mesh->get_entities(stk::mesh::Node, Amanzi::AmanziMesh::GHOST, e);
      // CHECK(e.size() <= 8);
      e.clear();

      // processes > 1 should have only 1 ghost face

      mesh->get_entities(stk::mesh::Face, Amanzi::AmanziMesh::GHOST, e);

      if (me == 0) {
        // CHECK(e.empty());
      } else {
        // CHECK_EQUAL(e.size(), 1);
      }
      e.clear();

      // the number of USED faces depends on the number of cells owned

      // mesh->get_entities(stk::mesh::Face, Amanzi::AmanziMesh::USED, e);

      int nface_expected(ncell*5+1);

      // CHECK(!e.empty());
      // CHECK_EQUAL(e.size(), nface_expected);
      e.clear();
        
      // processes should have at least 1 but at most 2 shared
      // cells, but it doesn't
            
      mesh->get_entities(stk::mesh::Element, Amanzi::AmanziMesh::GHOST, e);
      e.clear();
           
      mesh->get_entities(stk::mesh::Element, Amanzi::AmanziMesh::USED, e);
      e.clear();

      // CHECK(!e.empty());
      // CHECK(e.size() <= 2);

      // for (int p = 0; p < nproc; p++) {
      //     if (me == p) {
      //         Amanzi::AmanziMesh::STK::Entity_vector nodes;
      //         mesh->get_entities(stk::mesh::Node, USED, nodes);
      //         for (unsigned int i = 0; i < nodes.size(); i++) {
      //             unsigned int gid(nodes[i]->identifier());
      //             const double *coord = mesh->coordinates(gid);
      //             std::cerr << "Node " << gid << ": "
      //                       << coord[0] << ", "
      //                       << coord[1] << ", "
      //                       << coord[2] << std::endl;
      //         }
      //         std::cerr << ">>>>>> Process " << p << " End <<<<<< " << std::endl;
      //     }
      //     comm.Barrier();
      // }

    } else {

      mesh->get_entities(stk::mesh::Node, Amanzi::AmanziMesh::GHOST, e);
      CHECK(e.empty());

      mesh->get_entities(stk::mesh::Face, Amanzi::AmanziMesh::GHOST, e);
      CHECK(e.empty());
        
      mesh->get_entities(stk::mesh::Element, Amanzi::AmanziMesh::GHOST, e);
      CHECK(e.empty());
    }            

    Auditor audit("stk_mesh_hextest2_", mesh);
    audit();
  } 

  TEST (HexGenerator)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
        mesh_map(new Amanzi::AmanziMesh::Mesh_STK(comm, 10, 10, 10));
     
    Auditor audit("stk_mesh_generated_", mesh_map);
    audit();
  }

  TEST (HexGeneratorParam)
  {
    Teuchos::ParameterList parameter_list;
    parameter_list.set<int>("Numer of Cells in X", 10);
    parameter_list.set<int>("Numer of Cells in Y", 10);
    parameter_list.set<int>("Numer of Cells in Z", 10);
    
    parameter_list.set<double>("X_Min", 0);
    parameter_list.set<double>("X_Max", 1);
    
    parameter_list.set<double>("Y_Min", 0);
    parameter_list.set<double>("Y_Max", 1);
    
    parameter_list.set<double>("Z_Min", 0);
    parameter_list.set<double>("Z_Max", 1);

    parameter_list.set<int>("Number of mesh blocks",0);

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
        mesh_map(new Amanzi::AmanziMesh::Mesh_STK(parameter_list, &comm));
     
    Auditor audit("stk_mesh_generated_", mesh_map);
    audit();
  }

  TEST (HexPartition)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::AmanziMesh::Mesh_STK
        *mesh_stk = new Amanzi::AmanziMesh::Mesh_STK(comm, 4, 2, 2);

    Teuchos::RCP<Epetra_CrsGraph> cgraph = mesh_stk->cellgraph();
    cgraph->Print(std::cerr);

    mesh_stk->redistribute();

    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(mesh_stk);
    Auditor audit("stk_mesh_rpartitioned_", mesh);
    audit();
  }
}

