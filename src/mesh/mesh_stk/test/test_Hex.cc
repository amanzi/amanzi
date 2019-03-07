/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   test_Hex.cc
 * @author William A. Perkins
 * @date Wed Sep 28 09:03:25 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 18, 2010 by William A. Perkins
// Last Change: Wed Sep 28 09:03:25 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <UnitTest++.h>
#include <AmanziComm.hh>

#include "../Mesh_STK_Impl.hh"
#include "../Mesh_STK_factory.hh"
#include "../Data_structures.hh"
#include "HexMeshGenerator.hh"
#include "GenerationSpec.hh"
#include "Auditor.hh"

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"


SUITE (HexMesh)
{
  TEST (HexMesh)
  {
    const unsigned int isize(3), jsize(3), ksize(3);

    Epetra_MpiComm comm_(MPI_COMM_WORLD);
    const int nproc(comm_.NumProc());
    const int me(comm_.MyPID());

    stk::ParallelMachine pm(comm_.Comm());

    Amanzi::AmanziMesh::Data::HexMeshGenerator g(&comm_, isize, jsize, ksize);

    Teuchos::RCP<Amanzi::AmanziMesh::Data::Data> meshdata(g.generate());

    // Disable printing - must check/verify, not very useful to just print
    // for (int p = 0; p < nproc; p++) {
    //   if (me == p) {
    //     std::cerr << std::endl;
    //     std::cerr << ">>>>>> Process " << p << " Begin <<<<<<" << std::endl;
    //     meshdata->to_stream(std::cerr, true);
    //     std::cerr << ">>>>>> Process " << p << " End <<<<<<" << std::endl;
    //     std::cerr << std::endl;
    //   }
    //   comm_.Barrier();
    // }

    // need to have 1-based global indexes for stk::mesh
    Teuchos::RCP<Epetra_Map> cmap(g.cellmap(true));

    CHECK_EQUAL(cmap->NumGlobalElements(), isize*jsize*ksize);
    CHECK_EQUAL(cmap->MinAllGID(), 1);
    CHECK_EQUAL(cmap->MaxAllGID(), isize*jsize*ksize);

    // need to have 1-based global indexes for stk::mesh
    Teuchos::RCP<Epetra_Map> vmap(g.vertexmap(true));
    CHECK_EQUAL(vmap->MinAllGID(), 1);
    CHECK_EQUAL(vmap->MaxAllGID(), (isize+1)*(jsize+1)*(ksize+1));


    Amanzi::AmanziMesh::Framework::STK::Mesh_STK_factory mf(&comm_, 1000);
    Amanzi::AmanziMesh::Data::Fields nofields;
    Amanzi::AmanziMesh::Framework::STK::Mesh_STK_Impl_p 
      mesh(mf.build_mesh(*meshdata, *cmap, *vmap, nofields, Teuchos::null));

    CHECK_EQUAL (mesh->rank_id (), me);
        
    int lcount, gcount;
    lcount = mesh->count_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::CELL),Amanzi::AmanziMesh::Parallel_type::OWNED);
    comm_.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, isize*jsize*ksize);

    lcount = mesh->count_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::FACE), Amanzi::AmanziMesh::Parallel_type::OWNED);
    comm_.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, 
                 (isize  )*(jsize  )*(ksize+1) + 
                 (isize  )*(jsize+1)*(ksize ) + 
                 (isize+1)*(jsize  )*(ksize ));

    lcount = mesh->count_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::NODE), Amanzi::AmanziMesh::Parallel_type::OWNED);
    comm_.SumAll(&lcount, &gcount, 1);
    CHECK_EQUAL (gcount, (isize+1)*(jsize+1)*(ksize+1));



    // Check sets in a different test
    // stk::mesh::Part *side;

    // side = mesh->get_set("West", mesh->kind_to_rank(Amanzi::AmanziMesh::FACE));
    // lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::Parallel_type::OWNED);
    // comm_.SumAll(&lcount, &gcount, 1);
    // CHECK_EQUAL (gcount, isize*jsize);

    // side = mesh->get_set("East", mesh->kind_to_rank(Amanzi::AmanziMesh::FACE));
    // lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::Parallel_type::OWNED);
    // comm_.SumAll(&lcount, &gcount, 1);
    // CHECK_EQUAL (gcount, isize*jsize);

    // side = mesh->get_set("South", mesh->kind_to_rank(Amanzi::AmanziMesh::FACE));
    // lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::Parallel_type::OWNED);
    // comm_.SumAll(&lcount, &gcount, 1);
    // CHECK_EQUAL (gcount, isize*ksize);

    // side = mesh->get_set("North", mesh->kind_to_rank(Amanzi::AmanziMesh::FACE));
    // lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::Parallel_type::OWNED);
    // comm_.SumAll(&lcount, &gcount, 1);
    // CHECK_EQUAL (gcount, isize*ksize);

    // side = mesh->get_set("Bottom", mesh->kind_to_rank(Amanzi::AmanziMesh::FACE));
    // lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::Parallel_type::OWNED);
    // comm_.SumAll(&lcount, &gcount, 1);
    // CHECK_EQUAL (gcount, isize*jsize);

    // side = mesh->get_set("East", mesh->kind_to_rank(Amanzi::AmanziMesh::FACE));
    // lcount = mesh->count_entities(*side, Amanzi::AmanziMesh::Parallel_type::OWNED);
    // comm_.SumAll(&lcount, &gcount, 1);
    // CHECK_EQUAL (gcount, jsize*ksize);


    // Must check/verify, not just print
    //  mesh->summary(std::cerr);

    // Must disable until STKmesh can give us contiguous entity IDs

    // Auditor audit("stk_mesh_hextest1_", mesh);
    //audit();

  }

  TEST (HexGhosting)
  {
    const unsigned int isize(1), jsize(1), ksize(8);

    Epetra_MpiComm comm_(MPI_COMM_WORLD);
    const int nproc(comm_.NumProc());
    const int me(comm_.MyPID());

    stk::ParallelMachine pm(comm_.Comm());

    Amanzi::AmanziMesh::Data::HexMeshGenerator g(&comm_, isize, jsize, ksize);

    Teuchos::RCP<Amanzi::AmanziMesh::Data::Data> meshdata(g.generate());
    Teuchos::RCP<Epetra_Map> cmap(g.cellmap(true));
    Teuchos::RCP<Epetra_Map> vmap(g.vertexmap(true));

    Amanzi::AmanziMesh::Framework::STK::Mesh_STK_factory mf(&comm_, 1000);
    Amanzi::AmanziMesh::Data::Fields nofields;
    Amanzi::AmanziMesh::Framework::STK::Mesh_STK_Impl_p 
      mesh(mf.build_mesh(*meshdata, *cmap, *vmap, nofields, Teuchos::null));

    Amanzi::AmanziMesh::Framework::STK::Entity_vector e;

    int ncell(mesh->count_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::CELL), Amanzi::AmanziMesh::Parallel_type::OWNED));

    if (nproc > 1) {

      //      mesh->summary(std::cerr);

      // all processes should have at least 1 but at most 8 shared nodes

      mesh->get_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::NODE), Amanzi::AmanziMesh::Parallel_type::GHOST, e);
      // CHECK(e.size() <= 8);
      e.clear();

      // processes > 1 should have only 1 ghost face

      mesh->get_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::FACE), Amanzi::AmanziMesh::Parallel_type::GHOST, e);

      if (me == 0) {
        // CHECK(e.empty());
      } else {
        // CHECK_EQUAL(e.size(), 1);
      }
      e.clear();

      // the number of ALL faces depends on the number of cells owned

      // mesh->get_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::FACE), Amanzi::AmanziMesh::Parallel_type::ALL, e);

      int nface_expected(ncell*5+1);

      // CHECK(!e.empty());
      // CHECK_EQUAL(e.size(), nface_expected);
      e.clear();
        
      // processes should have at least 1 but at most 2 shared
      // cells, but it doesn't
            
      mesh->get_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::CELL), Amanzi::AmanziMesh::Parallel_type::GHOST, e);
      e.clear();
           
      mesh->get_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::CELL), Amanzi::AmanziMesh::Parallel_type::ALL, e);
      e.clear();

      // CHECK(!e.empty());
      // CHECK(e.size() <= 2);

      // for (int p = 0; p < nproc; p++) {
      //     if (me == p) {
      //         Amanzi::AmanziMesh::Framework::STK::Entity_vector nodes;
      //         mesh->get_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::NODE), Parallel_type::ALL, nodes);
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
      //     comm_.Barrier();
      // }

    } else {

      mesh->get_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::NODE), Amanzi::AmanziMesh::Parallel_type::GHOST, e);
      CHECK(e.empty());

      mesh->get_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::FACE), Amanzi::AmanziMesh::Parallel_type::GHOST, e);
      CHECK(e.empty());
        
      mesh->get_entities(mesh->kind_to_rank(Amanzi::AmanziMesh::CELL), Amanzi::AmanziMesh::Parallel_type::GHOST, e);
      CHECK(e.empty());
    }            

    // Must disable until STKmesh can give us contiguous entity IDs

    //    Auditor audit("stk_mesh_hextest2_", mesh);
    //    audit();
  } 

  TEST (HexGenerator)
  {
    Epetra_MpiComm comm_(MPI_COMM_WORLD);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
        mesh_map(new Amanzi::AmanziMesh::Mesh_STK(&comm_, 10, 10, 10));
    
    // Must disable until STKmesh can give us contiguous entity IDs
 
    // Auditor audit("stk_mesh_generated_", mesh_map);
    // audit();
  }

  TEST (HexGeneratorParam)
  {

    // FIXME: I DON'T KNOW HOW TO INITIALIZE ARRAYS USING TEUCHOS PARAMETER LISTS
    //    Teuchos::ParameterList parameter_list;
    //    parameter_list.set< Teuchos::Array<double> >("domain low coordinate",{0.0, 0.0, 0.0});
    //    parameter_list.set< Teuchos::Array<double> >("domain high coordinate",{1.0, 1.0, 1.0});
    //    parameter_list.set< Teuchos::Array<int> >("number of cells", {10, 10, 10});

    //    Amanzi::AmanziMesh::GenerationSpec gspec(parameter_list);

    //    Epetra_MpiComm comm_(MPI_COMM_WORLD);
    //    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
    //        mesh_map(new Amanzi::AmanziMesh::Mesh_STK(gspec, &comm_));
     
    //    Auditor audit("stk_mesh_generated_", mesh_map);
    //    audit();

  }

  TEST (HexPartition)
  {
    Epetra_MpiComm comm_(MPI_COMM_WORLD);
    Amanzi::AmanziMesh::Mesh_STK
        *mesh_stk = new Amanzi::AmanziMesh::Mesh_STK(&comm_, 4, 2, 2);

    Teuchos::RCP<Epetra_CrsGraph> cgraph = mesh_stk->cellgraph();
    //    cgraph->Print(std::cerr);  // must check, not just print!
    

    mesh_stk->redistribute();

    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(mesh_stk);

    // Disable until STKmesh can give us contiguous IDs

    //    Auditor audit("stk_mesh_rpartitioned_", mesh);
    //    audit();
  }
}

