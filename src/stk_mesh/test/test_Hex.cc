// -------------------------------------------------------------
/**
 * @file   test_Hex.cc
 * @author William A. Perkins
 * @date Mon Dec 13 09:42:30 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 18, 2010 by William A. Perkins
// Last Change: Mon Dec 13 09:42:30 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <UnitTest++.h>
#include <Epetra_MpiComm.h>

#include "../Mesh.hh"
#include "../Mesh_factory.hh"
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

        Mesh_data::HexMeshGenerator g(&comm, isize, jsize, ksize);

        Teuchos::RCP<Mesh_data::Data> meshdata(g.generate());

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


        STK_mesh::Mesh_factory mf(pm, 1000);
        Mesh_data::Fields nofields;
        STK_mesh::Mesh_p mesh(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));

        CHECK_EQUAL (mesh->rank_id (), me);
        
        int lcount, gcount;
        lcount = mesh->count_entities(stk::mesh::Element,OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*jsize*ksize);

        lcount = mesh->count_entities(stk::mesh::Face, OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, 
                     (isize  )*(jsize  )*(ksize+1) + 
                     (isize  )*(jsize+1)*(ksize ) + 
                     (isize+1)*(jsize  )*(ksize ));

        lcount = mesh->count_entities(stk::mesh::Node, OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, (isize+1)*(jsize+1)*(ksize+1));

        stk::mesh::Part *side;

        side = mesh->get_set("West", stk::mesh::Face);
        lcount = mesh->count_entities(*side, OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*jsize);

        side = mesh->get_set("East", stk::mesh::Face);
        lcount = mesh->count_entities(*side, OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*jsize);

        side = mesh->get_set("South", stk::mesh::Face);
        lcount = mesh->count_entities(*side, OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*ksize);

        side = mesh->get_set("North", stk::mesh::Face);
        lcount = mesh->count_entities(*side, OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*ksize);

        side = mesh->get_set("Bottom", stk::mesh::Face);
        lcount = mesh->count_entities(*side, OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*jsize);

        side = mesh->get_set("East", stk::mesh::Face);
        lcount = mesh->count_entities(*side, OWNED);
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

        Mesh_data::HexMeshGenerator g(&comm, isize, jsize, ksize);

        Teuchos::RCP<Mesh_data::Data> meshdata(g.generate());
        Teuchos::RCP<Epetra_Map> cmap(g.cellmap(true));
        Teuchos::RCP<Epetra_Map> vmap(g.vertexmap(true));

        STK_mesh::Mesh_factory mf(pm, 1000);
        Mesh_data::Fields nofields;
        STK_mesh::Mesh_p mesh(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));

        STK_mesh::Entity_vector e;

        int ncell(mesh->count_entities(stk::mesh::Element, OWNED));

        if (nproc > 1) {

            mesh->summary(std::cerr);

            // all processes should have at least 1 but at most 8 shared nodes

            mesh->get_entities(stk::mesh::Node, GHOST, e);
            // CHECK(e.size() <= 8);
            e.clear();

            // processes > 1 should have only 1 ghost face

            mesh->get_entities(stk::mesh::Face, GHOST, e);

            if (me == 0) {
                // CHECK(e.empty());
            } else {
                // CHECK_EQUAL(e.size(), 1);
            }
            e.clear();

            // the number of USED faces depends on the number of cells owned

            // mesh->get_entities(stk::mesh::Face, USED, e);

            int nface_expected(ncell*5+1);

            // CHECK(!e.empty());
            // CHECK_EQUAL(e.size(), nface_expected);
            e.clear();
        
            // processes should have at least 1 but at most 2 shared
            // cells, but it doesn't
            
            mesh->get_entities(stk::mesh::Element, GHOST, e);
            e.clear();
           
            mesh->get_entities(stk::mesh::Element, USED, e);
            e.clear();

            // CHECK(!e.empty());
            // CHECK(e.size() <= 2);

            // for (int p = 0; p < nproc; p++) {
            //     if (me == p) {
            //         STK_mesh::Entity_vector nodes;
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

            mesh->get_entities(stk::mesh::Node, GHOST, e);
            CHECK(e.empty());

            mesh->get_entities(stk::mesh::Face, GHOST, e);
            CHECK(e.empty());
        
            mesh->get_entities(stk::mesh::Element, GHOST, e);
            CHECK(e.empty());
        }            

        Auditor audit("stk_mesh_hextest2_", mesh);
        audit();
    } 

    TEST (HexGenerator)
    {
        Epetra_MpiComm comm(MPI_COMM_WORLD);
        Teuchos::RCP<Mesh_maps_base> 
            mesh_map(new STK_mesh::Mesh_maps_stk(comm, 10, 10, 10));
        
        Auditor audit("stk_mesh_generated_", mesh_map);
        audit();
    }
}

