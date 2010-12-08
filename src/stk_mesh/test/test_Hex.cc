// -------------------------------------------------------------
/**
 * @file   test_Hex.cc
 * @author William A. Perkins
 * @date Wed Dec  8 08:01:10 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 18, 2010 by William A. Perkins
// Last Change: Wed Dec  8 08:01:10 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <UnitTest++.h>
#include <Epetra_MpiComm.h>

#include "../Mesh.hh"
#include "../Mesh_factory.hh"
#include "../Data_structures.hh"
#include "../Mesh_maps_stk.hh"
#include "MeshAudit.hh"
#include "HexMeshGenerator.hh"


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
        lcount = mesh->count_entities(stk::mesh::Element, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*jsize*ksize);

        lcount = mesh->count_entities(stk::mesh::Face, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, 
                     (isize  )*(jsize  )*(ksize+1) + 
                     (isize  )*(jsize+1)*(ksize ) + 
                     (isize+1)*(jsize  )*(ksize ));

        lcount = mesh->count_entities(stk::mesh::Node, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, (isize+1)*(jsize+1)*(ksize+1));

        stk::mesh::Part *side;

        side = mesh->get_set("West", stk::mesh::Face);
        lcount = mesh->count_entities(*side, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*jsize);

        side = mesh->get_set("East", stk::mesh::Face);
        lcount = mesh->count_entities(*side, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*jsize);

        side = mesh->get_set("South", stk::mesh::Face);
        lcount = mesh->count_entities(*side, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*ksize);

        side = mesh->get_set("North", stk::mesh::Face);
        lcount = mesh->count_entities(*side, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*ksize);

        side = mesh->get_set("Bottom", stk::mesh::Face);
        lcount = mesh->count_entities(*side, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, isize*jsize);

        side = mesh->get_set("East", stk::mesh::Face);
        lcount = mesh->count_entities(*side, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, jsize*ksize);

        mesh->summary(std::cerr);

        // Teuchos::RCP<Mesh_maps_base> mesh_map(new STK_mesh::Mesh_maps_stk(mesh));

        // mesh_map->cell_map(false).Print(std::cout);

        // if (comm.NumProc() == 1) {
        //   MeshAudit audit(mesh_map);
        //   CHECK(audit.check_entity_counts());
        // } else {
        //   std::ostringstream ofile;
        //   ofile << "stk_mesh_audit_" << std::setfill('0') << std::setw(4) << comm.MyPID() << ".txt";
        //   std::ofstream ofs(ofile.str().c_str());
        //   if (comm.MyPID() == 0)
        //     std::cout << "Writing results to " << ofile.str() << ", etc." << std::endl;
        //   MeshAudit audit(mesh_map, ofs);
        //   CHECK(audit.check_entity_counts());
        // }
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

        int ncell(mesh->count_entities(stk::mesh::Element, STK_mesh::OWNED));

        if (nproc > 1) {

            mesh->summary(std::cerr);

            // all processes should have at least 1 but at most 8 shared nodes

            mesh->get_entities(stk::mesh::Node, STK_mesh::GHOST, e);
            CHECK(e.size() <= 8);
            e.clear();

            // processes > 1 should have only 1 ghost face

            mesh->get_entities(stk::mesh::Face, STK_mesh::GHOST, e);

            if (me == 0) {
                CHECK(e.empty());
            } else {
                CHECK_EQUAL(e.size(), 1);
            }
            e.clear();

            // the number of USED faces depends on the number of cells owned

            mesh->get_entities(stk::mesh::Face, STK_mesh::USED, e);

            int nface_expected(ncell*5+1);

            CHECK(!e.empty());
            CHECK_EQUAL(e.size(), nface_expected);
            e.clear();
        
            // processes should have at least 1 but at most 2 shared
            // cells, but it doesn't
            
            mesh->get_entities(stk::mesh::Element, STK_mesh::GHOST, e);
            e.clear();
           
            mesh->get_entities(stk::mesh::Element, STK_mesh::USED, e);
            e.clear();

            // CHECK(!e.empty());
            // CHECK(e.size() <= 2);

        } else {

            mesh->get_entities(stk::mesh::Node, STK_mesh::GHOST, e);
            CHECK(e.empty());

            mesh->get_entities(stk::mesh::Face, STK_mesh::GHOST, e);
            CHECK(e.empty());
        
            mesh->get_entities(stk::mesh::Element, STK_mesh::GHOST, e);
            CHECK(e.empty());
        }            
    }        
}

