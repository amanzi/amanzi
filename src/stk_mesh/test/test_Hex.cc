// -------------------------------------------------------------
/**
 * @file   test_Hex.cc
 * @author William A. Perkins
 * @date Mon Nov 29 12:16:14 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 18, 2010 by William A. Perkins
// Last Change: Mon Nov 29 12:16:14 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <UnitTest++.h>
#include <Epetra_MpiComm.h>

#include "../Mesh.hh"
#include "../Mesh_factory.hh"
#include "../Data_structures.hh"
#include "../Mesh_maps_stk.hh"
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

        // STK_mesh::Mesh_maps_stk mesh_map(mesh);

        // lcount = mesh_map.count_entities(Mesh_data::CELL, STK_mesh::OWNED);
        // comm.SumAll(&lcount, &gcount, 1);
        // CHECK_EQUAL (gcount, isize*jsize*ksize);

        // lcount = mesh_map.count_entities(Mesh_data::FACE, STK_mesh::OWNED);
        // comm.SumAll(&lcount, &gcount, 1);
        // CHECK_EQUAL (gcount, 
        //              (isize  )*(jsize  )*(ksize+1) + 
        //              (isize  )*(jsize+1)*(ksize ) + 
        //              (isize+1)*(jsize  )*(ksize ));
    }

}

