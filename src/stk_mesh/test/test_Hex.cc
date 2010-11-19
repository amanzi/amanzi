// -------------------------------------------------------------
/**
 * @file   test_Hex.cc
 * @author William A. Perkins
 * @date Fri Nov 19 07:36:30 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 18, 2010 by William A. Perkins
// Last Change: Fri Nov 19 07:36:30 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <UnitTest++.h>
#include <Epetra_MpiComm.h>

#include "../Mesh.hh"
#include "../Mesh_factory.hh"
#include "../Data_structures.hh"
#include "HexMeshGenerator.hh"


const unsigned int size(3);

SUITE (HexMesh)
{
    TEST (HexMesh)
    {
        Epetra_MpiComm comm(MPI_COMM_WORLD);
        const int nproc(comm.NumProc());
        const int me(comm.MyPID());

        stk::ParallelMachine pm(comm.Comm());

        Mesh_data::HexMeshGenerator g(&comm, size, size, size);

        Teuchos::RCP<Mesh_data::Data> meshdata(g.generate());

        // need to have 1-based global indexes for stk::mesh
        Teuchos::RCP<Epetra_Map> cmap(g.cellmap(true));

        CHECK_EQUAL(cmap->NumGlobalElements(), size*size*size);
        CHECK_EQUAL(cmap->MinAllGID(), 1);
        CHECK_EQUAL(cmap->MaxAllGID(), size*size*size);

        // need to have 1-based global indexes for stk::mesh
        Teuchos::RCP<Epetra_Map> vmap(g.vertexmap(true));
        CHECK_EQUAL(vmap->MinAllGID(), 1);
        CHECK_EQUAL(vmap->MaxAllGID(), (size+1)*(size+1)*(size+1));


        STK_mesh::Mesh_factory mf(pm, 1000);
        Mesh_data::Fields nofields;
        STK_mesh::Mesh_p mesh(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));

        CHECK_EQUAL (mesh->rank_id (), me);
        
        int lcount, gcount;
        lcount = mesh->count_entities(stk::mesh::Element, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, size*size*size);
        
        lcount = mesh->count_entities(stk::mesh::Node, STK_mesh::OWNED);
        comm.SumAll(&lcount, &gcount, 1);
        CHECK_EQUAL (gcount, (size+1)*(size+1)*(size+1));
    }

}

