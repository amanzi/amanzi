// -------------------------------------------------------------
/**
 * @file   test_Read.cc
 * @author William A. Perkins
 * @date Tue Nov 23 06:44:36 2010
 * 
 * @brief Some unit tests for reading a (serial) Exodus file and
 * building a STK_mesh::Mesh_maps_stk instance.
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 22, 2010 by William A. Perkins
// Last Change: Tue Nov 23 06:44:36 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <UnitTest++.h>
#include <Epetra_MpiComm.h>

#include "Exodus_Readers.hh"
#include "Parallel_Exodus_file.hh"
#include "../Mesh_maps_stk.hh"
#include "../Mesh_factory.hh"

SUITE (Exodus)
{
    TEST (Reader)
    {
        Epetra_MpiComm comm(MPI_COMM_WORLD);
        const int nproc(comm.NumProc());
        Teuchos::RCP<Mesh_data::Data> meshdata;
        STK_mesh::Mesh_factory mf(comm.Comm(), 1000);
        Mesh_data::Fields nofields;
        STK_mesh::Mesh_p mesh;

        if (nproc == 1) {
            meshdata.reset(ExodusII::read_exodus_file("../exodus/test_files/hex_2x2x2_ss.exo"));
            mesh.reset(mf.build_mesh(*meshdata, nofields));
        } else if (nproc <= 4) {
            ExodusII::Parallel_Exodus_file fileset(comm, "../exodus/test_files/split1/hex_11x11x11_ss.par");
            meshdata = fileset.read_mesh();
            Teuchos::RCP<Epetra_Map> cmap, vmap;
            mesh.reset(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));
        }

        // STK_mesh::Mesh_maps_stk maps(mesh);
          
    }
}
