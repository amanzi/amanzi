// -------------------------------------------------------------
/**
 * @file   test_Hex_mesh.cc
 * @author William A. Perkins
 * @date Wed Nov 24 08:52:28 2010
 * 
 * @brief  A set of tests for the HexMeshGenerator class
 * 
 * 
 */
// -------------------------------------------------------------


#include <UnitTest++.h>
#include <iostream>
#include <stdexcept>
#include "Epetra_MpiComm.h"
#include "../HexMeshGenerator.hh"

const static unsigned int size(2);

SUITE (HexMeshGenerator)
{
  TEST (Generation)
  {
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    int me(comm.MyPID());

    Mesh_data::HexMeshGenerator gen(&comm, size, size, size);
    Mesh_data::Data *mesh;
    mesh = gen.generate();

    mesh->to_stream(std::cout, true);

    // FIXME: do some checks on mesh

    try {
      std::auto_ptr<Epetra_Map> cmap(gen.cellmap(true));
      CHECK_EQUAL(cmap->NumGlobalElements(), size*size*size);
      CHECK_EQUAL(cmap->NumGlobalElements(), gen.cells());
      CHECK_EQUAL(cmap->NumMyElements(), gen.mycells());
      cmap->Print(std::cerr);   // ends up in the test log file
    } catch (int e) {
      std::cerr << "Cell Epetra_Map error: " << e << std::endl;
      throw e;
    }

    try {
      std::auto_ptr<Epetra_Map> cmap(gen.vertexmap(true));
      CHECK_EQUAL(cmap->MaxAllGID(), (size+1)*(size+1)*(size+1));
      CHECK_EQUAL(cmap->MinAllGID(), 1);
      CHECK_EQUAL(cmap->NumMyElements(), gen.myvertexes());
      cmap->Print(std::cerr);   // ends up in the test log file
    } catch (int e) {
      std::cerr << "Vertex Epetra_Map error: " << e << std::endl;
      throw e;
    }

    comm.Barrier();             // probably not necessary

    // if it runs, it passes the test, right?
  }
}
