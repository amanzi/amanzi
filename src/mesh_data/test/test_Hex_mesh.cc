// -------------------------------------------------------------
/**
 * @file   test_Hex_mesh.cc
 * @author William A. Perkins
 * @date Thu Apr 21 13:27:01 2011
 * 
 * @brief  A set of tests for the HexMeshGenerator class
 * 
 * 
 */
// -------------------------------------------------------------


#include <UnitTest++.h>
#include <iostream>
#include <stdexcept>
#include <boost/format.hpp>
#include "Epetra_MpiComm.h"
#include "../HexMeshGenerator.hh"
#include "errors.hh"

const static unsigned int size(4);

SUITE (HexMeshGenerator)
{
  TEST (Generation)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    int me(comm.MyPID());

    Mesh_data::HexMeshGenerator gen(&comm, size*size, size, size);
    Mesh_data::Data *mesh;
    mesh = gen.generate();

    mesh->to_stream(std::cout, true);

    // FIXME: do some checks on mesh

    try {
      std::auto_ptr<Epetra_Map> cmap(gen.cellmap(true));
      CHECK_EQUAL(cmap->NumGlobalElements(), (size*size)*size*size);
      CHECK_EQUAL(cmap->NumGlobalElements(), gen.cells());
      CHECK_EQUAL(cmap->NumMyElements(), gen.mycells());
      cmap->Print(std::cerr);   // ends up in the test log file
    } catch (int e) {
      
      // it appears that Epetra_Map's throw an int when in trouble,
      // let's make into an Amanzi exception

      std::string msg = 
        boost::str(boost::format("Cell Epetra_Map error: %d") % e);
      Exceptions::amanzi_throw(Errors::Message(msg.c_str()));
    }

    try {
      std::auto_ptr<Epetra_Map> cmap(gen.vertexmap(true));
      CHECK_EQUAL(cmap->MaxAllGID(), (size*size+1)*(size+1)*(size+1));
      CHECK_EQUAL(cmap->MinAllGID(), 1);
      CHECK_EQUAL(cmap->NumMyElements(), gen.myvertexes());
      cmap->Print(std::cerr);   // ends up in the test log file
    } catch (int e) {
      std::string msg = 
        boost::str(boost::format("Vertex Epetra_Map error: %d") % e);
      Exceptions::amanzi_throw(Errors::Message(msg.c_str()));
    }

    comm.Barrier();             // probably not necessary

    // if it runs, it passes the test, right?
  }
}
