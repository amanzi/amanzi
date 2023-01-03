/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  State

*/

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "UnitTest++.h"

#include "IO.hh"
#include "Checkpoint.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "State.hh"

SUITE(RESTART2)
{
  /* This tests the ability to read checkpoint capabilities when the underlying
  meshes live on subcommunicators (in this case COMM_SELF and not on
  COMM_WORLD.
*/
  TEST(HDF5_MPI_AND_SERIAL)
  {
    using namespace Amanzi;

    auto comm = getDefaultComm();
    auto comm_serial = getCommSelf();

    std::string hdf5_datafile1 = "new_data_mpi_serial_parallel";

    int size = comm->NumProc();
    int rank = comm->MyPID();
    auto domain_mesh_mstk =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4 * size, 1, 1, comm));
    auto serial_mesh_mstk =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(0, 0, 0, 1, 1, 1, 4, 1, 1, comm_serial));
    auto domain_mesh = Teuchos::rcp(new AmanziMesh::Mesh(domain_mesh_mstk,Teuchos::null)); 
    auto serial_mesh = Teuchos::rcp(new AmanziMesh::Mesh(serial_mesh_mstk,Teuchos::null)); 


    auto S = Teuchos::rcp(new State());
    S->RegisterDomainMesh(domain_mesh);
    std::stringstream serial_dname_ss;
    serial_dname_ss << "serial_" << rank;

    std::string serial_dname = serial_dname_ss.str();
    Key serial_fname = Keys::getKey(serial_dname, "field");
    S->RegisterMesh(serial_dname, serial_mesh);

    S->Require<CompositeVector, CompositeVectorSpace>("my_field", Tags::DEFAULT, "my_field")
      .SetMesh(domain_mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->Require<CompositeVector, CompositeVectorSpace>(serial_fname, Tags::DEFAULT, serial_fname)
      .SetMesh(serial_mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    S->Setup();
    S->set_cycle(0);
    S->set_time(0.0);

    // Set values
    int cell_index_list[] = { 0, 1, 2, 3 };
    double cell_values[] = { 10.0, 20.0, 30.0, 40.0 };
    for (int i = 0; i != 4; ++i) { cell_values[i] += rank * 40; }
    auto& mf =
      *S->GetW<CompositeVector>("my_field", Tags::DEFAULT, "my_field").ViewComponent("cell");
    mf(0)->ReplaceMyValues(4, cell_values, cell_index_list);

    auto& sf =
      *S->GetW<CompositeVector>(serial_fname, Tags::DEFAULT, serial_fname).ViewComponent("cell");
    sf(0)->ReplaceMyValues(4, cell_values, cell_index_list);

    // make sure we can write them all... for real!
    Teuchos::ParameterList chkplist("checkpoint");
    chkplist.set("single file checkpoint", false);
    Checkpoint chkp(chkplist, *S);
    chkp.Write(*S);

    // make sure we can read them all
    State S2;
    S2.RegisterDomainMesh(domain_mesh);
    S2.RegisterMesh(serial_dname, serial_mesh);

    S2.Require<CompositeVector, CompositeVectorSpace>("my_field", Tags::DEFAULT, "my_field")
      .SetMesh(domain_mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S2.Require<CompositeVector, CompositeVectorSpace>(serial_fname, Tags::DEFAULT, serial_fname)
      .SetMesh(domain_mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S2.Setup();

    ReadCheckpoint(comm, S2, "checkpoint00000");

    S2.GetW<CompositeVector>("my_field", Tags::DEFAULT, "my_field")
      .Update(-1.0, S->Get<CompositeVector>("my_field"), 1.0);
    double diff;
    S2.Get<CompositeVector>("my_field").Norm2(&diff);
    CHECK_CLOSE(0.0, diff, 1.0e-10);

    S2.GetW<CompositeVector>(serial_fname, Tags::DEFAULT, serial_fname)
      .Update(-1.0, S->Get<CompositeVector>(serial_fname), 1.0);
    diff = 0.0;
    S2.Get<CompositeVector>(serial_fname).Norm2(&diff);
    CHECK_CLOSE(0.0, diff, 1.0e-10);
  }
}
