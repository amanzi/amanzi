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

#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "UnitTest++.h"

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

    int size = comm->getSize();
    int rank = comm->getRank();

    AmanziMesh::MeshFactory meshfac(comm);
    AmanziMesh::MeshFactory meshfac_serial(comm);

    auto domain_mesh = meshfac.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4 * size, 1, 1);
    auto serial_mesh = meshfac_serial.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4 * size, 1, 1);

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
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->Require<CompositeVector, CompositeVectorSpace>(serial_fname, Tags::DEFAULT, serial_fname)
      .SetMesh(serial_mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    S->Setup();
    S->set_cycle(0);
    S->set_time(0.0);

    // Set values
    int cell_index_list[] = { 0, 1, 2, 3 };
    double cell_values[] = { 10.0, 20.0, 30.0, 40.0 };

    MultiVector_type& mf =
      *S->GetW<CompositeVector>("my_field", Tags::DEFAULT, "my_field").getComponent("cell");
    AMANZI_ASSERT(mf.getLocalLength() == 4);

    MultiVector_type& sf =
      *S->GetW<CompositeVector>(serial_fname, Tags::DEFAULT, serial_fname).getComponent("cell");
    AMANZI_ASSERT(sf.getLocalLength() == 4);


    for (int i = 0; i != 4; ++i) {
      cell_values[i] += rank * 40;
      mf.replaceLocalValue(cell_index_list[i], 0, cell_values[i]);
      sf.replaceLocalValue(cell_index_list[i], 0, cell_values[i]);
    }

    // make sure we can write them all... for real!
    Teuchos::ParameterList chkplist("checkpoint");
    chkplist.set("single file", false);
    Checkpoint chkp(chkplist, *S);
    chkp.write(*S);

    // make sure we can read them all
    State S2;
    S2.RegisterDomainMesh(domain_mesh);
    S2.RegisterMesh(serial_dname, serial_mesh);

    S2.Require<CompositeVector, CompositeVectorSpace>("my_field", Tags::DEFAULT, "my_field")
      .SetMesh(domain_mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S2.Require<CompositeVector, CompositeVectorSpace>(serial_fname, Tags::DEFAULT, serial_fname)
      .SetMesh(domain_mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S2.Setup();

    Checkpoint chkp2("checkpoint00000", comm);
    chkp2.read(S2);

    S2.GetW<CompositeVector>("my_field", Tags::DEFAULT, "my_field")
      .update(-1.0, S->Get<CompositeVector>("my_field"), 1.0);
    double diff = S2.Get<CompositeVector>("my_field").norm2();
    CHECK_CLOSE(0.0, diff, 1.0e-10);

    S2.GetW<CompositeVector>(serial_fname, Tags::DEFAULT, serial_fname)
      .update(-1.0, S->Get<CompositeVector>(serial_fname), 1.0);
    diff = S2.Get<CompositeVector>(serial_fname).norm2();
    CHECK_CLOSE(0.0, diff, 1.0e-10);
  }
}
