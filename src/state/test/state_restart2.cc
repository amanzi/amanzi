/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "UnitTest++.h"

#include "Checkpoint.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "State.hh"

SUITE(RESTART2) {

TEST(HDF5_MPI_AND_SERIAL) {
  using namespace Amanzi;

  auto comm = getDefaultComm();
  auto comm_serial = getCommSelf();

  std::string hdf5_datafile1 = "new_data_mpi_serial_parallel";

  int size = comm->NumProc();
  int rank = comm->MyPID();
  auto domain_mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4*size, 1, 1, comm));
  auto serial_mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(0,0,0,1,1,1,4,1,1,comm_serial));

  auto S = Teuchos::rcp(new State());
  S->RegisterDomainMesh(domain_mesh);
  std::stringstream serial_dname_ss;
  serial_dname_ss << "serial_" << rank;

  std::string serial_dname = serial_dname_ss.str();
  Key serial_fname = Keys::getKey(serial_dname, "field");
  S->RegisterMesh(serial_dname, serial_mesh);

  S->RequireField("my_field", "my_field")->SetMesh(domain_mesh)->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(serial_fname, serial_fname)->SetMesh(serial_mesh)->SetComponent("cell", AmanziMesh::CELL, 1);

  S->Setup();
  S->set_cycle(0);
  S->set_time(0.0);

  // Set values
  int cell_index_list[] = {0, 1, 2, 3};
  double cell_values[] = {10.0, 20.0, 30.0, 40.0};
  for (int i=0; i!=4; ++i) {
    cell_values[i] += rank * 40;
  }
  Epetra_MultiVector& mf = *S->GetFieldData("my_field", "my_field")->ViewComponent("cell", false);
  mf(0)->ReplaceMyValues(4, cell_values, cell_index_list);

  Epetra_MultiVector& sf = *S->GetFieldData(serial_fname, serial_fname)->ViewComponent("cell", false);
  sf(0)->ReplaceMyValues(4, cell_values, cell_index_list);

  // make sure we can write them all... for real!
  Teuchos::ParameterList chkplist("checkpoint");
  chkplist.set("single file checkpoint", false);
  Checkpoint chkp(chkplist, *S);
  chkp.Write(*S, 0.);

  // make sure we can read them all
  State S2;
  S2.RegisterDomainMesh(domain_mesh);
  S2.RegisterMesh(serial_dname, serial_mesh);
  S2.RequireField("my_field", "my_field")->SetMesh(domain_mesh)->SetComponent("cell", AmanziMesh::CELL, 1);
  S2.RequireField(serial_fname, serial_fname)->SetMesh(serial_mesh)->SetComponent("cell", AmanziMesh::CELL, 1);
  S2.Setup();

  ReadCheckpoint(S2, "checkpoint00000");

  S2.GetFieldData("my_field", "my_field")->Update(-1., *S->GetFieldData("my_field"), 1.);
  double diff;
  S2.GetFieldData("my_field")->Norm2(&diff);
  CHECK_CLOSE(0.0, diff, 1.e-10);

  S2.GetFieldData(serial_fname, serial_fname)->Update(-1., *S->GetFieldData(serial_fname), 1.);
  diff = 0.;
  S2.GetFieldData(serial_fname)->Norm2(&diff);
  CHECK_CLOSE(0.0, diff, 1.e-10);

}

}
