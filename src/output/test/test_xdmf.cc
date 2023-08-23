/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "UnitTest++.h"

#define MSTK_HAVE_MPI
#include "Mesh_MSTK.hh"

#include "OutputXDMF.hh"

TEST(XDMF)
{
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh_MSTK> Mesh_mstk =
    Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, comm));
  auto Mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh(Mesh_mstk, Teuchos::null));

  // unsigned int num_cells = Mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
  //         Amanzi::AmanziMesh::Parallel_type::OWNED);

  Teuchos::RCP<Epetra_Vector> node_quantity;
  Teuchos::RCP<Epetra_Vector> cell_quantity;
  Teuchos::RCP<Epetra_Vector> fake_pressure;

  // Setup node quantity
  int node_index_list[] = { 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
  double node_values[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120 };
  node_quantity =
    Teuchos::rcp(new Epetra_Vector(Mesh->getMap(Amanzi::AmanziMesh::Entity_kind::NODE, false)));
  node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);

  // Setup cell quantity
  int cell_index_list[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  double cell_values[] = { 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0 };
  cell_quantity =
    Teuchos::rcp(new Epetra_Vector(Mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false)));
  cell_quantity->ReplaceGlobalValues(8, cell_values, cell_index_list);

  // Setup second cell quantity -- called fake pressure
  double fake_values[] = { 9, 8, 7, 6 };
  fake_pressure =
    Teuchos::rcp(new Epetra_Vector(Mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false)));
  fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);

  // Write a file which contains both mesh and data.
  Teuchos::ParameterList plist;
  Amanzi::OutputXDMF io(plist, Mesh, true, false, true);

  double time = 0.0;

  for (int i = 0; i < 15; i++) {
    std::cout << "iteration... " << i << std::endl;

    cell_quantity->ReplaceGlobalValues(8, cell_values, cell_index_list);
    fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);
    node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);

    // write time step data
    io.InitializeCycle(time, i, "");
    io.WriteVector(*cell_quantity, "cell_quantity", Amanzi::AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*fake_pressure, "pressure", Amanzi::AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*node_quantity, "node_quantity", Amanzi::AmanziMesh::Entity_kind::CELL);

    // advance time and values
    time += 2.0;
    for (int j = 0; j < 8; j++) { cell_values[j] += 10.0; }
    for (int j = 0; j < 4; j++) { fake_values[j] += 1.0; }
    for (int j = 0; j < 12; j++) { node_values[j] += 10.0; }

    // close file
    io.FinalizeCycle();
  }

#endif
}
