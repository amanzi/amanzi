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

#include "../OutputSilo.hh"

TEST(SILO_STRUCTURED)
{
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh_MSTK> Mesh =
    Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, comm));
  int ncells =
    Mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  int nnodes =
    Mesh->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED);

  Teuchos::RCP<Epetra_Vector> node_quantity;
  Teuchos::RCP<Epetra_Vector> cell_quantity;
  Teuchos::RCP<Epetra_Vector> fake_pressure;

  // Setup node quantity
  node_quantity = Teuchos::rcp(new Epetra_Vector(Mesh->node_map(false)));
  for (int n = 0; n != nnodes; ++n) {
    Amanzi::AmanziGeometry::Point nc;
    Mesh->node_get_coordinates(n, &nc);
    (*node_quantity)[n] = nc[0];
  }

  // Setup cell quantity
  cell_quantity = Teuchos::rcp(new Epetra_Vector(Mesh->cell_map(false)));
  for (int c = 0; c != ncells; ++c) { (*cell_quantity)[c] = Mesh->cell_centroid(c)[0]; }

  // Setup second cell quantity -- called fake pressure
  fake_pressure = Teuchos::rcp(new Epetra_Vector(Mesh->cell_map(false)));
  for (int c = 0; c != ncells; ++c) {
    if (c % 2 == 0) {
      (*fake_pressure)[c] = Mesh->cell_centroid(c)[0];
    } else {
      (*fake_pressure)[c] = 0.;
    }
  }

  // Write a file which contains both mesh and data.
  Teuchos::ParameterList plist;
  std::string fnb = std::string("amanzi_vis_np") + std::to_string(comm->NumProc());
  plist.set("file name base", fnb);
  Amanzi::OutputSilo io(plist, Mesh, true, false);

  double time = 0.0;

  int NITS = 15;
  for (int i = 0; i < NITS; i++) {
    std::cout << "iteration... " << i << std::endl;

    // write time step data
    io.InitializeCycle(time, i, "");
    io.WriteVector(*cell_quantity, "cell_quantity", Amanzi::AmanziMesh::CELL);
    io.WriteVector(*fake_pressure, "pressure", Amanzi::AmanziMesh::CELL);
    io.WriteVector(*node_quantity, "node_quantity", Amanzi::AmanziMesh::NODE);

    // advance time and values
    time += 2.0;
    double scalar = ((double)i + 2) / (i + 1);
    node_quantity->Scale(scalar);
    cell_quantity->Scale(scalar);
    fake_pressure->Scale(scalar);
    fake_pressure->Print(std::cout);

    // close file
    io.FinalizeCycle();
  }
}


TEST(SILO_POLYGONAL)
{
  auto comm = Amanzi::getDefaultComm();
  if (comm->NumProc() > 1) return; // this test only in serial

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh_MSTK> Mesh =
    Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK("./test/four_polygon.exo", comm));

  Teuchos::RCP<Epetra_Vector> node_quantity;
  Teuchos::RCP<Epetra_Vector> cell_quantity;
  Teuchos::RCP<Epetra_Vector> fake_pressure;

  // Setup node quantity
  int node_index_list[] = { 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
  double node_values[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120 };
  node_quantity = Teuchos::rcp(new Epetra_Vector(Mesh->node_map(false)));
  node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);

  // Setup cell quantity
  cell_quantity = Teuchos::rcp(new Epetra_Vector(Mesh->cell_map(false)));
  for (int c = 0; c != cell_quantity->MyLength(); ++c) { (*cell_quantity)[c] = 10.0 * c; }

  // Setup second cell quantity -- called fake pressure
  int cell_index_list[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  double cell_values[] = { 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0 };
  double fake_values[] = { 9, 8, 7, 6 };
  fake_pressure = Teuchos::rcp(new Epetra_Vector(Mesh->cell_map(false)));
  fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);

  // Write a file which contains both mesh and data.
  Teuchos::ParameterList plist;
  std::string fnb = std::string("polygonal_vis_np") + std::to_string(comm->NumProc());
  plist.set("file name base", fnb);
  Amanzi::OutputSilo io(plist, Mesh, true, false);

  double time = 0.0;

  for (int i = 0; i < 15; i++) {
    std::cout << "iteration... " << i << std::endl;

    cell_quantity->ReplaceGlobalValues(8, cell_values, cell_index_list);
    fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);
    node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);

    // write time step data
    io.InitializeCycle(time, i, "");
    io.WriteVector(*cell_quantity, "cell_quantity", Amanzi::AmanziMesh::CELL);
    io.WriteVector(*fake_pressure, "pressure", Amanzi::AmanziMesh::CELL);
    io.WriteVector(*node_quantity, "node_quantity", Amanzi::AmanziMesh::NODE);

    // advance time and values
    time += 2.0;
    for (int j = 0; j < 8; j++) { cell_values[j] += 10.0; }
    for (int j = 0; j < 4; j++) { fake_values[j] += 1.0; }
    for (int j = 0; j < 12; j++) { node_values[j] += 10.0; }

    // close file
    io.FinalizeCycle();
  }
}
