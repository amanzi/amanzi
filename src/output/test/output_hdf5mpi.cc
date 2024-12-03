/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "UnitTest++.h"
#include "HDF5_MPI.hh"
#include "Mesh_MSTK.hh"
TEST(HDF5_MPI)
{
  auto comm = Amanzi::getDefaultComm();

  std::string hdf5_meshfile = "new_mesh_mpi";
  std::string hdf5_datafile1 = "new_data_mpi";
  std::string hdf5_datafile2 = "new_restart_mpi";

  auto Mesh_mstk =
    Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, comm));

  auto Mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh(
    Mesh_mstk, Teuchos::rcp(new Amanzi::AmanziMesh::MeshAlgorithms()), Teuchos::null));

  // Setup node quantity
  int node_index_list[] = { 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
  double node_values[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120 };
  Epetra_Vector node_quantity(Mesh->getMap(Amanzi::AmanziMesh::Entity_kind::NODE, false));
  node_quantity.ReplaceGlobalValues(12, node_values, node_index_list);

  // Setup cell quantity
  int cell_index_list[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  double cell_values[] = { 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0 };
  Epetra_Vector cell_quantity(Mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false));
  cell_quantity.ReplaceGlobalValues(8, cell_values, cell_index_list);

  // Setup second cell quantity -- called fake pressure
  double fake_values[] = { 9, 8, 7, 6 };
  Epetra_Vector fake_pressure(Mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false));
  fake_pressure.ReplaceGlobalValues(4, fake_values, cell_index_list);

  double time = 0.0;
  int cycle = 0;

  {
    // Mock visualization write
    // Write a file which contains both mesh and data.
    Amanzi::HDF5_MPI viz_output(comm);
    viz_output.setTrackXdmf(true);
    viz_output.createMeshFile(Mesh, hdf5_meshfile);
    viz_output.createDataFile(hdf5_datafile1);

    for (int i = 0; i < 15; i++) {
      std::cout << "iteration... " << i << std::endl;

      cell_quantity.ReplaceGlobalValues(8, cell_values, cell_index_list);
      fake_pressure.ReplaceGlobalValues(4, fake_values, cell_index_list);
      node_quantity.ReplaceGlobalValues(12, node_values, node_index_list);

      // write time step data
      viz_output.open_h5file();
      viz_output.createTimestep(time, i, "");
      viz_output.writeCellDataReal(cell_quantity, "cell_quantity");
      viz_output.writeCellDataReal(fake_pressure, "pressure");
      viz_output.writeNodeDataReal(node_quantity, "node_quantity");

      // advance time and values
      time += 2.0;
      cycle += 1;
      for (int j = 0; j < 8; j++) { cell_values[j] += 10.0; }
      for (int j = 0; j < 4; j++) { fake_values[j] += 1.0; }
      for (int j = 0; j < 12; j++) { node_values[j] += 10.0; }

      // close file
      viz_output.endTimestep();
      viz_output.close_h5file();
    }
  }

  {
    // Mock checkpoint write
    Amanzi::HDF5_MPI restart_output(comm);
    restart_output.setTrackXdmf(false);
    restart_output.createDataFile(hdf5_datafile2);
    // You can add mesh data to restart file, but is not necessary for valid restart
    restart_output.createMeshFile(Mesh, hdf5_datafile2);

    // write out restart
    restart_output.open_h5file();
    restart_output.writeAttrReal(time, "time");
    restart_output.writeAttrInt(cycle, "cycle");
    restart_output.writeAttrString("string value", "attr name");
    restart_output.writeCellDataReal(cell_quantity, "cell_quantity");
    restart_output.writeCellDataReal(fake_pressure, "pressure");
    restart_output.writeNodeDataReal(node_quantity, "node_quantity");
    restart_output.close_h5file();

    // write out string dataset
    /*
      int num_wstrs = 5;
      char **strArray;
      strArray = (char**) malloc(5*sizeof(char*));
      for (int i=0; i<num_wstrs; i++) {
      strArray[i] = (char *)malloc(MAX_STRING_LENGTH*sizeof(char));
      }
      sprintf(strArray[0], "Calcium");
      sprintf(strArray[1], "Magnesium");
      sprintf(strArray[2], "Uranium");
      sprintf(strArray[3], "Unobtainium");
      sprintf(strArray[4], "My Favorite Mineral in the Whole World");

      restart_output.writeDataString(strArray,num_wstrs,"string_dataset");
    */
  }

  {
    // test reading data back
    Amanzi::HDF5_MPI restart_input(comm, hdf5_datafile2 + ".h5");
    std::cout << hdf5_datafile2 + ".h5" << std::endl;

    restart_input.open_h5file();

    double newtime;
    int newcycle;
    std::string newstring;

    // read and compare attributes
    restart_input.readAttrReal(newtime, "time");
    CHECK_EQUAL(time, newtime);

    restart_input.readAttrInt(newcycle, "cycle");
    CHECK_EQUAL(cycle, newcycle);

    restart_input.readAttrString(newstring, "attr name");
    CHECK(std::string("string value") == newstring);

    // read and compare vectors
    Epetra_Vector read_quantity(Mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false));
    restart_input.readData(read_quantity, "cell_quantity");
    read_quantity.Update(-1, cell_quantity, 1);

    double norm(1.0);
    read_quantity.NormInf(&norm);
    CHECK_CLOSE(0.0, norm, 1.e-10);

    // make sure some basic things throw
    // -- read int as real
    newtime = 0.0;
    CHECK_THROW(restart_input.readAttrReal(newtime, "cycle"), Errors::Message);

    // -- read nonexistent attr
    CHECK_THROW(restart_input.readAttrReal(newtime, "nonexistent"), Errors::Message);

    // -- read nonexistent string
    CHECK_THROW(restart_input.readAttrString(newstring, "nonexistent_string"), Errors::Message);

    // -- read incorrect length vector
    CHECK_THROW(restart_input.readData(node_quantity, "cell_quantity"), Errors::Message);

    // -- read nonexistent vector
    CHECK_THROW(restart_input.readData(node_quantity, "nonexistent_data"), Errors::Message);

    restart_input.close_h5file();
  }
}
