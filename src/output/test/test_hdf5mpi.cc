/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#define MSTK_HAVE_MPI 1

#include "UnitTest++.h"
#include "HDF5_MPI.hh"
#include "Mesh_MSTK.hh"
TEST(HDF5_MPI)
{
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();

  std::string hdf5_meshfile = "new_mesh_mpi";
  std::string hdf5_datafile1 = "new_data_mpi";
  std::string hdf5_datafile2 = "new_restart_mpi";

  auto Mesh =
    Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, comm));

  Teuchos::RCP<Epetra_Vector> node_quantity;
  Teuchos::RCP<Epetra_Vector> cell_quantity;
  Teuchos::RCP<Epetra_Vector> fake_pressure;

  // Setup node quantity
  int node_index_list[] = { 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
  double node_values[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120 };
  node_quantity = Teuchos::rcp(new Epetra_Vector(Mesh->node_map(false)));
  node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);

  // Setup cell quantity
  int cell_index_list[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  double cell_values[] = { 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0 };
  cell_quantity = Teuchos::rcp(new Epetra_Vector(Mesh->cell_map(false)));
  cell_quantity->ReplaceGlobalValues(8, cell_values, cell_index_list);

  // Setup second cell quantity -- called fake pressure
  double fake_values[] = { 9, 8, 7, 6 };
  fake_pressure = Teuchos::rcp(new Epetra_Vector(Mesh->cell_map(false)));
  fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);

  // Write a file which contains both mesh and data.
  Amanzi::HDF5_MPI* viz_output = new Amanzi::HDF5_MPI(comm);
  viz_output->setTrackXdmf(true);
  viz_output->createMeshFile(Mesh, hdf5_meshfile);
  viz_output->createDataFile(hdf5_datafile1);

  Amanzi::HDF5_MPI* restart_output = new Amanzi::HDF5_MPI(comm);
  restart_output->setTrackXdmf(false);
  restart_output->createDataFile(hdf5_datafile2);
  // You can add mesh data to restart file, but is not necessary for valid restart
  restart_output->createMeshFile(Mesh, hdf5_datafile2);

  double time = 0.0;
  int cycle = 0;

  for (int i = 0; i < 15; i++) {
    std::cout << "iteration... " << i << std::endl;

    cell_quantity->ReplaceGlobalValues(8, cell_values, cell_index_list);
    fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);
    node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);

    // write time step data
    viz_output->open_h5file();
    viz_output->createTimestep(time, i, "");
    viz_output->writeCellDataReal(*cell_quantity, "cell_quantity");
    viz_output->writeCellDataReal(*fake_pressure, "pressure");
    viz_output->writeNodeDataReal(*node_quantity, "node_quantity");

    // advance time and values
    time += 2.0;
    cycle += 1;
    for (int j = 0; j < 8; j++) { cell_values[j] += 10.0; }
    for (int j = 0; j < 4; j++) { fake_values[j] += 1.0; }
    for (int j = 0; j < 12; j++) { node_values[j] += 10.0; }

    // close file
    viz_output->endTimestep();
    viz_output->close_h5file();
  }


  // write out restart

  restart_output->open_h5file();
  restart_output->writeAttrReal(time, "time");
  restart_output->writeAttrInt(cycle, "cycle");
  restart_output->writeAttrString("string value", "attr name");
  restart_output->writeCellDataReal(*cell_quantity, "cell_quantity");
  restart_output->writeCellDataReal(*fake_pressure, "pressure");
  restart_output->writeNodeDataReal(*node_quantity, "node_quantity");
  restart_output->close_h5file();

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

  restart_output->writeDataString(strArray,num_wstrs,"string_dataset");
  */

  std::cout << "O.K." << std::endl;

  delete viz_output;
  delete restart_output;

  // test reading data back
  std::cout << "E>> create restart_input with file " << hdf5_datafile2 << ".h5" << std::endl;
  Amanzi::HDF5_MPI* restart_input = new Amanzi::HDF5_MPI(comm, hdf5_datafile2 + ".h5");
  std::cout << hdf5_datafile2 + ".h5" << std::endl;

  restart_input->open_h5file();

  double newtime;
  int newcycle;
  std::string newstring;

  restart_input->readAttrReal(newtime, "time");
  std::cout << "E>> read back attribute time = " << newtime << std::endl;
  restart_input->readAttrInt(newcycle, "cycle");
  std::cout << "E>> read back attribute cycle = " << newcycle << std::endl;
  restart_input->readAttrString(newstring, "attr name");
  std::cout << "E>> read back attribute string = " << newstring << std::endl;
  std::cout << "E>> compare results" << std::endl;
  std::cout << "E>> original:" << std::endl << *cell_quantity;

  Teuchos::RCP<Epetra_Vector> read_quantity;
  read_quantity = Teuchos::rcp(new Epetra_Vector(Mesh->cell_map(false)));
  std::cout << std::endl;
  restart_input->readData(*read_quantity, "cell_quantity");

  std::cout << "E>> read back:" << std::endl << *read_quantity;
  std::cout << "E>> cell map:" << std::endl << Mesh->cell_map(false);

  // reading back string dataset
  /*
  char **strBack;
  int num_rstrs = 0;
  restart_input->readDataString(&strBack, &num_rstrs, "string_dataset");
  std::cout << "E>> reading back string dataset["<<num_rstrs<<"]: " << std::endl;
  for (int i=0 ; i<num_rstrs; i++) {
    std::cout << "    " << strBack[i] << std::endl;
  }
  */

  restart_input->close_h5file();

  delete restart_input;


#endif
}
