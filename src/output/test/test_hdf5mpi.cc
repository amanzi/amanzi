#include "UnitTest++.h"
#include "../hdf5mpi_mesh.hh"
#if HAVE_STK_MESH
#include "Mesh_STK.hh"
#endif
TEST(HDF5_MPI) {
  
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
//#else
//  Epetra_SerialComm *comm = new Epetra_SerialComm();
//#endif
    
  std::string hdf5_meshfile  = "new_mesh";
  std::string hdf5_datafile1 = "new_data";
  std::string hdf5_datafile2 = "new_restart";
  
  //Teuchos::RCP<Amanzi::AmanziMesh::Mesh_STK> 
  //  Mesh(new Amanzi::AmanziMesh::Mesh_STK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 1, 1,
  //                                        comm));
  Amanzi::AmanziMesh::Mesh_STK Mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, comm);

  unsigned int num_nodes = Mesh.count_entities(Amanzi::AmanziMesh::NODE, 
                                                Amanzi::AmanziMesh::OWNED);
  unsigned int num_cells = Mesh.count_entities(Amanzi::AmanziMesh::CELL, 
                                                Amanzi::AmanziMesh::OWNED);

  //Teuchos::RCP<Mesh_maps_base> Mesh(new STK_mesh::Mesh_maps_stk(0.0, 0.0, 0.0,
  //			            1.0, 1.0, 1.0, 4, 1, 1, comm));
  //STK_mesh::Mesh_maps_stk Mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 1, 1, comm);
  //unsigned int num_nodes = Mesh.count_entities(Mesh_data::NODE, OWNED);
  //unsigned int num_cells = Mesh.count_entities(Mesh_data::CELL, OWNED);

  Teuchos::RCP<Epetra_Vector> node_quantity;
  Teuchos::RCP<Epetra_Vector> cell_quantity;
  Teuchos::RCP<Epetra_Vector> fake_pressure;

  // Setup node quantity
  int node_index_list[] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  double node_values[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120};
  node_quantity = Teuchos::rcp(new Epetra_Vector(Mesh.node_map(false)));
  node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);

  // Setup cell quantity
  int cell_index_list[] = {0, 1, 2, 3, 4, 5, 6, 7};
  double cell_values[] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};
  cell_quantity = Teuchos::rcp(new Epetra_Vector(Mesh.cell_map(false)));
  cell_quantity->ReplaceGlobalValues(8, cell_values, cell_index_list);

  // Setup second cell quantity -- called fake pressure
  double fake_values[] = {9, 8, 7, 6};
  fake_pressure = Teuchos::rcp(new Epetra_Vector(Mesh.cell_map(false)));
  fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);

  // Write a file which contains both mesh and data.
  Amanzi::HDF5_MPI *viz_output = new Amanzi::HDF5_MPI(*comm);
  viz_output->setTrackXdmf(true);
  viz_output->createMeshFile(Mesh, hdf5_meshfile);
  viz_output->createDataFile(hdf5_datafile1);
  
  Amanzi::HDF5_MPI *restart_output = new Amanzi::HDF5_MPI(*comm);
  restart_output->setTrackXdmf(false);
  restart_output->createDataFile(hdf5_datafile2);
  restart_output->createMeshFile(Mesh, hdf5_datafile2);

  double time = 0.0;
  int cycle = 0;
  for (int i = 0; i < 15; i++) {
    
    cell_quantity->ReplaceGlobalValues(8, cell_values, cell_index_list);
    fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);
    node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);
    
    // write time step data
    viz_output->createTimestep(time, i);
    viz_output->writeCellDataReal(*cell_quantity, "cell_quantity");
    viz_output->writeCellDataReal(*fake_pressure, "pressure");
    viz_output->writeNodeDataReal(*node_quantity, "node_quantity");

    // advance time and values
    time += 2.0;
    cycle += 1;
    for (int j = 0; j < 8; j++) {
      cell_values[j] += 10.0;
    }
    for (int j = 0; j < 4; j++) {
      fake_values[j] += 1.0;
    }
    for (int j = 0; j < 12; j++) {
      node_values[j] += 10.0;
    }

    // close file
    viz_output->endTimestep();
  }
   
  // write out restart
  restart_output->writeAttrReal(time, "time");
  restart_output->writeAttrInt(cycle, "cycle");
  restart_output->writeCellDataReal(*cell_quantity, "cell_quantity");
  restart_output->writeCellDataReal(*fake_pressure, "pressure");
  restart_output->writeNodeDataReal(*node_quantity, "node_quantity");
  
  // test reading data back
  double newtime;
  int newcycle;
  restart_output->readAttrReal(newtime,"time");
  cout << "E>> read back attribute time = " << newtime << endl;
  restart_output->readAttrInt(newcycle,"cycle");
  cout << "E>> read back attribute cycle = " << newcycle << endl;
  cout << "E>> compare results" << endl;
  cout << "E>> original:" << endl << *cell_quantity;
  Teuchos::RCP<Epetra_Vector> read_quantity;
  read_quantity = Teuchos::rcp(new Epetra_Vector(Mesh.cell_map(false)));
  cout << endl;
  restart_output->readData(*read_quantity, "cell_quantity");
  
  cout << "E>> read back:" << endl << *read_quantity;
  cout << "E>> cell map:" << endl << Mesh.cell_map(false);
  
  delete viz_output;
  delete restart_output;

#endif
}

