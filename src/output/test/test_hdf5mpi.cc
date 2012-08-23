#include "UnitTest++.h"
#include "../hdf5mpi_mesh.hh"
#if HAVE_STK_MESH
#include "Mesh_STK.hh"
#endif
TEST(HDF5_MPI)  {
  
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
//#else
//  Epetra_SerialComm *comm = new Epetra_SerialComm();
//#endif
    
  std::string hdf5_meshfile  = "new_mesh_mpi";
  std::string hdf5_datafile1 = "new_data_mpi";
  std::string hdf5_datafile2 = "new_restart_mpi";
  
  //Teuchos::RCP<Amanzi::AmanziMesh::Mesh_STK> 
  //  Mesh(new Amanzi::AmanziMesh::Mesh_STK(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 1,
  //                                        1, comm));
  Amanzi::AmanziMesh::Mesh_STK Mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1, comm);

  unsigned int num_nodes = Mesh.num_entities(Amanzi::AmanziMesh::NODE, 
                                                Amanzi::AmanziMesh::OWNED);
  unsigned int num_cells = Mesh.num_entities(Amanzi::AmanziMesh::CELL, 
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

  // Setup up mesh region
  Epetra_Map regMap(3, 0, *comm);
  double region_cells[] = {0,2,4};
  double region_cells2[] = {4,5,7};
  int region_index_list[] = {0,1,2};
  Teuchos::RCP<Epetra_Vector> mesh_region1, mesh_region2;
  mesh_region1 = Teuchos::rcp(new Epetra_Vector(regMap,false));
  mesh_region1->ReplaceGlobalValues(3, region_cells, region_index_list);
  std::string region_name1, region_name2;
  region_name1 = "Region1";
  mesh_region2 = Teuchos::rcp(new Epetra_Vector(regMap,false));
  mesh_region2->ReplaceGlobalValues(3, region_cells2, region_index_list);
  region_name2 = "Region2";

  // Write a file which contains both mesh and data.
  Amanzi::HDF5_MPI *viz_output = new Amanzi::HDF5_MPI(*comm);
  viz_output->setTrackXdmf(true);
  viz_output->createMeshFile(Mesh, hdf5_meshfile);
  viz_output->createDataFile(hdf5_datafile1);
  viz_output->writeMeshRegion(Mesh, *mesh_region1, region_name1);
  viz_output->writeMeshRegion(Mesh, *mesh_region2, region_name2);
  
  // Create restart file
  Amanzi::HDF5_MPI *restart_output = new Amanzi::HDF5_MPI(*comm);
  restart_output->setTrackXdmf(false);
  restart_output->createDataFile(hdf5_datafile2);
  // You can add mesh data to restart file, but is not necessary for valid restart
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
  restart_output->writeAttrString("string value", "attr name");
  restart_output->writeCellDataReal(*cell_quantity, "cell_quantity");
  restart_output->writeCellDataReal(*fake_pressure, "pressure");
  restart_output->writeNodeDataReal(*node_quantity, "node_quantity");
  
  // write out string dataset
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
  
  delete viz_output;
  delete restart_output;
  
  // test reading data back
  cout << "E>> create restart_input with file " << hdf5_datafile2 << ".h5" << endl;
  Amanzi::HDF5_MPI *restart_input = new Amanzi::HDF5_MPI(*comm,hdf5_datafile2+".h5");
  double newtime;
  int newcycle;
  std::string newstring;
  restart_input->readAttrReal(newtime,"time");
  cout << "E>> read back attribute time = " << newtime << endl;
  restart_input->readAttrInt(newcycle,"cycle");
  cout << "E>> read back attribute cycle = " << newcycle << endl;
  restart_input->readAttrString(newstring,"attr name");
  cout << "E>> read back attribute string = " << newstring.c_str() << endl;
  cout << "E>> compare results" << endl;
  cout << "E>> original:" << endl << *cell_quantity;
  Teuchos::RCP<Epetra_Vector> read_quantity;
  read_quantity = Teuchos::rcp(new Epetra_Vector(Mesh.cell_map(false)));
  cout << endl;
  restart_input->readData(*read_quantity, "cell_quantity");
  
  cout << "E>> read back:" << endl << *read_quantity;
  cout << "E>> cell map:" << endl << Mesh.cell_map(false);

  // reading back string dataset
  char **strBack;
  int num_rstrs = 0;
  restart_input->readDataString(&strBack, &num_rstrs, "string_dataset");
  cout << "E>> reading back string dataset["<<num_rstrs<<"]: " << endl;
  for (int i=0 ; i<num_rstrs; i++) {
    cout << "    " << strBack[i] << endl;
  }
  
  delete restart_input;

#endif
}

