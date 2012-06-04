#include "UnitTest++.h"
#include "../cgns_mesh.hh"
#include "Mesh_simple.hh"


TEST(CGNS) {

  using namespace std;

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  string cgns_meshfile = "test_mesh.cgns";
  string cgns_datafile1 = "test_data.cgns";
  string cgns_fullfile = "test_cgns_full.cgns";

  Amanzi::AmanziMesh::Mesh_simple Mesh (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 1, 1, comm);
    
  unsigned int num_nodes = Mesh.num_entities(Amanzi::AmanziMesh::NODE, 
                                               Amanzi::AmanziMesh::OWNED);
  unsigned int num_cells = Mesh.num_entities(Amanzi::AmanziMesh::CELL, 
                                               Amanzi::AmanziMesh::OWNED);

  Teuchos::RCP<Epetra_Vector> node_quantity;
  Teuchos::RCP<Epetra_Vector> cell_quantity;
  Teuchos::RCP<Epetra_Vector> fake_pressure;

  // Setup node quantity
  int node_index_list[] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  double node_values[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120};
  node_quantity = Teuchos::rcp( new Epetra_Vector(Mesh.node_map(false)));
  node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);
  
  // Setup cell quantity
  int cell_index_list[] = {0, 1, 2, 3};
  double cell_values[] = {10, 20, 30, 40};
  cell_quantity = Teuchos::rcp( new Epetra_Vector(Mesh.cell_map(false))); 
  cell_quantity->ReplaceGlobalValues(4, cell_values, cell_index_list);
  
  // Setup second cell quantity -- called fake pressure
  double fake_values[] = {9, 8, 7, 6};
  fake_pressure = Teuchos::rcp( new Epetra_Vector(Mesh.cell_map(false))); 
  fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);  

  // Write a file which contains both mesh and data.
  Amanzi::CGNS::create_mesh_file(Mesh, cgns_fullfile);
  Amanzi::CGNS::open_data_file(cgns_fullfile);
     
  double time = 0.0;
  for (int i=0; i<15; i++) {
      Amanzi::CGNS::create_timestep(time, i, Amanzi::AmanziMesh::CELL);
      Amanzi::CGNS::write_field_data(*cell_quantity, "cell_quantity");
      Amanzi::CGNS::write_field_data(*fake_pressure, "pressure");
      
      // advance time and values
      time += 1.0;
      for (int j=0; j<4; j++) {
	  cell_values[j] += 10.0;
	  fake_values[j] += 1.0;
      }
      cell_quantity->ReplaceGlobalValues(4, cell_values, cell_index_list); 
      fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list); 
  }
    
  // close file
  Amanzi::CGNS::close_data_file();

}

