#include "UnitTest++.h"
#include "../GMVMesh.hh"
#include "MeshFactory.hh"


TEST(GMV)
{
  using namespace std;

#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  string gmv_meshfile = "test_mesh.gmv";
  string gmv_datafile1 = "test_gmv1.gmv";
  string gmv_fullfile = "test_gmv_full.gmv";

  Amanzi::AmanziMesh::Preference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> Mesh =
    meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 1, 1);

  unsigned int num_nodes =
    Mesh->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED);
  unsigned int num_cells =
    Mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  Teuchos::RCP<Epetra_Vector> node_quantity;
  Teuchos::RCP<Epetra_Vector> cell_quantity;
  Teuchos::RCP<Epetra_Vector> fake_pressure;

  // Setup node quantity
  int node_index_list[] = { 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
  double node_values[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120 };
  node_quantity = Teuchos::rcp(new Epetra_Vector(Mesh->node_map(false)));
  node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);

  // Setup cell quantity
  int cell_index_list[] = { 1, 2, 3, 4 };
  double cell_values[] = { 10, 20, 30, 40 };
  cell_quantity = Teuchos::rcp(new Epetra_Vector(Mesh->cell_map(false)));
  cell_quantity->ReplaceGlobalValues(4, cell_values, cell_index_list);

  // Setup second cell quantity -- called fake pressure
  double fake_values[] = { 9, 8, 7, 6 };
  fake_pressure = Teuchos::rcp(new Epetra_Vector(Mesh->cell_map(false)));
  fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);


  Amanzi::GMV::create_mesh_file(*(Mesh.get()), gmv_meshfile);
  Amanzi::GMV::open_data_file(gmv_meshfile, gmv_datafile1, num_nodes, num_cells);
  Amanzi::GMV::start_data();
  Amanzi::GMV::write_node_data(*node_quantity, "node_quantity");
  Amanzi::GMV::write_cell_data(*cell_quantity, "cell_quantity");
  Amanzi::GMV::write_cell_data(*fake_pressure, "pressure");
  Amanzi::GMV::close_data_file();

  // Write a file which contains both mesh and data.
  Amanzi::GMV::open_data_file(*(Mesh.get()), gmv_fullfile);
  Amanzi::GMV::start_data();
  Amanzi::GMV::write_node_data(*node_quantity, "node_quantity");
  Amanzi::GMV::write_cell_data(*cell_quantity, "cell_quantity");
  Amanzi::GMV::write_cell_data(*fake_pressure, "pressure");
  Amanzi::GMV::close_data_file();
}
