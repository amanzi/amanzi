#include <iostream>
#include "stdlib.h"
#include "../gmv_mesh.hh"
#include "Setup_Epetra_example.hh"
#include <stk_util/parallel/Parallel.hpp>

stk::ParallelMachine parallel_machine;

using namespace std;

int test(int argc, char *argv[])
{
  string gmv_meshfile = "test_mesh.gmv";
  string gmv_datafile1 = "test_gmv1.gmv";
  string gmv_fullfile = "test_gmv_full.gmv";
  struct Map_setup stk_map_setup;
  struct Epetra_mesh_setup epetra_setup;

  // Write a mesh file and a data file which uses the mesh file.
  unsigned int num_nodes = stk_map_setup.mesh_map.count_entities(Mesh_data::NODE, STK_mesh::OWNED);
  unsigned int num_cells = stk_map_setup.mesh_map.count_entities(Mesh_data::CELL, STK_mesh::OWNED);
  GMV::create_mesh_file(stk_map_setup.mesh_map, gmv_meshfile);
  GMV::open_data_file(gmv_meshfile, gmv_datafile1, num_nodes, num_cells);
  GMV::write_node_data(*epetra_setup.node_quantity, "node_quantity");
  GMV::write_cell_data(*epetra_setup.cell_quantity, "cell_quantity");
  GMV::write_cell_data(*epetra_setup.fake_pressure, "pressure");
  GMV::close_data_file();

  // Write a file which contains both mesh and data.
  GMV::open_data_file(stk_map_setup.mesh_map, gmv_fullfile);
  GMV::write_node_data(*epetra_setup.node_quantity, "node_quantity");
  GMV::write_cell_data(*epetra_setup.cell_quantity, "cell_quantity");
  GMV::write_cell_data(*epetra_setup.fake_pressure, "pressure");
  GMV::close_data_file();

}

int main(int argc, char *argv[])
{
  parallel_machine = stk::parallel_machine_init (&argc, &argv);
  test(argc, argv);
  stk::parallel_machine_finalize ();
}
