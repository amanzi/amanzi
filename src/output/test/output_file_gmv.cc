/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "UnitTest++.h"
#include "../GMVMesh.hh"
#include "MeshFactory.hh"
#include "AmanziVector.hh"

TEST(GMV)
{
  using namespace Amanzi;
  auto comm = Amanzi::getDefaultComm();

  std::string gmv_meshfile = "test_mesh.gmv";
  std::string gmv_datafile1 = "test_gmv1.gmv";
  std::string gmv_fullfile = "test_gmv_full.gmv";

  Amanzi::AmanziMesh::Preference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);

  auto mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 1, 1);

  unsigned int num_nodes =
    mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE, Amanzi::AmanziMesh::Parallel_kind::OWNED);
  unsigned int num_cells =
    mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_kind::OWNED);

  // Setup node quantity
  std::vector<int> node_index_list{ 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
  std::vector<double> node_values{ 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120 };
  auto node_quantity = Teuchos::rcp(new Vector_type(mesh->getMap(Amanzi::AmanziMesh::Entity_kind::NODE,false)));
  for (int i=0; i!=node_index_list.size(); ++i)
    node_quantity->replaceGlobalValue(node_index_list[i], node_values[i]);

  // Setup cell quantity
  std::vector<int> cell_index_list{ 0, 1, 2, 3 };
  std::vector<double> cell_values{ 10, 20, 30, 40 };
  auto cell_quantity = Teuchos::rcp(new Vector_type(mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL,false)));
  for (int i=0; i!=cell_index_list.size(); ++i)
    cell_quantity->replaceGlobalValue(cell_index_list[i], cell_values[i]);

  // Setup second cell quantity -- called fake pressure
  std::vector<double> fake_values{ 9, 8, 7, 6 };
  auto fake_pressure = Teuchos::rcp(new Vector_type(mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL,false)));
  for (int i=0; i!=cell_index_list.size(); ++i)
    fake_pressure->replaceGlobalValue(cell_index_list[i], fake_values[i]);


  Amanzi::GMV::create_mesh_file(*(mesh.get()), gmv_meshfile);
  Amanzi::GMV::open_data_file(gmv_meshfile, gmv_datafile1, num_nodes, num_cells);
  Amanzi::GMV::start_data();
  Amanzi::GMV::write_node_data(*node_quantity, "node_quantity");
  Amanzi::GMV::write_cell_data(*cell_quantity, "cell_quantity");
  Amanzi::GMV::write_cell_data(*fake_pressure, "pressure");
  Amanzi::GMV::close_data_file();

  // Write a file which contains both mesh and data.
  Amanzi::GMV::open_data_file(*(mesh.get()), gmv_fullfile);
  Amanzi::GMV::start_data();
  Amanzi::GMV::write_node_data(*node_quantity, "node_quantity");
  Amanzi::GMV::write_cell_data(*cell_quantity, "cell_quantity");
  Amanzi::GMV::write_cell_data(*fake_pressure, "pressure");
  Amanzi::GMV::close_data_file();
}
