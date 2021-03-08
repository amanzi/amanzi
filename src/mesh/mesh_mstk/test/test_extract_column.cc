#include <UnitTest++.h>
#include <fstream>

#include "../Mesh_MSTK.hh"

#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"

// Extract a column from a mesh
TEST(Extract_Column_MSTK)
{
  auto comm = Amanzi::getDefaultComm();
  Teuchos::ParameterList reg_spec; // no regions declared here
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Generate a mesh consisting of 3x3x3 elements 
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0,0,0,1,1,1,3,3,3,comm,gm));

  CHECK_EQUAL(1, mesh->build_columns());
  CHECK_EQUAL(9, mesh->num_columns());

  int cell0 = 0;
  int colid = mesh->column_ID(cell0);
  Amanzi::AmanziMesh::Entity_ID_List const& cell_list = mesh->cells_of_column(colid);

  CHECK_EQUAL(3,cell_list.size());

  // check we are not doubling up on some cells (catches previous bug)
  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      if (i != j) {
        CHECK(cell_list[i] != cell_list[j]);
      }
    }
  }

  Amanzi::AmanziMesh::Entity_ID_List const& face_list = mesh->faces_of_column(colid);
  CHECK_EQUAL(4,face_list.size());

  // check we are not doubling up on some faces (catches previous bug)
  for (int i=0; i!=4; ++i) {
    for (int j=0; j!=4; ++j) {
      if (i != j) {
        CHECK(face_list[i] != face_list[j]);
      }
    }
  }

  // construct a column mesh by extracting from mesh
  Amanzi::AmanziMesh::Mesh_MSTK column_mesh(mesh,cell_list, Amanzi::AmanziMesh::Entity_kind::CELL,false,Amanzi::getCommSelf());
  
  // Number of cells in column mesh
  int ncells_col = column_mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                            Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(3,ncells_col);

  // Number of faces in the column mesh
  int nfaces_col = column_mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE,
                                            Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(16,nfaces_col);
  
  // Check that their parents are as expected
  for (int i = 0; i < ncells_col; ++i) {
    int parent_cell = column_mesh.getEntityParent(Amanzi::AmanziMesh::Entity_kind::CELL,i);
    CHECK_EQUAL(cell_list[i], parent_cell);
  }
}


// Do the same but for an exodus mesh, with sets.
TEST(Extract_Column_MSTK_SETS)
{
  std::string filename("test/hex_3x3x3_sets.exo");
  auto comm = Amanzi::getDefaultComm();
  Teuchos::ParameterList parameterlist;
 
  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved
  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions"); 
  
  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def = top_surface.sublist("region: labeled set");
  top_surface_def.set<std::string>("label","106");
  top_surface_def.set<std::string>("file",filename.c_str());
  top_surface_def.set<std::string>("format","Exodus II");
  top_surface_def.set<std::string>("entity","face");

  Teuchos::ParameterList& side_surface = reg_spec.sublist("Side Surface");
  Teuchos::ParameterList& side_surface_def = side_surface.sublist("region: labeled set");
  side_surface_def.set<std::string>("label","102");
  side_surface_def.set<std::string>("file",filename.c_str());
  side_surface_def.set<std::string>("format","Exodus II");
  side_surface_def.set<std::string>("entity","face");

  Teuchos::ParameterList& r1_surface = reg_spec.sublist("Region 1");
  Teuchos::ParameterList& r1_surface_def = r1_surface.sublist("region: labeled set");
  r1_surface_def.set<std::string>("label","30000");
  r1_surface_def.set<std::string>("file",filename.c_str());
  r1_surface_def.set<std::string>("format","Exodus II");
  r1_surface_def.set<std::string>("entity","cell");
  
  Teuchos::ParameterList& r2_surface = reg_spec.sublist("Region 2");
  Teuchos::ParameterList& r2_surface_def = r2_surface.sublist("region: labeled set");
  r2_surface_def.set<std::string>("label","20000");
  r2_surface_def.set<std::string>("file",filename.c_str());
  r2_surface_def.set<std::string>("format","Exodus II");
  r2_surface_def.set<std::string>("entity","cell");

  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(filename.c_str(), comm, gm));
  CHECK_EQUAL(1, mesh->build_columns());
  CHECK_EQUAL(9, mesh->num_columns());

  // make sure we can get sets on the mesh
  Amanzi::AmanziMesh::Entity_ID_List set_ids;
  mesh->get_set_entities("Region 1", Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_type::ALL, &set_ids);
  CHECK_EQUAL(9, set_ids.size());

  int ncells = 3;
  Amanzi::AmanziMesh::Entity_ID_List const& cell_list = mesh->cells_of_column(0);
  CHECK_EQUAL(ncells,cell_list.size());

  Amanzi::AmanziMesh::Entity_ID_List const& face_list = mesh->faces_of_column(0);
  CHECK_EQUAL(ncells+1,face_list.size());

  // construct a column mesh by extracting from mesh
  Amanzi::AmanziMesh::Mesh_MSTK column_mesh(mesh, cell_list, Amanzi::AmanziMesh::Entity_kind::CELL,false,Amanzi::getCommSelf());
  
  // Number of cells in column mesh
  int ncells_col = column_mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                            Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(ncells,ncells_col);

  // Number of faces in the column mesh
  int nfaces_col = column_mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE,
                                            Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(5*ncells + 1, nfaces_col);
  
  // Check that their parents are as expected
  for (int i = 0; i < ncells_col; ++i) {
    int parent_cell = column_mesh.getEntityParent(Amanzi::AmanziMesh::Entity_kind::CELL,i);
    CHECK_EQUAL(cell_list[i], parent_cell);
  }

  // check we can still get sets
  Amanzi::AmanziMesh::Entity_ID_List set_ids2;
  bool is_valid = column_mesh.valid_set_name("Region 1", Amanzi::AmanziMesh::Entity_kind::CELL);
  CHECK(is_valid);
  column_mesh.get_set_entities("Region 1", Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_type::ALL, &set_ids2);
  CHECK_EQUAL(1, set_ids2.size());

  set_ids2.clear();
  is_valid = column_mesh.valid_set_name("Top Surface", Amanzi::AmanziMesh::Entity_kind::FACE);
  CHECK(is_valid);
  column_mesh.get_set_entities("Top Surface", Amanzi::AmanziMesh::Entity_kind::FACE, Amanzi::AmanziMesh::Parallel_type::ALL, &set_ids2);
  CHECK_EQUAL(1, set_ids2.size());

  set_ids2.clear();
  is_valid = column_mesh.valid_set_name("Side Surface", Amanzi::AmanziMesh::Entity_kind::FACE);
  CHECK(is_valid);
  column_mesh.get_set_entities("Side Surface", Amanzi::AmanziMesh::Entity_kind::FACE, Amanzi::AmanziMesh::Parallel_type::ALL, &set_ids2);
  CHECK_EQUAL(3, set_ids2.size());
}

