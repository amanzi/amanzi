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
  Amanzi::AmanziMesh::Mesh_MSTK column_mesh(mesh,cell_list, Amanzi::AmanziMesh::CELL,false,Amanzi::getCommSelf());
  
  // Number of cells in column mesh
  int ncells_col = column_mesh.num_entities(Amanzi::AmanziMesh::CELL,
                                            Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(3,ncells_col);

  // Number of faces in the column mesh
  int nfaces_col = column_mesh.num_entities(Amanzi::AmanziMesh::FACE,
                                            Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(16,nfaces_col);
  
  // Check that their parents are as expected
  for (int i = 0; i < ncells_col; ++i) {
    int parent_cell = column_mesh.entity_get_parent(Amanzi::AmanziMesh::CELL,i);
    CHECK_EQUAL(cell_list[i], parent_cell);
  }
}


// Do the same but for an exodus mesh, with sets.
TEST(Extract_Column_MSTK_SETS)
{
  // construct a column mesh
  auto comm = Amanzi::getDefaultComm();

  std::string infilename = "test/test_extract_column.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK("test/four_polygons.exo",comm,gm));
  CHECK_EQUAL(1, mesh->build_columns());
  CHECK_EQUAL(4, mesh->num_columns());

  // make sure we can get sets on the mesh
  Amanzi::AmanziMesh::Entity_ID_List set_ids;
  mesh->get_set_entities("organic", Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL, &set_ids);
  CHECK_EQUAL(40, set_ids.size());
  
  int ncells = 118;
  
  Amanzi::AmanziMesh::Entity_ID_List const& cell_list = mesh->cells_of_column(0);
  CHECK_EQUAL(ncells,cell_list.size());

  // check we are not doubling up on some cells (catches previous bug)
  for (int i=0; i!=ncells; ++i) {
    for (int j=0; j!=ncells; ++j) {
      if (i != j) {
        CHECK(cell_list[i] != cell_list[j]);
      }
    }
  }

  Amanzi::AmanziMesh::Entity_ID_List const& face_list = mesh->faces_of_column(0);
  CHECK_EQUAL(ncells+1,face_list.size());

  // check we are not doubling up on some faces (catches previous bug)
  for (int i=0; i!=ncells+1; ++i) {
    for (int j=0; j!=ncells+1; ++j) {
      if (i != j) {
        CHECK(face_list[i] != face_list[j]);
      }
    }
  }

  // construct a column mesh by extracting from mesh
  Amanzi::AmanziMesh::Mesh_MSTK column_mesh(mesh,cell_list, Amanzi::AmanziMesh::CELL,false,Amanzi::getCommSelf());
  
  // Number of cells in column mesh
  int ncells_col = column_mesh.num_entities(Amanzi::AmanziMesh::CELL,
                                            Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(ncells,ncells_col);

  // Number of faces in the column mesh
  int nfaces_col = column_mesh.num_entities(Amanzi::AmanziMesh::FACE,
                                            Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(7*ncells + 1,nfaces_col);
  
  // Check that their parents are as expected
  for (int i = 0; i < ncells_col; ++i) {
    int parent_cell = column_mesh.entity_get_parent(Amanzi::AmanziMesh::CELL,i);
    CHECK_EQUAL(cell_list[i], parent_cell);
  }

  // check we can still get sets
  Amanzi::AmanziMesh::Entity_ID_List set_ids2;
  column_mesh.get_set_entities("organic", Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::ALL, &set_ids2);
  CHECK_EQUAL(10, set_ids2.size());
  
}

