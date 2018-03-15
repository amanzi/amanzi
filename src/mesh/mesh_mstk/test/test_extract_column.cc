#include <UnitTest++.h>
#include <iostream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"


// Extract some surfaces as-is from 3D mesh

TEST(Extract_Column_MSTK)
{

  Teuchos::RCP<Epetra_MpiComm> comm_(new Epetra_MpiComm(MPI_COMM_WORLD));

  Teuchos::ParameterList reg_spec; // no regions declared here
  
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, comm_.get()));

  // Generate a mesh consisting of 3x3x3 elements 
  Amanzi::AmanziMesh::Mesh_MSTK mesh(0,0,0,1,1,1,3,3,3,comm_.get(),gm);

  CHECK_EQUAL(1, mesh.build_columns());
  CHECK_EQUAL(9,mesh.num_columns());

  int cell0 = 0;
  int colid = mesh.column_ID(cell0);
  Amanzi::AmanziMesh::Entity_ID_List const& cell_list = mesh.cells_of_column(colid);

  CHECK_EQUAL(3,cell_list.size());

  // check we are not doubling up on some cells (catches previous bug)
  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      if (i != j) {
        CHECK(cell_list[i] != cell_list[j]);
      }
    }
  }


  Amanzi::AmanziMesh::Entity_ID_List const& face_list = mesh.faces_of_column(colid);
  CHECK_EQUAL(4,face_list.size());

  // check we are not doubling up on some faces (catches previous bug)
  for (int i=0; i!=4; ++i) {
    for (int j=0; j!=4; ++j) {
      if (i != j) {
        CHECK(face_list[i] != face_list[j]);
      }
    }
  }
  
  Amanzi::AmanziMesh::Mesh_MSTK column_mesh(mesh,cell_list,
                                            Amanzi::AmanziMesh::CELL,
                                            false,false);


  
  // Number of cells in column mesh

  int ncells_col = column_mesh.num_entities(Amanzi::AmanziMesh::CELL,
                                            Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(3,ncells_col);


  // Check that their parents are as expected

  for (int i = 0; i < ncells_col; ++i) {
    int parent_cell = column_mesh.entity_get_parent(Amanzi::AmanziMesh::CELL,i);
    CHECK_EQUAL(cell_list[i], parent_cell);
  }

}

