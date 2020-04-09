/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>
#include <fstream>

#include "../Mesh_MSTK.hh"


#include "AmanziMap.hh"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"


// Extract some surfaces as-is from 3D mesh

TEST(Extract_Column_MSTK)
{
  auto comm = Amanzi::getDefaultComm();

  Teuchos::ParameterList reg_spec; // no regions declared here

  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Generate a mesh consisting of 3x3x3 elements
  auto mesh = Teuchos::rcp(
    new Amanzi::AmanziMesh::Mesh_MSTK(0, 0, 0, 1, 1, 1, 3, 3, 3, comm, gm));

  CHECK_EQUAL(1, mesh->build_columns());
  CHECK_EQUAL(9, mesh->num_columns());

  int cell0 = 0;
  int colid = mesh->column_ID(cell0);
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cell_list =
    mesh->cells_of_column(colid);

  CHECK_EQUAL(3, cell_list.extent(0));

  // check we are not doubling up on some cells (catches previous bug)
  for (int i = 0; i != 3; ++i) {
    for (int j = 0; j != 3; ++j) {
      if (i != j) { CHECK(cell_list(i) != cell_list(j)); }
    }
  }


  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> face_list =
    mesh->faces_of_column(colid);
  CHECK_EQUAL(4, face_list.extent(0));

  // check we are not doubling up on some faces (catches previous bug)
  for (int i = 0; i != 4; ++i) {
    for (int j = 0; j != 4; ++j) {
      if (i != j) { CHECK(face_list(i) != face_list(j)); }
    }
  }

  Amanzi::AmanziMesh::Mesh_MSTK column_mesh(
    mesh, cell_list, Amanzi::AmanziMesh::CELL, false, Amanzi::getCommSelf());


  // Number of cells in column mesh

  int ncells_col = column_mesh.num_entities(
    Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(3, ncells_col);


  // Check that their parents are as expected

  for (int i = 0; i < ncells_col; ++i) {
    int parent_cell =
      column_mesh.entity_get_parent(Amanzi::AmanziMesh::CELL, i);
    CHECK_EQUAL(cell_list[i], parent_cell);
  }
}
