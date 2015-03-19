/*
  This is the mimetic discretization component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"
#include "UnitTest++.h"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "mfd3d_electromagnetics.hh"
#include "tensor.hh"


/* **************************************************************** */
TEST(MASS_MATRIX_3D) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Mass matrix for edge elements in 3D" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);

  bool request_faces(true), request_edges(true);

  // RCP<Mesh> mesh = meshfactory("test/dodecahedron.exo", NULL, request_faces, request_edges); 
  RCP<Mesh> mesh = meshfactory("test/one_cell.exo", NULL, request_faces, request_edges); 
  // RCP<Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 3, NULL, true, true); 
 
  MFD3D_Electromagnetics mfd(mesh);

  int cell = 0;
  AmanziMesh::Entity_ID_List edges;
  mesh->cell_get_edges(cell, &edges);

  int nedges = edges.size();
  int nrows = nedges;

  Tensor T(3, 2);
  T(0, 0) = 2.0;
  T(1, 1) = 1.0;
  T(0, 1) = 1.0;
  T(1, 0) = 1.0;
  T(2, 2) = 1.0;

  DenseMatrix M(nrows, nrows);
  mfd.MassMatrix(cell, T, M);

  printf("Stiffness matrix for cell %3d\n", cell);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nrows; j++ ) printf("%8.4f ", M(i, j)); 
    printf("\n");
  }

  // verify SPD propery
  for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

  // verify exact integration property
  double xi, yi, xj;
  double vxx = 0.0, vxy = 0.0, volume = mesh->cell_volume(cell); 
  for (int i = 0; i < nedges; i++) {
    int e1 = edges[i];
    const AmanziGeometry::Point& t1 = mesh->edge_vector(e1);
    double a1 = mesh->edge_length(e1);

    xi = t1[0] / a1;
    yi = t1[1] / a1;
    for (int j = 0; j < nedges; j++) {
      int e2 = edges[j];
      const AmanziGeometry::Point& t2 = mesh->edge_vector(e2);
      double a2 = mesh->edge_length(e2);
      xj = t2[0] / a2;
      vxx += M(i, j) * xi * xj;
      vxy += M(i, j) * yi * xj;
    }
  }
  CHECK_CLOSE(volume, vxx, 1e-10);
  CHECK_CLOSE(-volume, vxy, 1e-10);

  delete comm;
}


