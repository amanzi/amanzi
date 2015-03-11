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

  bool request_faces=true, request_edges=true;

  RCP<Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, NULL, 
			       request_faces, request_edges); 
  // RCP<Mesh> mesh = meshfactory("test/one_cell.exo"); 
 
  MFD3D_Electromagnetics mfd(mesh);

  int nrows = 12, nedges = 12, cell = 0;
  Tensor T(3, 1);
  T(0, 0) = 1;

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
  AmanziMesh::Entity_ID_List edges;
  mesh->cell_get_edges(cell, &edges);
    
  int d = mesh->space_dimension();
  Point p(d);

  double xi, yi, xj;
  double vxx = 0.0, vxy = 0.0, volume = mesh->cell_volume(cell); 
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh->edge_vector(e);
    xi = tau[0];
    yi = tau[1];
    for (int j = 0; j < nedges; j++) {
      e = edges[j];
      const AmanziGeometry::Point& tau = mesh->edge_vector(e);
      xj = p[0];
      vxx += M(i, j) * xi * xj;
      vxy += M(i, j) * yi * xj;
    }
  }
  // CHECK_CLOSE(vxx, volume, 1e-10);
  // CHECK_CLOSE(vxy, 0.0, 1e-10);

  delete comm;
}


