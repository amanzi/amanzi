/*
  WhetStone

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "DenseMatrix.hh"
#include "dg.hh"


/* ****************************************************************
* Test of DG mass matrices
**************************************************************** */
TEST(DG_MASS_MATRIX) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: DG mass matrices" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1); 
 
  DG dg(mesh);

  for (int k = 0; k < 3; k++) {
    int nk = (k + 1) * (k + 2) / 2;
    DenseMatrix M(nk, nk);

    dg.TaylorMassMatrix(0, k, M);
    std::cout << M << std::endl;

    CHECK_CLOSE(M(0, 0), 1.0, 1e-12);
    if (k > 0) {
      CHECK_CLOSE(M(1, 1), 1.0 / 12, 1e-12);
    }
    if (k > 1) {
      CHECK_CLOSE(M(3, 3), 1.0 / 80, 1e-12);
      CHECK_CLOSE(M(4, 4), 1.0 / 144, 1e-12);
    }
  }
}


