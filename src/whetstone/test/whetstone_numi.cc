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
#include <vector>

#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Point.hh"

#include "NumericalIntegration.hh"
#include "Polynomial.hh"


/* **************************************************************** */
TEST(NUMI_CELL_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Numerical Integartion over 2D cells" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_pentagon.exo", gm, true, true); 
 
  NumericalIntegration numi(mesh);

  int cell(0);
  double val;

  // 0th-order polynomial
  Polynomial poly(2, 0);
  poly(0, 0) = 1.0;
  val = numi.IntegratePolynomialCell(cell, poly);

  printf("order=0  value=%10.6g\n", val);
  CHECK_CLOSE(val, mesh->cell_volume(cell), 1e-10);
 
  // 1st-order polynomial
  poly.Reshape(2, 1);
  poly(1, 0) = 2.0;
  poly(1, 1) = 3.0;
  poly.set_origin(mesh->cell_centroid(cell));
  val = numi.IntegratePolynomialCell(cell, poly);

  printf("order=1  value=%10.6g\n", val);
  CHECK_CLOSE(val, mesh->cell_volume(cell), 1e-10);
 
  delete comm;
}

