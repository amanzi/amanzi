/*
  This is the operators component of the Amanzi code. 

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

#include "UnitTest++.h"

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "OperatorDefs.hh"
#include "ReconstructionCell.hh"


/* *****************************************************************
* Exactness on linear functions.
* **************************************************************** */
TEST(RECONSTRUCTION_LINEAR) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Exactness on linear functions." << std::endl;

  // create rectangular mesh
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);

  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 7, 7);

  // create and initialize cell-based field 
  const Epetra_Map& cmap = mesh->cell_map(false);
  Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(cmap, 1));
  Epetra_MultiVector grad_exact(cmap, 2);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*field)[0][c] = xc[0] + 2 * xc[1];
    grad_exact[0][c] = 1.0;
    grad_exact[1][c] = 2.0;
  }

  // Compute reconstruction
  ReconstructionCell lifting(mesh);
  lifting.Init(field);
  lifting.Compute(); 

  // calculate gradient error
  const Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");
  grad_exact.Update(-1.0, grad_computed, 1.0);

  double error;
  grad_exact.Norm2(&error);
  CHECK_CLOSE(error, 0.0, 1.0e-14);
  
  printf("error=%8.4f\n", error);
}


