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
***************************************************************** */
TEST(RECONSTRUCTION_LINEAR) {
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

  Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 7, 7);

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
  Teuchos::ParameterList plist;
  ReconstructionCell lifting(mesh);
  lifting.Init(field, plist);
  lifting.Compute(); 

  // calculate gradient error
  const Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");
  int ierr = grad_exact.Update(-1.0, grad_computed, 1.0);

  double error[2];
  grad_exact.Norm2(error);
  CHECK_CLOSE(0.0, error[0], 1.0e-12);
  CHECK_CLOSE(0.0, error[1], 1.0e-12);
  
  printf("errors: %8.4f %8.4f\n", error[0], error[1]);
}


TEST(RECONSTRUCTION_LINEAR_LIMITER) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Exactness on linear functions with limiters." << std::endl;

  // create rectangular mesh
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);

  Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 7, 7);

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

  // create and initialize flux
  const Epetra_Map& fmap = mesh->face_map(true);
  Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(fmap, 1));

  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  std::vector<int> bc_model(nfaces_wghost, 0);
  std::vector<double> bc_value(nfaces_wghost, 0.0);

  AmanziGeometry::Point velocity(1.0, 2.0);
  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = (velocity * normal) / mesh->face_area(f);

    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[1]) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = xf[0] + 2 * xf[1];
    }
  }

  // Compute reconstruction
  Teuchos::ParameterList plist;
  plist.set<std::string>("limiter", "Barth-Jespersen");

  ReconstructionCell lifting(mesh);
  lifting.Init(field, plist);
  lifting.Compute(); 

  // Apply limiter
  lifting.InitLimiter(flux);
  lifting.ApplyLimiter(bc_model, bc_value);

  // calculate gradient error
  const Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");
  int ierr = grad_exact.Update(-1.0, grad_computed, 1.0);

  double error[2];
  grad_exact.Norm2(error);
  CHECK_CLOSE(0.0, error[0], 1.0e-12);
  CHECK_CLOSE(0.0, error[1], 1.0e-12);
  
  printf("errors: %8.4f %8.4f\n", error[0], error[1]);
}

