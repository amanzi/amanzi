/*
  Operators

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

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"

// Amanzi::Operators
#include "ErrorAnalysis.hh"
#include "OperatorDefs.hh"
#include "ReconstructionCellGrad.hh"
#include "ReconstructionCellPoly.hh"
#include "SmoothnessIndicatorShu.hh"


/* *****************************************************************
* Exactness on linear functions in two dimensions
***************************************************************** */
TEST(RECONSTRUCTION_LINEAR_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Exactness on linear functions in 2D." << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 7, 7);

  // create and initialize cell-based field 
  Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
  Epetra_MultiVector grad_exact(mesh->cell_map(false), 2);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*field)[0][c] = xc[0] + 2 * xc[1];
    if (c < ncells_owned) {
      grad_exact[0][c] = 1.0;
      grad_exact[1][c] = 2.0;
    }
  }

  // Compute reconstruction
  Teuchos::ParameterList plist;
  plist.set<std::string>("limiter", "tensorial");
  plist.set<int>("polynomial_order", 1);
  plist.set<bool>("limiter extension for transport", false);

  ReconstructionCellGrad lifting(mesh);
  lifting.Init(plist);
  lifting.Compute(field); 

  // calculate gradient error
  double err_int, err_glb, gnorm;
  auto& grad_computed = *lifting.data()->ViewComponent("cell");

  ComputePolyError(mesh, grad_computed, grad_exact, err_int, err_glb, gnorm);
  CHECK_CLOSE(0.0, err_int + err_glb, 1.0e-12);

  if (MyPID == 0) printf("errors (interior & global): %8.4f %8.4f\n", err_int, err_glb);
}


/* *****************************************************************
* Exactness on linear functions in three dimensions
***************************************************************** */
TEST(RECONSTRUCTION_LINEAR_3D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Exactness on linear functions in 3D." << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 7, 6, 5);

  // create and initialize cell-based field 
  Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
  Epetra_MultiVector grad_exact(mesh->cell_map(false), 3);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*field)[0][c] = xc[0] + 2 * xc[1] + 3 * xc[2];
    if (c < ncells_owned) {
      grad_exact[0][c] = 1.0;
      grad_exact[1][c] = 2.0;
      grad_exact[2][c] = 3.0;
    }
  }

  // Compute reconstruction
  Teuchos::ParameterList plist;
  plist.set<std::string>("limiter", "tensorial");
  plist.set<int>("polynomial_order", 1);
  plist.set<bool>("limiter extension for transport", false);

  ReconstructionCellGrad lifting(mesh);
  lifting.Init(plist);
  lifting.Compute(field);

  // calculate gradient error
  double err_int, err_glb, gnorm;
  auto& grad_computed = *lifting.data()->ViewComponent("cell");

  ComputePolyError(mesh, grad_computed, grad_exact, err_int, err_glb, gnorm);
  CHECK_CLOSE(0.0, err_int + err_glb, 1.0e-12);

  if (MyPID == 0) printf("errors (interior & global): %8.4f %8.4f\n", err_int, err_glb);
}


/* *****************************************************************
* Exactness on quadratic functions in two dimensions
***************************************************************** */
TEST(RECONSTRUCTION_QUADRATIC_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Exactness on quadratic functions in 2D." << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 7, 7);

  // create and initialize cell-based field 
  Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
  Epetra_MultiVector poly_exact(mesh->cell_map(false), 5);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // -- cell values
  double h(1.0 / 7);
  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*field)[0][c] = xc[0] + 2 * xc[1] 
                   + 3 * xc[0] * xc[0] + 4 * xc[0] * xc[1] + 5 * xc[1] * xc[1]
                   + 2.0 * h * h / 3;
    if (c < ncells_owned) {
      poly_exact[0][c] = 1.0 + 6 * xc[0] + 4 * xc[1];
      poly_exact[1][c] = 2.0 + 10 * xc[1] + 4 * xc[0];
      poly_exact[2][c] = 3.0;
      poly_exact[3][c] = 4.0;
      poly_exact[4][c] = 5.0;
    }
  }

  // -- boundary values
  AmanziMesh::Entity_ID_List cells;

  auto bcs = Teuchos::rcp(new Operators::BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR)); 
  std::vector<int>& bc_model = bcs->bc_model();
  std::vector<double>& bc_value = bcs->bc_value();

  for (int f = 0; f < nfaces_wghost; ++f) {
    mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    if (cells.size() == 1) {
      const auto& xf = mesh->face_centroid(f);
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = xf[0] + 2 * xf[1] + 3 * xf[0] * xf[0] + 4 * xf[0] * xf[1] + 5 * xf[1] * xf[1];
    }
  }

  // Compute reconstruction
  Teuchos::ParameterList plist;
  plist.set<std::string>("limiter", "tensorial");
  plist.set<int>("polynomial_order", 1);
  plist.set<bool>("limiter extension for transport", false);

  auto lifting = Teuchos::rcp(new ReconstructionCellPoly(mesh));
  lifting->Init(plist);
  lifting->Compute(field, 0, bcs);

  // -- calculate polynomial error
  double err_int, err_glb, gnorm;
  auto& poly_computed = *lifting->data()->ViewComponent("cell");

  ComputePolyError(mesh, poly_computed, poly_exact, err_int, err_glb, gnorm);
  CHECK_CLOSE(0.0, err_int + err_glb, 1.0e-12);

  if (MyPID == 0) printf("errors (interior & global): %8.4f %8.4f\n", err_int, err_glb);

  // Analyze smoothness of the reconstruction
  lifting->data()->ScatterMasterToGhosted("cell");

  SmoothnessIndicatorShu indicator(mesh);
  indicator.Compute(lifting);
  auto measure = indicator.get_measure();

  for (int c = 0; c < measure.size(); ++c) {
    CHECK_CLOSE(measure[c], 1.0, 1e-12);
  }
}

