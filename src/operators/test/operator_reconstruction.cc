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

const std::string LIMITERS[3] = {"B-J", "Tensorial", "Kuzmin"};

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
  Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
  Epetra_MultiVector grad_exact(mesh->cell_map(false), 2);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

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
  ReconstructionCell lifting(mesh);
  lifting.Init(field, plist);
  lifting.Compute(); 

  // calculate gradient error
  const Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");
  int ierr = grad_exact.Update(-1.0, grad_computed, 1.0);
  CHECK(!ierr);

  double error[2];
  grad_exact.Norm2(error);
  CHECK_CLOSE(0.0, error[0], 1.0e-12);
  CHECK_CLOSE(0.0, error[1], 1.0e-12);
  
  if (MyPID == 0) printf("errors: %8.4f %8.4f\n", error[0], error[1]);
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
  Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
  Epetra_MultiVector grad_exact(mesh->cell_map(false), 2);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*field)[0][c] = xc[0] + 2 * xc[1];
    if (c < ncells_owned) {
      grad_exact[0][c] = 1.0;
      grad_exact[1][c] = 2.0;
    }
  }

  // create and initialize flux
  // Since limiters do not allow maximum on the outflow bounadry, 
  // we use this trick: re-entering flow everywhere.
  const Epetra_Map& fmap = mesh->face_map(true);
  Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(fmap, 1));

  int dir;
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
  AmanziGeometry::Point velocity(1.0, 2.0), center(0.5, 0.5);
  Amanzi::AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
    velocity = center - xf;
    mesh->face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);
    const AmanziGeometry::Point& normal = mesh->face_normal(f, false, cells[0], &dir);
    (*flux)[0][f] = (velocity * normal) / mesh->face_area(f);
  }

  for (int i = 0; i < 3; i++) {
    std::vector<int> bc_model;
    std::vector<double> bc_value;
    Teuchos::ParameterList plist;
    if (i == 0) {
      plist.set<std::string>("limiter", "Barth-Jespersen");
    } else if (i == 1) {
      plist.set<std::string>("limiter", "tensorial");
    } else if (i == 2) {
      plist.set<std::string>("limiter", "Kuzmin");
    }

    if (i < 2) {
      bc_model.assign(nfaces_wghost, 0);
      bc_value.assign(nfaces_wghost, 0.0);

      for (int f = 0; f < nfaces_wghost; f++) {
        const AmanziGeometry::Point& xf = mesh->face_centroid(f);
        if (fabs(xf[0]) < 1e-6 || fabs(1.0 - xf[0]) < 1e-6 ||
            fabs(xf[1]) < 1e-6 || fabs(1.0 - xf[1]) < 1e-6) {
          bc_model[f] = OPERATOR_BC_DIRICHLET;
          bc_value[f] = xf[0] + 2 * xf[1];
        }
      }
    } else {
      bc_model.assign(nnodes_wghost, 0);
      bc_value.assign(nnodes_wghost, 0.0);
      AmanziGeometry::Point xv(2);

      for (int v = 0; v < nnodes_wghost; v++) {
        mesh->node_get_coordinates(v, &xv);
        if (fabs(xv[0]) < 1e-6 || fabs(1.0 - xv[0]) < 1e-6 ||
            fabs(xv[1]) < 1e-6 || fabs(1.0 - xv[1]) < 1e-6) {
          bc_model[v] = OPERATOR_BC_DIRICHLET;
          bc_value[v] = xv[0] + 2 * xv[1];
        }
      }
    }

    // Compute reconstruction
    ReconstructionCell lifting(mesh);
    lifting.Init(field, plist);
    lifting.Compute(); 

    // Apply limiter
    lifting.InitLimiter(flux);
    lifting.ApplyLimiter(bc_model, bc_value);

    // calculate gradient error
    Epetra_MultiVector grad_computed(*lifting.gradient()->ViewComponent("cell"));
    int ierr = grad_computed.Update(-1.0, grad_exact, 1.0);

    double error[2];
    grad_computed.Norm2(error);
    CHECK_CLOSE(0.0, error[0], 1.0e-12);
    CHECK_CLOSE(0.0, error[1], 1.0e-12);
  
    if (MyPID == 0)
        printf("%9s: errors: %8.4f %8.4f\n", LIMITERS[i].c_str(), error[0], error[1]);
  }
}


TEST(RECONSTRUCTION_SMOOTH_FIELD) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: Accuracy on a smooth field." << std::endl;

  // create rectangular mesh
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);

  for (int n = 14; n < 100; n*=2) { 
    Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, n, n);

    // create and initialize cell-based field ussing f(x,y) = x^2 y + 2 x y^3
    Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
    Epetra_MultiVector grad_exact(mesh->cell_map(false), 2);

    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

    for (int c = 0; c < ncells_wghost; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      double x = xc[0], y = xc[1];
      (*field)[0][c] = x*x*y + 2*x*y*y*y;
      if (c < ncells_owned) {
        grad_exact[0][c] = 2*x*y + 2*y*y*y;
        grad_exact[1][c] = x*x + 6*x*y*y;
      }
    }

    // create and initialize flux
    // Since limiters do not allow maximum on the outflow bounadry, 
    // we use this trick: re-entering flow everywhere.
    const Epetra_Map& fmap = mesh->face_map(true);
    Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(fmap, 1));

    int dir;
    int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
    AmanziGeometry::Point velocity(1.0, 2.0), center(0.5, 0.5);
    Amanzi::AmanziMesh::Entity_ID_List cells;

    for (int f = 0; f < nfaces_wghost; f++) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);
      velocity = center - xf;
      mesh->face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);
      const AmanziGeometry::Point& normal = mesh->face_normal(f, false, cells[0], &dir);
      (*flux)[0][f] = (velocity * normal) / mesh->face_area(f);
    }

    for (int i = 0; i < 3; i++) {
      std::vector<int> bc_model;
      std::vector<double> bc_value;
      Teuchos::ParameterList plist;
      if (i == 0) {
        plist.set<std::string>("limiter", "Barth-Jespersen");
      } else if (i == 1) {
        plist.set<std::string>("limiter", "tensorial");
      } else if (i == 2) {
        plist.set<std::string>("limiter", "Kuzmin");
      }

      if (i < 2) {
        bc_model.assign(nfaces_wghost, 0);
        bc_value.assign(nfaces_wghost, 0.0);

        for (int f = 0; f < nfaces_wghost; f++) {
          const AmanziGeometry::Point& xf = mesh->face_centroid(f);
          double x = xf[0], y = xf[1];
          if (fabs(xf[0]) < 1e-6 || fabs(1.0 - xf[0]) < 1e-6 ||
              fabs(xf[1]) < 1e-6 || fabs(1.0 - xf[1]) < 1e-6) {
            bc_model[f] = OPERATOR_BC_DIRICHLET;
            bc_value[f] = x*x*y + 2*x*y*y*y;
          }
        }
      } else {
        bc_model.assign(nnodes_wghost, 0);
        bc_value.assign(nnodes_wghost, 0.0);
        AmanziGeometry::Point xv(2);

        for (int v = 0; v < nnodes_wghost; v++) {
          mesh->node_get_coordinates(v, &xv);
          double x = xv[0], y = xv[1];
          if (fabs(xv[0]) < 1e-6 || fabs(1.0 - xv[0]) < 1e-6 ||
              fabs(xv[1]) < 1e-6 || fabs(1.0 - xv[1]) < 1e-6) {
            bc_model[v] = OPERATOR_BC_DIRICHLET;
            bc_value[v] = x*x*y + 2*x*y*y*y;
          }
        }
      } 

      // Compute reconstruction
      ReconstructionCell lifting(mesh);
      lifting.Init(field, plist);
      lifting.Compute(); 

      // Apply limiter
      lifting.InitLimiter(flux);
      lifting.ApplyLimiter(bc_model, bc_value);

      // calculate gradient error
      Epetra_MultiVector grad_computed(*lifting.gradient()->ViewComponent("cell"));
      int ierr = grad_computed.Update(-1.0, grad_exact, 1.0);
      CHECK(!ierr);

      double error[2], norms[2];
      grad_computed.Norm2(error);
      grad_exact.Norm2(norms);

      error[0] /= norms[0];
      error[1] /= norms[1];
      CHECK(error[0] < 1.0 / n);
      CHECK(error[1] < 1.0 / n);
  
      if (MyPID == 0)
          printf("%9s: rel errors: %10.6f %10.6f\n", LIMITERS[i].c_str(), error[0], error[1]);
    }
  }
}


