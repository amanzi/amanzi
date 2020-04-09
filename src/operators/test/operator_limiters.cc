/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Tests for limiters.
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
#include "LimiterCell.hh"
#include "OperatorDefs.hh"
#include "ReconstructionCell.hh"

const std::string LIMITERS[7] = {"B-J", "Tensorial", "Tens. c2c", "Kuzmin", "B-J c2c", "B-J all", "M-G all"};

/* *****************************************************************
* Limiters must be 1 on linear functions in two dimensions
***************************************************************** */
TEST(LIMITER_LINEAR_FUNCTION_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: Limiters for linear functions in 2D." << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 7, 7);

  // create and initialize cell-based field 
  Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
  Epetra_MultiVector grad_exact(mesh->cell_map(false), 2);

  int ncells_owned  = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*field)[0][c] = xc[0] + 2 * xc[1];
    if (c < ncells_owned) {
      grad_exact[0][c] = 1.0;
      grad_exact[1][c] = 2.0;
    }
  }

  for (int i = 0; i < 7; i++) {
    std::vector<int> bc_model;
    std::vector<double> bc_value;
    Teuchos::ParameterList plist;
    plist.set<int>("polynomial_order", 1);
    plist.set<bool>("limiter extension for transport", false);

    if (i == 0) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "face to cells");
    } else if (i == 1) {
      plist.set<std::string>("limiter", "tensorial");
    } else if (i == 2) {
      plist.set<std::string>("limiter", "tensorial")
           .set<std::string>("limiter stencil", "cell to closest cells");
    } else if (i == 3) {
      plist.set<std::string>("limiter", "Kuzmin");
    } else if (i == 4) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "cell to closest cells");
    } else if (i == 5) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "cell to all cells");
    } else if (i == 6) {
      plist.set<std::string>("limiter", "Michalak-Gooch")
           .set<std::string>("limiter stencil", "cell to all cells");
    }

    if (i != 3) {
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
    lifting.ComputeGradient(); 

    // Apply limiter
    LimiterCell limiter(mesh);
    limiter.Init(plist);
    limiter.ApplyLimiter(field, 0, lifting.gradient(), bc_model, bc_value);

    // calculate gradient error.
    double err_int, err_glb, gnorm;
    Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");

    ComputeGradError(mesh, grad_computed, grad_exact, err_int, err_glb, gnorm);
    // Michalak-Gooch limiter is not linearity preserving near boundary
    CHECK_CLOSE(0.0, err_int, 1.0e-12);
    if (i < 6) CHECK_CLOSE(0.0, err_glb, 1.0e-12);

    if (MyPID == 0)
        printf("%9s: errors: %8.4f %8.4f\n", LIMITERS[i].c_str(), err_int, err_glb);
  }
}


/* *****************************************************************
* Limiters must be 1 on linear functions in three dimensions.
***************************************************************** */
TEST(LIMITER_LINEAR_FUNCTION_3D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: Limiters for linear functions in 3D." << std::endl;

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

  // create and initialize flux
  // Since limiters do not allow maximum on the outflow bounadry, 
  // we use this trick: re-entering flow everywhere.
  const Epetra_Map& fmap = mesh->face_map(true);
  Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(fmap, 1));

  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  AmanziGeometry::Point velocity(3), center(0.5, 0.5, 0.5);

  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
    velocity = center - xf;
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = (velocity * normal) / mesh->face_area(f);
  }

  for (int i = 0; i < 4; i++) {
    std::vector<int> bc_model;
    std::vector<double> bc_value;
    Teuchos::ParameterList plist;
    plist.set<int>("polynomial_order", 1);
    plist.set<bool>("limiter extension for transport", false);

    if (i == 0) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "face to cells");
    } else if (i == 1) {
      plist.set<std::string>("limiter", "tensorial");
    } else if (i == 2) {
      plist.set<std::string>("limiter", "tensorial")
           .set<std::string>("limiter stencil", "cell to closest cells");
    } else if (i == 3) {
      plist.set<std::string>("limiter", "Kuzmin");
    }

    if (i != 3) {
      bc_model.assign(nfaces_wghost, 0);
      bc_value.assign(nfaces_wghost, 0.0);

      for (int f = 0; f < nfaces_wghost; f++) {
        const AmanziGeometry::Point& xf = mesh->face_centroid(f);
        if (fabs(xf[0]) < 1e-6 || fabs(1.0 - xf[0]) < 1e-6 ||
            fabs(xf[1]) < 1e-6 || fabs(1.0 - xf[1]) < 1e-6 ||
            fabs(xf[2]) < 1e-6 || fabs(1.0 - xf[2]) < 1e-6) {
          bc_model[f] = OPERATOR_BC_DIRICHLET;
          bc_value[f] = xf[0] + 2 * xf[1] + 3 * xf[2];
        }
      }
    } else {
      bc_model.assign(nnodes_wghost, 0);
      bc_value.assign(nnodes_wghost, 0.0);
      AmanziGeometry::Point xv(3);

      for (int v = 0; v < nnodes_wghost; v++) {
        mesh->node_get_coordinates(v, &xv);
        if (fabs(xv[0]) < 1e-6 || fabs(1.0 - xv[0]) < 1e-6 ||
            fabs(xv[1]) < 1e-6 || fabs(1.0 - xv[1]) < 1e-6 ||
            fabs(xv[2]) < 1e-6 || fabs(1.0 - xv[2]) < 1e-6) {
          bc_model[v] = OPERATOR_BC_DIRICHLET;
          bc_value[v] = xv[0] + 2 * xv[1] + 3 * xv[2];
        }
      }
    }

    // Compute reconstruction
    ReconstructionCell lifting(mesh);
    lifting.Init(field, plist);
    lifting.ComputeGradient(); 

    // Apply limiter
    LimiterCell limiter(mesh);
    limiter.Init(plist, flux);
    limiter.ApplyLimiter(field, 0, lifting.gradient(), bc_model, bc_value);

    // calculate gradient error
    double err_int, err_glb, gnorm;
    Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");

    ComputeGradError(mesh, grad_computed, grad_exact, err_int, err_glb, gnorm);
    CHECK_CLOSE(0.0, err_int + err_glb, 1.0e-12);

    if (MyPID == 0)
        printf("%9s: errors: %8.4f %8.4f\n", LIMITERS[i].c_str(), err_int, err_glb);
  }
}


/* *****************************************************************
* Convergence of limited functions in two dimensions.
***************************************************************** */
TEST(LIMITER_SMOOTH_FIELD_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: Accuracy on a smooth field in 2D." << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  for (int n = 14; n < 100; n*=2) { 
    Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, n, n - 1);

    // create and initialize cell-based field ussing f(x,y) = x^2 y + 2 x y^3
    Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
    Epetra_MultiVector grad_exact(mesh->cell_map(false), 2);

    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

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

    int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
    AmanziGeometry::Point velocity(1.0, 2.0), center(0.5, 0.5);

    for (int f = 0; f < nfaces_wghost; f++) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);
      velocity = center - xf;
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      (*flux)[0][f] = (velocity * normal) / mesh->face_area(f);
    }

    for (int i = 0; i < 7; i++) {
      std::vector<int> bc_model;
      std::vector<double> bc_value;
      Teuchos::ParameterList plist;
      plist.set<int>("polynomial_order", 1);
      plist.set<bool>("limiter extension for transport", false);

      if (i == 0) {
        plist.set<std::string>("limiter", "Barth-Jespersen")
             .set<std::string>("limiter stencil", "face to cells");
      } else if (i == 1) {
        plist.set<std::string>("limiter", "tensorial");
      } else if (i == 2) {
        plist.set<std::string>("limiter", "tensorial")
             .set<std::string>("limiter stencil", "cell to closest cells");
      } else if (i == 3) {
        plist.set<std::string>("limiter", "Kuzmin");
      } else if (i == 4) {
        plist.set<std::string>("limiter", "Barth-Jespersen")
             .set<std::string>("limiter stencil", "cell to closest cells");
      } else if (i == 5) {
        plist.set<std::string>("limiter", "Barth-Jespersen")
             .set<std::string>("limiter stencil", "cell to all cells");
      } else if (i == 6) {
        plist.set<std::string>("limiter", "Michalak-Gooch")
             .set<std::string>("limiter stencil", "cell to all cells");
      }

      if (i != 3) {
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
      lifting.ComputeGradient(); 

      // Apply limiter
      LimiterCell limiter(mesh);
      limiter.Init(plist, flux);
      limiter.ApplyLimiter(field, 0, lifting.gradient(), bc_model, bc_value);

      // calculate gradient error
      double err_int, err_glb, gnorm;
      Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");

      ComputeGradError(mesh, grad_computed, grad_exact, err_int, err_glb, gnorm);

      if (MyPID == 0)
          printf("%9s: rel errors: %9.5f %9.5f\n", LIMITERS[i].c_str(), err_int, err_glb);

      CHECK(err_int + err_glb < 1.0 / n);
    }
  }
}


/* *****************************************************************
* Convergence of limited functions in three dimensions.
***************************************************************** */
TEST(LIMITER_SMOOTH_FIELD_3D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: Accuracy on a smooth field in 3D" << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  for (int n = 14; n < 50; n*=2) { 
    Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n - 2, n - 1);

    // create and initialize cell-based field f(x,y,z) = x^2 y z^2 + 2 x y^3 z
    Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
    Epetra_MultiVector grad_exact(mesh->cell_map(false), 3);

    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

    for (int c = 0; c < ncells_wghost; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      double x = xc[0], y = xc[1], z = xc[2];
      (*field)[0][c] = x*x*y*z*z + 2*x*y*y*y*z;
      if (c < ncells_owned) {
        grad_exact[0][c] = 2*x*y*z*z + 2*y*y*y*z;
        grad_exact[1][c] = x*x*z*z + 6*x*y*y*z;
        grad_exact[2][c] = 2*x*x*y*z + 2*x*y*y*y;
      }
    }

    // create and initialize flux
    // Since limiters do not allow maximum on the outflow bounadry, 
    // we use this trick: re-entering flow everywhere.
    const Epetra_Map& fmap = mesh->face_map(true);
    Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(fmap, 1));

    int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
    AmanziGeometry::Point velocity(3), center(0.5, 0.5, 0.5);

    for (int f = 0; f < nfaces_wghost; f++) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);
      velocity = center - xf;
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      (*flux)[0][f] = (velocity * normal) / mesh->face_area(f);
    }

    for (int i = 0; i < 3; i++) {
      std::vector<int> bc_model;
      std::vector<double> bc_value;
      Teuchos::ParameterList plist;
      plist.set<int>("polynomial_order", 1);
      plist.set<bool>("limiter extension for transport", false);

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
          double x = xf[0], y = xf[1], z = xf[2];
          if (fabs(xf[0]) < 1e-6 || fabs(1.0 - xf[0]) < 1e-6 ||
              fabs(xf[1]) < 1e-6 || fabs(1.0 - xf[1]) < 1e-6 ||
              fabs(xf[2]) < 1e-6 || fabs(1.0 - xf[2]) < 1e-6) {
            bc_model[f] = OPERATOR_BC_DIRICHLET;
            bc_value[f] = x*x*y*z*z + 2*x*y*y*y*z;
          }
        }
      } else {
        bc_model.assign(nnodes_wghost, 0);
        bc_value.assign(nnodes_wghost, 0.0);
        AmanziGeometry::Point xv(3);

        for (int v = 0; v < nnodes_wghost; v++) {
          mesh->node_get_coordinates(v, &xv);
          double x = xv[0], y = xv[1], z = xv[2];
          if (fabs(xv[0]) < 1e-6 || fabs(1.0 - xv[0]) < 1e-6 ||
              fabs(xv[1]) < 1e-6 || fabs(1.0 - xv[1]) < 1e-6 ||
              fabs(xv[2]) < 1e-6 || fabs(1.0 - xv[2]) < 1e-6) {
            bc_model[v] = OPERATOR_BC_DIRICHLET;
            bc_value[v] = x*x*y*z*z + 2*x*y*y*y*z;
          }
        }
      } 

      // Compute reconstruction
      ReconstructionCell lifting(mesh);
      lifting.Init(field, plist);
      lifting.ComputeGradient(); 

      // Apply limiter
      LimiterCell limiter(mesh);
      limiter.Init(plist, flux);
      limiter.ApplyLimiter(field, 0, lifting.gradient(), bc_model, bc_value);

      // calculate gradient error
      double err_int, err_glb, gnorm;
      Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");

      ComputeGradError(mesh, grad_computed, grad_exact, err_int, err_glb, gnorm);

      if (MyPID == 0)
          printf("n=%d  %9s: rel errors: %9.5f %9.5f\n", n, LIMITERS[i].c_str(), err_int, err_glb);

      CHECK(err_int + err_glb < 1.0 / n);
    }
  }
}


/* *****************************************************************
* Convergece of limited functions in two dimensions.
***************************************************************** */
void SmoothField2DPoly(double extension)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: smooth field on a polygonal mesh, extension=" << extension << std::endl;

  // load polygonal mesh
  Teuchos::ParameterList region_list;
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // create and initialize cell-based field ussing f(x,y) = x^2 y + 2 x y^3
  Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
  Epetra_MultiVector grad_exact(mesh->cell_map(false), 2);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double x = xc[0], y = xc[1];
    (*field)[0][c] = x*x*y + 2*x*y*y*y;
    if (c < ncells_owned) {
      grad_exact[0][c] = 2*x*y + 2*y*y*y;
      grad_exact[1][c] = x*x + 6*x*y*y;
    }
  }

  // Create and initialize flux
  // Since limiters do not allow maximum on the outflow bounadry, 
  // we use this trick: re-entering flow everywhere.
  const Epetra_Map& fmap = mesh->face_map(true);
  Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(fmap, 1));

  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  AmanziGeometry::Point velocity(1.0, 2.0), center(0.5, 0.5);

  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
    velocity = center - xf;
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = (velocity * normal) / mesh->face_area(f);
  }

  for (int i = 0; i < 7; i++) {
    std::vector<int> bc_model;
    std::vector<double> bc_value;
    Teuchos::ParameterList plist;
    plist.set<int>("polynomial_order", 1);
    plist.set<bool>("limiter extension for transport", extension);

    if (i == 0) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "face to cells");
    } else if (i == 1) {
      plist.set<std::string>("limiter", "tensorial");
    } else if (i == 2) {
      plist.set<std::string>("limiter", "tensorial")
           .set<std::string>("limiter stencil", "cell to closest cells");
    } else if (i == 3) {
      plist.set<std::string>("limiter", "Kuzmin");
    } else if (i == 4) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "cell to closest cells");
    } else if (i == 5) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "cell to all cells");
    } else if (i == 6) {
      plist.set<std::string>("limiter", "Michalak-Gooch")
           .set<std::string>("limiter stencil", "cell to all cells");
    }

    if (i != 3) {
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
    lifting.ComputeGradient(); 

    // Apply limiter
    LimiterCell limiter(mesh);
    limiter.Init(plist, flux);
    limiter.ApplyLimiter(field, 0, lifting.gradient(), bc_model, bc_value);

    // calculate gradient error
    double err_int, err_glb, gnorm;
    Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");

    ComputeGradError(mesh, grad_computed, grad_exact, err_int, err_glb, gnorm);

    if (MyPID == 0)
        printf("%9s: rel errors: %9.5f %9.5f\n", LIMITERS[i].c_str(), err_int, err_glb);
  }
}

TEST(LIMITER_SMOOTH_FIELD_POLYMESH) {
  SmoothField2DPoly(false);
  SmoothField2DPoly(true);
}


/* *****************************************************************
* Limiters must be 1 on linear functions in three dimensions.
***************************************************************** */
TEST(LIMITER_LINEAR_FUNCTION_FRACTURES) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: Limiters for linear functions in fractures." << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/fractures.exo");

  // create and initialize cell-based field 
  Teuchos::RCP<Epetra_MultiVector> field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));
  Epetra_MultiVector grad_exact(mesh->cell_map(false), 3);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    if (fabs(xc[2] - 0.5) < 1e-10) {
      (*field)[0][c] = xc[0] + 2 * xc[1];
      if (c < ncells_owned) {
        grad_exact[0][c] = 1.0;
        grad_exact[1][c] = 2.0;
        grad_exact[2][c] = 0.0;
      }
    } else if (fabs(xc[1] - 0.5) < 1e-10) {
      (*field)[0][c] = xc[0] + 3 * xc[2];
      if (c < ncells_owned) {
        grad_exact[0][c] = 1.0;
        grad_exact[1][c] = 0.0;
        grad_exact[2][c] = 3.0;
      }
    }
  }

  // create and initialize flux
  // Since limiters do not allow maximum on the outflow bounadry, 
  // we use this trick: re-entering flow everywhere.
  const Epetra_Map& fmap = mesh->face_map(true);
  Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(fmap, 1));

  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  AmanziGeometry::Point velocity(3), center(0.5, 0.5, 0.5);

  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
    velocity = center - xf;
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = (velocity * normal) / mesh->face_area(f);
  }

  for (int i = 1; i < 2; i++) {
    std::vector<int> bc_model;
    std::vector<double> bc_value;
    Teuchos::ParameterList plist;
    plist.set<int>("polynomial_order", 1);
    plist.set<bool>("limiter extension for transport", false);

    if (i == 1) {
      plist.set<std::string>("limiter", "tensorial");
    }

    if (i < 2) {
      bc_model.assign(nfaces_wghost, 0);
      bc_value.assign(nfaces_wghost, 0.0);

      for (int f = 0; f < nfaces_wghost; f++) {
        const AmanziGeometry::Point& xf = mesh->face_centroid(f);
        if (fabs(xf[0]) < 1e-6 || fabs(1.0 - xf[0]) < 1e-6 ||
            fabs(xf[1]) < 1e-6 || fabs(1.0 - xf[1]) < 1e-6 ||
            fabs(xf[2]) < 1e-6 || fabs(1.0 - xf[2]) < 1e-6) {
          bc_model[f] = OPERATOR_BC_DIRICHLET;
          bc_value[f] = xf[0] + 2 * xf[1] + 3 * xf[2];
        }
      }
    }

    // Compute reconstruction
    ReconstructionCell lifting(mesh);
    lifting.Init(field, plist);
    lifting.ComputeGradient(); 

    // calculate gradient error
    double err_int, err_glb, gnorm;
    Epetra_MultiVector& grad_computed = *lifting.gradient()->ViewComponent("cell");

    ComputeGradError(mesh, grad_computed, grad_exact, err_int, err_glb, gnorm);
    CHECK_CLOSE(0.0, err_int + err_glb, 1.0e-12);

    if (MyPID == 0)
        printf("%9s: errors: %8.4f %8.4f\n", LIMITERS[i].c_str(), err_int, err_glb);
  }
}



