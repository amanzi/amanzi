/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Tests for limiters for DG schemes.
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
#include "CompositeVector.hh"
#include "GMVMesh.hh"
#include "DG_Modal.hh"
#include "MeshFactory.hh"

// Amanzi::Operators
#include "AnalyticDG08b.hh"
#include "ErrorAnalysis.hh"
#include "LimiterCell.hh"
#include "OperatorDefs.hh"

const std::string LIMITERS[3] = {"B-J", "B-J c2c", "B-J all"};

/* *****************************************************************
* Limiters must be localized based on a smoothness indicator.
***************************************************************** */
void RunTest(std::string filename)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: Smoothness indicator and limiters for DG" << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Mesh> mesh = meshfactory(filename, Teuchos::null);

  // create and initialize cell-based field 
  int nk(6), dim(2);
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  auto field = Teuchos::rcp(new CompositeVector(cvs1));
  auto field_c = field->ViewComponent("cell", true);

  int order = 2;
  std::string basis("normalized");
  WhetStone::DG_Modal dg(order, mesh, basis);
  AnalyticDG08b ana(mesh, order, true);
  ana.set_shapes(false, true, false);

  ana.InitialGuess(dg, *field_c, 0.0);
  field->ScatterMasterToGhosted("cell");

  int ncells_owned  = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // memory for gradient
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, dim);
  auto grad = Teuchos::rcp(new CompositeVector(cvs2));
  Epetra_MultiVector& grad_c = *grad->ViewComponent("cell");

  for (int i = 0; i < 3; i++) {
    Teuchos::ParameterList plist;
    plist.set<int>("polynomial_order", 2);
    plist.set<bool>("limiter extension for transport", false);

    if (i == 0) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "face to cells");
    } else if (i == 1) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "cell to closest cells");
    } else if (i == 2) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
           .set<std::string>("limiter stencil", "cell to all cells");
    }

    std::vector<int> bc_model(nfaces_wghost, 0);
    std::vector<double> bc_value(nfaces_wghost, 0.0);

    const auto& fmap = mesh->face_map(true);
    const auto& bmap = mesh->exterior_face_map(true);
    for (int bf = 0; bf < bmap.NumMyElements(); ++bf) {
      int f = fmap.LID(bmap.GID(bf));
      const auto& xf = mesh->face_centroid(f);
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.SolutionExact(xf, 0.0);
    }

    // create gradient for limiting
    WhetStone::DenseVector data(nk), data2(nk);

    for (int c = 0; c < ncells_owned; ++c) {
      for (int i = 0; i < nk; ++i) data(i) = (*field_c)[i][c];
      auto poly = dg.cell_basis(c).CalculatePolynomial(mesh, c, order, data);
      auto pgrad = WhetStone::Gradient(poly);
      for (int i = 0; i < dim; ++i) grad_c[i][c] = pgrad[i](0);
    }

    Epetra_MultiVector grad_exact(grad_c);

    // create list of cells where to apply limiter
    double threshold = -4 * std::log10((double) order) - 4.6;
    AmanziMesh::Entity_ID_List ids;

    for (int c = 0; c < ncells_owned; ++c) {
      double honorm(0.0);
      for (int i = dim + 1; i < nk; ++i)
        honorm += (*field_c)[i][c] * (*field_c)[i][c];

      double unorm = honorm;
      for (int i = 0; i <= dim; ++i)
        unorm += (*field_c)[i][c] * (*field_c)[i][c];

      if (unorm > 0.0 && std::log10(honorm / unorm) > threshold) {
        ids.push_back(c);
      }
    }

    // Apply limiter
    LimiterCell limiter(mesh);
    limiter.Init(plist);
    limiter.ApplyLimiter(ids, field_c, 0, grad, bc_model, bc_value);

    // calculate gradient error
    double err_int, err_glb, gnorm;
    ComputeGradError(mesh, grad_c, grad_exact, err_int, err_glb, gnorm);
    // CHECK_CLOSE(0.0, err_int + err_glb, 1.0e-12);

    int nids, itmp = ids.size();
    mesh->get_comm()->SumAll(&itmp, &nids, 1);
    double fraction = 100.0 * nids / grad_c.GlobalLength();
    if (MyPID == 0) 
      printf("%9s: errors: %10.6f %10.6f  ||grad||=%8.4f  indicator=%5.1f%%\n",
          LIMITERS[i].c_str(), err_int, err_glb, gnorm, fraction);
    CHECK(fraction < 15.0);

    // calculate error of limited gradient
    double err(0.0), unorm(0.0);
    for (int n = 0; n < ids.size(); ++n) {
      int c = ids[n];
      double volume = mesh->cell_volume(c);

      for (int i = 0; i < nk; ++i) data2(i) = (*field_c)[i][c];

      data(0) = (*field_c)[0][c];
      for (int i = 0; i < dim; ++i) data(i + 1) = grad_c[i][c];
      for (int i = dim + 1; i < nk; ++i) data(i) = 0.0;
      dg.cell_basis(c).ChangeBasisNaturalToMy(data);

      data2 -= data;
      err += (data2 * data2) * volume;
      unorm += (data * data) * volume;
    }

    double tmp = err;
    mesh->get_comm()->SumAll(&tmp, &err, 1);
    tmp = unorm;
    mesh->get_comm()->SumAll(&tmp, &unorm, 1);
    if (MyPID == 0) 
      printf("       sol errors: %10.6f   ||u||=%8.4f\n", std::pow(err, 0.5), std::pow(unorm, 0.5));
    CHECK(err < 0.02);
  }
}


TEST(LIMITER_SMOOTHNESS_INDICATOR_2D) {
  RunTest("test/circle_quad10.exo");
  // RunTest("test/circle_quad20.exo");
  // RunTest("test/circle_quad40.exo");
  // RunTest("test/circle_quad80.exo");
  // RunTest("test/circle_quad160.exo");
}

