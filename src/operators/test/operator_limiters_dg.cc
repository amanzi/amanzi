/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

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
#include "LimiterCellDG.hh"
#include "OperatorDefs.hh"
#include "ReconstructionCellLinear.hh"

const std::string LIMITERS[3] = { "B-J", "B-J c2c", "B-J all" };

/* *****************************************************************
* Limiters must be localized based on a smoothness indicator.
***************************************************************** */
void
RunTest(std::string filename, std::string basis, double& l2norm)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0)
    std::cout << "\nTest: Smoothness indicator and limiters for DG, basis=" << basis << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(filename);

  // create and initialize cell-based field
  int nk(6), dim(2);
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, nk);
  auto field = Teuchos::rcp(new CompositeVector(cvs1));
  auto field_c = field->ViewComponent("cell", true);

  int order = 2;
  Teuchos::ParameterList dglist;
  dglist.set<int>("method order", order).set<std::string>("dg basis", basis);

  WhetStone::DG_Modal dg(dglist, mesh);
  AnalyticDG08b ana(mesh, order, true);
  ana.set_shapes(false, true, false);

  ana.InitialGuess(dg, *field_c, 0.0);
  field->ScatterMasterToGhosted("cell");

  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  // memory for gradient
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, dim);
  auto grad = Teuchos::rcp(new CompositeVector(cvs2));
  Epetra_MultiVector& grad_c = *grad->ViewComponent("cell");

  for (int tst = 0; tst < 3; ++tst) {
    Teuchos::ParameterList plist;
    plist.set<int>("polynomial_order", 2);
    plist.set<bool>("limiter extension for transport", false);

    if (tst == 0) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
        .set<std::string>("limiter stencil", "face to cells");
    } else if (tst == 1) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
        .set<std::string>("limiter stencil", "cell to closest cells");
    } else if (tst == 2) {
      plist.set<std::string>("limiter", "Barth-Jespersen")
        .set<std::string>("limiter stencil", "cell to all cells");
    }

    std::vector<int> bc_model(nfaces_wghost, 0);
    std::vector<double> bc_value(nfaces_wghost, 0.0);

    const auto& fmap = mesh->getMap(AmanziMesh::Entity_kind::FACE,true);
    const auto& bmap = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
    for (int bf = 0; bf < bmap.NumMyElements(); ++bf) {
      int f = fmap.LID(bmap.GID(bf));
      const auto& xf = mesh->getFaceCentroid(f);
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.SolutionExact(xf, 0.0);
    }

    // create gradient for limiting
    WhetStone::DenseVector data(nk), data2(nk), data3(nk);

    for (int c = 0; c < ncells_owned; ++c) {
      for (int i = 0; i < nk; ++i) data(i) = (*field_c)[i][c];
      auto poly = dg.cell_basis(c).CalculatePolynomial(mesh, c, order, data);
      auto pgrad = WhetStone::Gradient(poly);
      for (int i = 0; i < dim; ++i) grad_c[i][c] = pgrad[i](0);
    }

    Epetra_MultiVector grad_exact(grad_c);

    // create list of cells where to apply limiter
    double L(0.0);
    double threshold = -4 * std::log10((double)order) - L;
    AmanziMesh::Entity_ID_View ids("ids", ncells_owned);
    std::size_t ids_count = 0; 

    for (int c = 0; c < ncells_owned; ++c) {
      double honorm(0.0);
      for (int i = dim + 1; i < nk; ++i) honorm += (*field_c)[i][c] * (*field_c)[i][c];

      double unorm = honorm;
      for (int i = 0; i <= dim; ++i) unorm += (*field_c)[i][c] * (*field_c)[i][c];

      if (unorm > 0.0 && std::log10(honorm / unorm) > threshold) { ids[ids_count++] = c; }
    }
    Kokkos::resize(ids,ids_count); 

    // Apply limiter
    auto lifting = Teuchos::rcp(new ReconstructionCellLinear(mesh, grad));
    LimiterCell limiter(mesh);
    limiter.Init(plist);
    limiter.ApplyLimiter(ids, field_c, 0, lifting, bc_model, bc_value);

    // calculate gradient error
    double err_int, err_glb, gnorm;
    ComputePolyError(mesh, grad_c, grad_exact, err_int, err_glb, gnorm);
    // CHECK_CLOSE(0.0, err_int + err_glb, 1.0e-12);

    int nids, itmp = ids.size();
    mesh->getComm()->SumAll(&itmp, &nids, 1);
    double fraction = 100.0 * nids / grad_c.GlobalLength();
    if (MyPID == 0)
      printf("%9s: errors: %10.6f %10.6f  ||grad||=%8.4f  indicator=%5.1f%%\n",
             LIMITERS[tst].c_str(),
             err_int,
             err_glb,
             gnorm,
             fraction);
    CHECK(fraction < 15.0);

    // calculate true L2 error of limited gradient
    WhetStone::Tensor K(dim, 1);
    K(0, 0) = 1.0;

    double err(0.0);
    l2norm = 0.0;
    for (int n = 0; n < ids.size(); ++n) {
      int c = ids[n];

      for (int i = 0; i < nk; ++i) data2(i) = (*field_c)[i][c];

      data(0) = (*field_c)[0][c];
      for (int i = 0; i < dim; ++i) data(i + 1) = grad_c[i][c];
      for (int i = dim + 1; i < nk; ++i) data(i) = 0.0;
      dg.cell_basis(c).ChangeBasisNaturalToMy(data);

      WhetStone::DenseMatrix M;
      dg.MassMatrix(c, K, M);

      M.Multiply(data2, data3, false);
      l2norm += data2 * data3;

      data2 -= data;
      M.Multiply(data2, data3, false);
      err += data2 * data3;
    }

    double tmp = err;
    mesh->getComm()->SumAll(&tmp, &err, 1);
    tmp = l2norm;
    mesh->getComm()->SumAll(&tmp, &l2norm, 1);
    if (MyPID == 0)
      printf(
        "       sol errors: %10.6f   ||u||=%8.4f\n", std::pow(err, 0.5), std::pow(l2norm, 0.5));
    CHECK(err < 0.02);
  }
}


TEST(LIMITER_SMOOTHNESS_INDICATOR_2D)
{
  double l2norm1;
  RunTest("test/circle_quad10.exo", "orthonormalized", l2norm1);
  // RunTest("test/circle_quad20.exo");
  // RunTest("test/circle_quad40.exo");
  // RunTest("test/circle_quad80.exo");
  // RunTest("test/circle_quad160.exo");
}


/* *****************************************************************
* New limiters
***************************************************************** */
void
RunTestGaussPoints(const std::string& limiter_name)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0)
    std::cout << "\nTest: Limiters at Gauss points, method=" << limiter_name << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/circle_quad10.exo");

  // create and initialize cell-based field
  int nk(6), dim(2);
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, nk);
  auto field = Teuchos::rcp(new CompositeVector(cvs1));
  auto field_c = field->ViewComponent("cell", true);

  int order = 2;
  Teuchos::ParameterList dglist;
  dglist.set<int>("method order", order).set<std::string>("dg basis", "orthonormalized");

  WhetStone::DG_Modal dg(dglist, mesh);
  AnalyticDG08b ana(mesh, order, true);
  ana.set_shapes(true, false, false); // cone, hump, cylinder

  ana.InitialGuess(dg, *field_c, 0.0);
  field->ScatterMasterToGhosted("cell");

  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  // memory for gradient
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, dim);

  Teuchos::ParameterList plist;
  plist.set<int>("polynomial_order", 2)
    .set<bool>("limiter extension for transport", false)
    .set<std::string>("limiter", limiter_name)
    .set<std::string>("limiter stencil", "cell to all cells")
    .set<int>("limiter points", 3);

  // boundary data
  std::vector<int> bc_model(nfaces_wghost, 0);
  std::vector<double> bc_value(nfaces_wghost, 0.0);

  const auto& fmap = mesh->getMap(AmanziMesh::Entity_kind::FACE,true);
  const auto& bmap = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
  for (int bf = 0; bf < bmap.NumMyElements(); ++bf) {
    int f = fmap.LID(bmap.GID(bf));
    const auto& xf = mesh->getFaceCentroid(f);
    bc_model[f] = OPERATOR_BC_DIRICHLET;
    bc_value[f] = ana.SolutionExact(xf, 0.0);
  }

  // Apply limiter
  LimiterCellDG limiter(mesh);
  limiter.Init(plist);
  limiter.ApplyLimiterDG(field_c, dg, bc_model, bc_value);
  const auto& factor = *limiter.limiter();
  for (int c = 0; c < ncells_owned; ++c)
    for (int i = 1; i < nk; ++i) (*field_c)[i][c] *= factor[c];

  double minlim, avglim, maxlim;
  limiter.limiter()->MinValue(&minlim);
  limiter.limiter()->MeanValue(&avglim);
  limiter.limiter()->MaxValue(&maxlim);

  // calculate extrema of the limited function
  double umin(1.0e+12), umax(-1.0e+12);
  WhetStone::DenseVector data(nk);

  for (int c = 0; c < ncells_owned; ++c) {
    umin = std::min(umin, (*field_c)[0][c]);
    umax = std::max(umax, (*field_c)[0][c]);
  }

  double tmp = umin;
  mesh->getComm()->MinAll(&tmp, &umin, 1);
  tmp = umax;
  mesh->getComm()->MaxAll(&tmp, &umax, 1);
  if (MyPID == 0) {
    printf("function min/max: %10.6f %10.6f\n", umin, umax);
    printf("limiter min/avg/max: %10.6f %10.6f %10.6f\n", minlim, avglim, maxlim);
  }

  CHECK(umax > 0.84);
  CHECK(avglim > 0.97);
}


TEST(LIMITER_GAUSS_POINTS)
{
  RunTestGaussPoints("Barth-Jespersen dg");
  RunTestGaussPoints("Barth-Jespersen dg hierarchical");
}
