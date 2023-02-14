/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Tests for FCT.
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "LeastSquare.hh"
#include "MeshFactory.hh"

// Amanzi::Operators
#include "FCT.hh"
#include "LimiterCell.hh"
#include "OperatorDefs.hh"
#include "ReconstructionCellLinear.hh"

double
fun_field(const Amanzi::AmanziGeometry::Point& p)
{
  int d = p.dim();
  double x = p[0];
  double y = p[1];
  if (d == 2) return x * x * y + 2 * x * y * y * y;

  double z = p[2];
  return x * x * y + 2 * x * y * y * y + 3 * x * x * y * y * z * z;
}

Amanzi::AmanziGeometry::Point
fun_velocity(const Amanzi::AmanziGeometry::Point& p)
{
  int d = p.dim();
  double x = p[0];
  double y = p[1];
  if (d == 2) return Amanzi::AmanziGeometry::Point(2 * x, -2 * y); // divergence-free

  double z = p[2];
  return Amanzi::AmanziGeometry::Point(2 * x, 2 * y, -4 * z);
}


/* *****************************************************************
* Flux corrected transport in 2D
***************************************************************** */
std::pair<double, double>
RunTest(int n, int d)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = (d == 2) ?
                                    meshfactory.create(0.0, 0.0, 1.0, 1.0, n, n) :
                                    meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n, n);

  int nfaces_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int ncells_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  // initialize data
  // -- field
  auto field = Teuchos::rcp(new Epetra_MultiVector(mesh->getMap(AmanziMesh::Entity_kind::CELL,true), 1));

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    (*field)[0][c] = fun_field(xc);
  }

  // -- flux
  int dir;

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh)->SetGhosted(true)->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  auto velocity = Teuchos::rcp(new CompositeVector(cvs));
  auto& velocity_f = *velocity->ViewComponent("face");

  auto flux_lo = Teuchos::rcp(new CompositeVector(cvs));
  auto flux_ho = Teuchos::rcp(new CompositeVector(cvs));

  auto flux_exact = Teuchos::rcp(new CompositeVector(cvs));
  auto flux_numer = Teuchos::rcp(new CompositeVector(cvs));

  auto& flux_lo_f = *flux_lo->ViewComponent("face");
  auto& flux_ho_f = *flux_ho->ViewComponent("face");
  auto& flux_exact_f = *flux_exact->ViewComponent("face");

  double dt(0.03 / n);
  for (int f = 0; f < nfaces_owned; f++) {
    auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);

    if (cells.size() == 2) {
      const auto& xf = mesh->getFaceCentroid(f);
      const auto& normal = mesh->getFaceNormal(f, cells[0], &dir);
      const auto& xc = (dir == 1) ? mesh->getCellCentroid(cells[0]) : mesh->getCellCentroid(cells[1]);

      // low-order and high-order fluxes are from 1st to 2nd cell
      double tmp = (fun_velocity(xf) * normal) * dt;
      velocity_f[0][f] = tmp * dir;

      flux_lo_f[0][f] = tmp * fun_field(xc); // default upwind
      flux_ho_f[0][f] = tmp * fun_field(xf);

      flux_exact_f[0][f] = flux_ho_f[0][f];
    }
  }

  // -- boundary conditions
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  auto& bc_model = bc->bc_model();
  auto& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(1.0 - xf[0]) < 1e-6 || fabs(xf[1]) < 1e-6 ||
        fabs(1.0 - xf[1]) < 1e-6 || fabs(xf[d]) < 1e-6 || fabs(1.0 - xf[d]) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = fun_field(xf);
    }
  }

  // compute reconstruction
  Teuchos::ParameterList plist;
  plist.set<int>("polynomial_order", 1)
    .set<bool>("limiter extension for transport", false)
    .set<std::string>("limiter", "Barth-Jespersen")
    .set<std::string>("limiter stencil", "cell to all cells");

  auto lifting = Teuchos::rcp(new ReconstructionCellLinear(mesh));
  lifting->Init(plist);
  lifting->Compute(field);

  // create limiter (Do we need it?)
  auto limiter = Teuchos::rcp(new LimiterCell(mesh));
  limiter->Init(plist, velocity->ViewComponent("face"));

  // FCT
  FCT fct(mesh, mesh, limiter, field);
  fct.Compute(*flux_lo, *flux_ho, *bc, *flux_numer);

  // Error calculations
  double err[1];
  std::pair<double, double> out;

  flux_lo->Update(1.0, *flux_exact, -1.0);
  flux_lo->Norm2(&err[0]);
  out.first = err[0] * n;

  flux_numer->Update(1.0, *flux_exact, -1.0);
  flux_numer->Norm2(&err[0]);
  out.second = err[0] * n;

  if (comm->MyPID() == 0) {
    std::cout << "error lower order: " << out.first << std::endl;
    std::cout << "error high order: " << out.second << std::endl;
    std::cout << "mean limiter: " << fct.get_alpha_mean() << std::endl;
  }

  return out;
}


TEST(FCT_2D)
{
  auto a1 = RunTest(10, 2);
  auto a2 = RunTest(20, 2);
  auto a3 = RunTest(40, 2);

  std::vector<double> h({ 1.0 / 10, 1.0 / 20, 1.0 / 40 });
  std::vector<double> err1({ a1.first, a2.first, a3.first });
  double rate = Amanzi::Utils::bestLSfit(h, err1);
  CHECK(rate < 1.0);
  std::cout << "\nError convergence rate: " << rate << std::endl;

  double err2 = std::max({ a1.second, a2.second, a3.second });
  CHECK(err2 < 2.0e-4);
  std::cout << "Deviation from high-order flux: " << err2 << std::endl << std::endl;
}


TEST(FCT_3D)
{
  auto a1 = RunTest(8, 3);
  auto a2 = RunTest(16, 3);
  auto a3 = RunTest(32, 3);

  std::vector<double> h({ 1.0 / 8, 1.0 / 16, 1.0 / 32 });
  std::vector<double> err1({ a1.first, a2.first, a3.first });
  double rate = Amanzi::Utils::bestLSfit(h, err1);
  CHECK(rate > 1.0);
  std::cout << "\nError convergence rate: " << rate << std::endl;

  double err2 = std::max({ a1.second, a2.second, a3.second });
  CHECK(err2 < 2.0e-4);
  std::cout << "Deviation from high-order flux: " << err2 << std::endl;
}
