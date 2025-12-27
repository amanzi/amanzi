/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "InverseFactory.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"
#include "State.hh"
#include "TreeVector.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "OperatorDefs.hh"
#include "PDE_DiffusionMultiMesh.hh"

#include "Analytic00.hh"
#include "Analytic01.hh"
#include "Analytic02.hh"
#include "Analytic08.hh"

/* *****************************************************************
* Modify meshes
* **************************************************************** */
double random(double range) {
  return 2 * range * (double(rand()) / RAND_MAX - 0.5);
}

void
ModifyMesh1(Amanzi::AmanziMesh::Mesh& mesh,
           const Amanzi::AmanziGeometry::Point& scale,
           const Amanzi::AmanziGeometry::Point& shift)
{
  int d = shift.dim();
  int nnodes_owned = mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE,
                                         Amanzi::AmanziMesh::Parallel_kind::OWNED);

  for (int n = 0; n < nnodes_owned; ++n) {
    auto xp = mesh.getNodeCoordinate(n);
    for (int i = 0; i < d; ++i) {
      xp[i] = scale[i] * xp[i] + shift[i];
    }
    mesh.setNodeCoordinate(n, xp);
  }
  mesh.recacheGeometry();
}

void
ModifyMesh2(Amanzi::AmanziMesh::Mesh& mesh)
{
  int nnodes_owned = mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE,
                                         Amanzi::AmanziMesh::Parallel_kind::OWNED);

  for (int n = 0; n < nnodes_owned; ++n) {
    auto xp = mesh.getNodeCoordinate(n);
    // xp[1] = xp[1] * (2.0 - xp[1]);
    // xp[1] = xp[1] + 0.06 * std::sin(4.0 * M_PI * xp[1]);
    if (std::fabs(xp[1]) > 1e-5 && std::fabs(1 - xp[1]) > 1e-5 &&  std::fabs(0.5 - xp[0]) > 1e-5)
      xp[1] += random(std::sqrt(0.03 / nnodes_owned));
    mesh.setNodeCoordinate(n, xp);
  }
  mesh.recacheGeometry();
}


void ModifyMesh3(Amanzi::AmanziMesh::Mesh& mesh, int nx, int id)
{
  int nnodes_owned = mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE,
                                         Amanzi::AmanziMesh::Parallel_kind::OWNED);

  double amp(0.04);

  for (int n = 0; n < nnodes_owned; ++n) {
    int i = n % (nx + 1);
    auto xp = mesh.getNodeCoordinate(n);
    double s = 0.1 * std::sin(2 * M_PI * xp[1]);

    if (id == 1) {
      double xmax = 0.5 + amp * std::sin(4 * M_PI * xp[1]);
      xp[0] = i * xmax / nx;
    } else if (id == 2) {
      double xmax = 0.5 + amp * std::sin(4 * M_PI * xp[1]);
      xp[0] = xmax + i * (1.0 - xmax) / nx;
    }
    mesh.setNodeCoordinate(n, xp);
  }
  mesh.recacheGeometry();
}


/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* **************************************************************** */
template<class Analytic>
void
TestDiffusionMultiMesh(int d, 
                       double tol, int icase, int level = 1, 
                       const std::string& filename1 = "",
                       const std::string& filename2 = "")
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: " << d << "D multi mesh problem\n";

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_multi_mesh.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  int n = plist.sublist("PK operator").sublist("diffusion operator").get<int>("number of particles");
  n *= level;
  plist.sublist("PK operator").sublist("diffusion operator").set<int>("number of particles", n);

  ParameterList region_list = plist.sublist("regions" + std::to_string(d));
  auto gm = Teuchos::rcp(new GeometricModel(d, region_list, *comm));

  // create meshese
  n = level;
  int n1x(4 * n), n1y(8 * n), n1z(6 * n);
  int n2x(5 * n), n2y(10 * n), n2z(6 * n);

  auto mlist = Teuchos::rcp(new Teuchos::ParameterList(plist.sublist("mesh")));
  MeshFactory meshfactory(comm, gm, mlist);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh1, mesh2;
  if (filename1 == "") {
    if (d == 2) {
      mesh1 = meshfactory.create(0.0, 0.0, 0.5, 1.0, n1x, n1y);
      mesh2 = meshfactory.create(0.5, 0.0, 1.0, 1.0, n2x, n2y);
    } else {
      mesh1 = meshfactory.create(0.0, 0.0, 0.0, 0.5, 1.0, 1.0, n1x, n1y, n1z);
      mesh2 = meshfactory.create(0.5, 0.0, 0.0, 1.0, 1.0, 1.0, n2x, n2y, n2z);
    }
  } else {
    mesh1 = meshfactory.create(filename1);
    mesh2 = meshfactory.create(filename2);
  }

  int ncells1_owned = mesh1->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces1_wghost = mesh1->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  int ncells2_owned = mesh2->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces2_wghost = mesh2->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  // modify meshes 
  cEntity_ID_View block1, block2;
  if (filename1 != "") {
    ModifyMesh1(*mesh1, AmanziGeometry::Point(0.5, 1.0), AmanziGeometry::Point(0.0, 0.0));
    ModifyMesh1(*mesh2, AmanziGeometry::Point(0.5, 1.0), AmanziGeometry::Point(0.5, 0.0));

    block1 = mesh1->getSetEntities("cut1", AmanziMesh::Entity_kind::FACE, Parallel_kind::OWNED);
    block2 = mesh2->getSetEntities("cut1", AmanziMesh::Entity_kind::FACE, Parallel_kind::OWNED);
  } else if (icase == 1) {
    ModifyMesh2(*mesh2);

    block1 = mesh1->getSetEntities("cut1", AmanziMesh::Entity_kind::FACE, Parallel_kind::OWNED);
    block2 = mesh2->getSetEntities("cut1", AmanziMesh::Entity_kind::FACE, Parallel_kind::OWNED);
  } else if (icase == 2) {
    block1 = mesh1->getSetEntities("cut1", AmanziMesh::Entity_kind::FACE, Parallel_kind::OWNED);
    block2 = mesh2->getSetEntities("cut1", AmanziMesh::Entity_kind::FACE, Parallel_kind::OWNED);

    ModifyMesh3(*mesh1, n1x, 1);
    ModifyMesh3(*mesh2, n2x, 2);
  }

  // populate diffusion coefficient
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator");

  Analytic ana1(mesh1, 3), ana2(mesh2, 3);

  auto Kc1 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells1_owned; ++c) {
    const Point& xc = mesh1->getCellCentroid(c);
    auto Kc = ana1.TensorDiffusivity(xc, 0.0);
    Kc1->push_back(Kc);
  }

  auto Kc2 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells2_owned; ++c) {
    const Point& xc = mesh2->getCellCentroid(c);
    auto Kc = ana2.TensorDiffusivity(xc, 0.0);
    Kc2->push_back(Kc);
  }

  // populate boundary data
  auto bc1 = Teuchos::rcp(new BCs(mesh1, Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc1_model = bc1->bc_model();
  std::vector<double>& bc1_value = bc1->bc_value();

  for (int f = 0; f < nfaces1_wghost; ++f) {
    const Point& xf = mesh1->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || 
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[d - 1]) < 1e-6 || fabs(xf[d - 1] - 1.0) < 1e-6) {
      bc1_model[f] = OPERATOR_BC_DIRICHLET;
      bc1_value[f] = ana1.pressure_exact(xf, 0.0);
    }
  }

  auto bc2 = Teuchos::rcp(new BCs(mesh2, Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc2_model = bc2->bc_model();
  std::vector<double>& bc2_value = bc2->bc_value();

  for (int f = 0; f < nfaces2_wghost; ++f) {
    const Point& xf = mesh2->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || 
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[d - 1]) < 1e-6 || fabs(xf[d - 1] - 1.0) < 1e-6) {
      bc2_model[f] = OPERATOR_BC_DIRICHLET;
      bc2_value[f] = ana2.pressure_exact(xf, 0.0);
    }
  }

  // create state and PDE
  auto S = Teuchos::rcp(new State());
  S->RegisterMesh("mesh1", mesh1);
  S->RegisterMesh("mesh2", mesh2);

  auto pde = Teuchos::rcp(new PDE_DiffusionMultiMesh(op_list));
  pde->SetVariableTensorCoefficient("mesh1", Kc1);
  pde->SetVariableTensorCoefficient("mesh2", Kc2);

  pde->Init(S);

  // create a source terms
  auto op = pde->get_matrix();
  TreeVector rhs(op->DomainMap());
  rhs.SubVector(0)->SetData(op->get_operator_block(0, 0)->rhs());
  rhs.SubVector(1)->SetData(op->get_operator_block(1, 1)->rhs());

  auto& rhs1 = *op->get_operator_block(0, 0)->rhs()->ViewComponent("cell");
  for (int c = 0; c < ncells1_owned; c++) {
    const Point& xc = mesh1->getCellCentroid(c);
    rhs1[0][c] = ana1.source_exact(xc, 0.0) * mesh1->getCellVolume(c);
  }

  auto& rhs2 = *op->get_operator_block(1, 1)->rhs()->ViewComponent("cell");
  for (int c = 0; c < ncells2_owned; c++) {
    const Point& xc = mesh2->getCellCentroid(c);
    rhs2[0][c] = ana2.source_exact(xc, 0.0) * mesh2->getCellVolume(c);
  }

  // create diffusion problem
  pde->SetBCs("mesh1", bc1, bc1);
  pde->SetBCs("mesh2", bc2, bc2);

  pde->UpdateMatrices(Teuchos::null, Teuchos::null);
  pde->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");

  auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(
    "Hypre AMG", plist.sublist("preconditioners"), "PCG", plist.sublist("solvers"));
  op->set_inverse_parameters(inv_list);
  op->InitializeInverse();
  op->ComputeInverse();
  // std::cout << setw(16) << setprecision(12) << *op->A() << std::endl; exit(0);

  auto sol = Teuchos::rcp(new TreeVector(op->DomainMap()));
  op->ApplyInverse(rhs, *sol);

  auto flux = Teuchos::rcp(new TreeVector(op->DomainMap()));
  pde->UpdateFlux(sol.ptr(), flux.ptr());

  // error calculation
  auto& p1 = *sol->SubVector(0)->Data()->ViewComponent("cell");
  auto& p2 = *sol->SubVector(1)->Data()->ViewComponent("cell");

  double pnorm1, l2_err1, inf_err1, pnorm2, l2_err2, inf_err2, err;
  ana1.ComputeCellError(p1, 0.0, pnorm1, l2_err1, inf_err1);
  ana2.ComputeCellError(p2, 0.0, pnorm2, l2_err2, inf_err2);
  CHECK(l2_err1 < tol);

  double p1min, p1max, p2min, p2max;
  p1.MinValue(&p1min);
  p2.MinValue(&p2min);

  p1.MaxValue(&p1max);
  p2.MaxValue(&p2max);

  if (MyPID == 0) {
    err = std::sqrt((l2_err1 * l2_err1 + l2_err2 * l2_err2) / (pnorm1 * pnorm1 + pnorm2 * pnorm2));
    printf("Total: L2(p) =%10.7f  Inf =%10.7f\n", err, std::max(inf_err1, inf_err2));

    l2_err1 /= pnorm1;
    l2_err2 /= pnorm2;
    printf("L2(p) =%10.7f  Inf =%9.6f  min/max = %9.6f %9.6f  #cells=%d\n",
           l2_err1,
           inf_err1,
           p1min,
           p1max,
           ncells1_owned);
    printf("L2(p) =%10.7f  Inf =%9.6f  min/max = %9.6f %9.6f  #cells=%d\n",
           l2_err2,
           inf_err2,
           p2min,
           p2max,
           ncells2_owned);
  }

  auto& flux1 = *flux->SubVector(0)->Data()->ViewComponent("face", true);
  auto& flux2 = *flux->SubVector(1)->Data()->ViewComponent("face", true);
  
  double unorm1, unorm2;
  ana1.ComputeFaceError(flux1, 0.0, unorm1, l2_err1, inf_err1);
  ana2.ComputeFaceError(flux2, 0.0, unorm2, l2_err2, inf_err2);

  if (MyPID == 0) {
    err = std::sqrt((l2_err1 * l2_err1 + l2_err2 * l2_err2) / (unorm1 * unorm1 + unorm2 * unorm2));
    printf("Total: L2(u) =%10.7f  Inf =%10.7f\n", err, std::max(inf_err1, inf_err2));

    l2_err1 /= unorm1;
    l2_err2 /= unorm2;
    printf("L2(u) =%10.7f  Inf =%9.6f  |u|=%9.6f\n", l2_err1, inf_err1, unorm1);
    printf("L2(u) =%10.7f  Inf =%9.6f  |u|=%9.6f\n", l2_err2, inf_err2, unorm2);
  }

  int dir;
  double tot_flux1(0.0), tot_flux2(0.0);
  for (int f : block1) {
    int c = getFaceOnBoundaryInternalCell(*mesh1, f);
    mesh1->getFaceNormal(f, c, &dir);
    tot_flux1 += flux1[0][f] * dir;
  }

  for (int f : block2) {
    int c = getFaceOnBoundaryInternalCell(*mesh2, f);
    mesh2->getFaceNormal(f, c, &dir);
    tot_flux2 += flux2[0][f] * dir;
  }
  printf("Total interface flux: %12.9f, err =%12.9f\n", tot_flux1, tot_flux1 + tot_flux2);
  CHECK(std::fabs(tot_flux1 + tot_flux2) < 1e-11);

  // -- lemma 3.3
  double sum1(0.0), sum2(0.0);
  auto interface_weights = pde->get_interface_weights();

  auto& p1_f = *sol->SubVector(0)->Data()->ViewComponent("face");
  auto& p2_f = *sol->SubVector(1)->Data()->ViewComponent("face");

  for (int f : block1) {
    int c = getFaceOnBoundaryInternalCell(*mesh1, f);
    mesh1->getFaceNormal(f, c, &dir);
    double area = (d == 3) ? mesh1->getFaceArea(f) : 1.0;
   
    double pavg(p1_f[0][f]), jump(p1_f[0][f]);
    for (auto data : interface_weights[0][f]) {
      int f2 = data.first;
      pavg += data.second * p2_f[0][f2];
      jump -= data.second * p2_f[0][f2];
    }
    sum1 += pavg * flux1[0][f] * dir;
    sum2 += jump * jump * area;
  }

  for (int f : block2) {
    int c = getFaceOnBoundaryInternalCell(*mesh2, f);
    mesh2->getFaceNormal(f, c, &dir);
    double area = (d == 3) ? mesh2->getFaceArea(f) : 1.0;

    double pavg(p2_f[0][f]), jump(p2_f[0][f]);
    for (auto data : interface_weights[1][f]) {
      int f1 = data.first;
      pavg += data.second * p1_f[0][f1];
      jump -= data.second * p1_f[0][f1];
    }
    sum1 += pavg * flux2[0][f] * dir;
    sum2 += jump * jump * area;
  }
  printf("Lemma 3.3, %12.9f >= 0.0, jump = %12.9f \n", sum1 / 2, sum2);
  CHECK(sum1 >= -1e-12);
  CHECK(std::fabs(sum1 / 2 - sum2) < 1e-12);

  // convergence theorem
  sum1 = 0.0;
  sum2 = 0.0;
  for (int f : block1) {
    int c = getFaceOnBoundaryInternalCell(*mesh1, f);
    mesh1->getFaceNormal(f, c, &dir);
   
    const auto& xf1 = mesh1->getFaceCentroid(f);
    double pavg(ana1.pressure_exact(xf1, 0.0));
    double jump(pavg);
    for (auto data : interface_weights[0][f]) {
      int f2 = data.first;
      const auto& xf2 = mesh2->getFaceCentroid(f2);
      pavg += data.second * ana2.pressure_exact(xf2, 0.0);
      jump -= data.second * ana2.pressure_exact(xf2, 0.0);
    }
    sum1 += pavg * flux1[0][f] * dir;
    sum2 += std::fabs(jump);
  }

  for (int f : block2) {
    int c = getFaceOnBoundaryInternalCell(*mesh2, f);
    mesh2->getFaceNormal(f, c, &dir);

    const auto& xf2 = mesh2->getFaceCentroid(f);
    double pavg(ana2.pressure_exact(xf2, 0.0));
    double jump(pavg);
    for (auto data : interface_weights[1][f]) {
      int f1 = data.first;
      const auto& xf1 = mesh1->getFaceCentroid(f1);
      pavg += data.second * ana1.pressure_exact(xf1, 0.0);
      jump -= data.second * ana1.pressure_exact(xf1, 0.0);
    }
    sum1 += pavg * flux2[0][f] * dir;
    sum2 += std::fabs(jump);
  }
  printf("Theorem, integral = %12.9f \n", sum1  / sum2);

  // interface pressure error
  double l2_err(0.0), pnorm(0.0);
  for (int f : block1) {
    const auto& xf = mesh1->getFaceCentroid(f);
    double area = mesh1->getFaceArea(f);

    double tmp = ana1.pressure_exact(xf, 0.0);
    l2_err += std::pow(tmp - p1_f[0][f], 2) * area;
    pnorm += tmp * tmp * area;
  }
  for (int f : block2) {
    const auto& xf = mesh2->getFaceCentroid(f);
    double area = mesh2->getFaceArea(f);

    double tmp = ana2.pressure_exact(xf, 0.0);
    l2_err += std::pow(tmp - p2_f[0][f], 2) * area;
    pnorm += tmp * tmp * area;
  }
  printf("Interface absolute error: L2(p)=%12.9f\n", std::sqrt(l2_err));

  // i/o
  Teuchos::ParameterList iolist;
  iolist.set<std::string>("file name base", "plot1");
  OutputXDMF io1(iolist, mesh1, true, false);

  iolist.set<std::string>("file name base", "plot2");
  OutputXDMF io2(iolist, mesh2, true, false);

  int loop(0);
  double t(0.0);
  io1.InitializeCycle(t, loop, "");
  io1.WriteVector(*p1(0), "solution1", AmanziMesh::Entity_kind::CELL);
  io1.FinalizeCycle();

  io2.InitializeCycle(t, loop, "");
  io2.WriteVector(*p2(0), "solution2", AmanziMesh::Entity_kind::CELL);
  io2.FinalizeCycle();
}


TEST(OPERATOR_DIFFUSION_TWO_MESH_PROBLEM)
{
  TestDiffusionMultiMesh<Analytic01>(2, 2e-1, 1, 1, "test/median8_filtered.exo", "test/poly8.exo");
  TestDiffusionMultiMesh<Analytic08>(2, 7e-3, 2);
  TestDiffusionMultiMesh<Analytic00>(2, 1e-2, 1);
  TestDiffusionMultiMesh<Analytic00>(3, 1e-2, 2);
}
