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
#include "Analytic01b.hh"

/* *****************************************************************
* Modify meshes
* **************************************************************** */
double random(double range) {
  return 2 * range * (double(rand()) / RAND_MAX - 0.5);
}


void
ModifyMesh(Amanzi::AmanziMesh::Mesh& mesh,
           const Amanzi::AmanziGeometry::Point& scale,
           const Amanzi::AmanziGeometry::Point& shift)
{
  int d = shift.dim();
  int nnodes_owned = mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE,
                                         Amanzi::AmanziMesh::Parallel_kind::OWNED);

  for (int n = 0; n < nnodes_owned; ++n) {
    auto xp = mesh.getNodeCoordinate(n);
    if (std::fabs(xp[0]) > 1e-5 && std::fabs(1 - xp[0]) > 1e-5 &&
        std::fabs(xp[1]) > 1e-5 && std::fabs(1 - xp[1]) > 1e-5 &&
        std::fabs(xp[d - 1]) > 1e-5 && std::fabs(1 - xp[d - 1]) > 1e-5) {
      xp[0] += random(0.15 / std::pow((double)nnodes_owned, 1.0 / d));
      xp[1] += random(0.15 / std::pow((double)nnodes_owned, 1.0 / d));
      if (d == 3) xp[2] += random(0.15 / std::pow((double)nnodes_owned, 1.0 / d));
    }
    for (int i = 0; i < d; ++i) {
      xp[i] = scale[i] * xp[i] + shift[i];
    }
    mesh.setNodeCoordinate(n, xp);
  }
  mesh.recacheGeometry();
}


void
ModifyMesh1(Amanzi::AmanziMesh::Mesh& mesh)
{
  int nnodes_owned = mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE,
                                         Amanzi::AmanziMesh::Parallel_kind::OWNED);

  for (int n = 0; n < nnodes_owned; ++n) {
    auto xp = mesh.getNodeCoordinate(n);
    double s = 0.06 * std::sin(M_PI * xp[0]) * std::sin(M_PI * xp[1]);
    xp[0] += s;
    xp[1] += s;
    mesh.setNodeCoordinate(n, xp);
  }
  mesh.recacheGeometry();
}


void
ModifyMesh2(Amanzi::AmanziMesh::Mesh& mesh)
{
  double a(0.5);
  int nnodes_owned = mesh.getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE,
                                         Amanzi::AmanziMesh::Parallel_kind::OWNED);

  for (int n = 0; n < nnodes_owned; ++n) {
    auto xp = mesh.getNodeCoordinate(n);
    auto yp = xp;
    xp[0] = yp[0] * (1.0 + a * yp[0]) / (1.0 + a);
    xp[1] = yp[1] * (1.0 + a * yp[1]) / (1.0 + a);
    mesh.setNodeCoordinate(n, xp);
  }
  mesh.recacheGeometry();
}


/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* **************************************************************** */
template<class Analytic>
void
TestDiffusionMultiMesh(int d, double tol,
                       int level = 1,
                       const std::string& filename = "",
                       const std::string& filename3 = "")
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: " << d << "D multi mesh problem: three meshes, file=" << filename << "\n";

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_multi_mesh_cuts.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  ParameterList region_list = plist.sublist("regions" + to_string(d));
  auto gm = Teuchos::rcp(new GeometricModel(d, region_list, *comm));

  // create meshes
  int n = level;
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh1, mesh2, mesh3;
  if (d == 2) {
    if (filename == "") {
      mesh1 = meshfactory.create(0.0, 0.0, 1.0, 1.0, 4 * n, 5 * n);
      mesh2 = meshfactory.create(0.0, 0.0, 1.0, 1.0, 5 * n, 4 * n);
      mesh3 = meshfactory.create(0.0, 0.0, 1.0, 1.0, 6 * n, 8 * n);
    } else if (filename != "") {
      mesh1 = meshfactory.create(filename);
      mesh2 = meshfactory.create(filename);
      mesh3 = meshfactory.create(filename3);
    }
    // ModifyMesh2(*mesh2);
    ModifyMesh(*mesh1, AmanziGeometry::Point(0.5, 0.5), AmanziGeometry::Point(0.0, 0.5));
    ModifyMesh(*mesh2, AmanziGeometry::Point(0.5, 0.5), AmanziGeometry::Point(0.0, 0.0));
    ModifyMesh(*mesh3, AmanziGeometry::Point(0.5, 1.0), AmanziGeometry::Point(0.5, 0.0));
  } else {
    mesh1 = meshfactory.create(0.0, 0.5, 0.0, 0.5, 1.0, 1.0, 5 * n, 5 * n, 5 * n);
    mesh2 = meshfactory.create(0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 7 * n, 7 * n, 5 * n);
    mesh3 = meshfactory.create(0.5, 0.0, 0.0, 1.0, 1.0, 1.0, 9 * n, 8 * n, 7 * n);
  }

  // extract mesh sets and continue with curvilinear modifications
  if (d == 2) {
    mesh1->getSetEntities("cut1top", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    mesh3->getSetEntities("cut1top", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    mesh2->getSetEntities("cut1bottom", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    mesh3->getSetEntities("cut1bottom", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    mesh1->getSetEntities("cut2", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    mesh2->getSetEntities("cut2", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

    ModifyMesh1(*mesh1);
    ModifyMesh1(*mesh2);
    ModifyMesh1(*mesh3);
  }

  int ncells1_owned = mesh1->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces1_wghost = mesh1->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  int ncells2_owned = mesh2->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces2_wghost = mesh2->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  int ncells3_owned = mesh3->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces3_wghost = mesh3->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  // populate diffusion coefficient
  WhetStone::Tensor Knull;
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator");

  Analytic ana1(mesh1, 0), ana2(mesh2, 0), ana3(mesh3, 0);

  auto Kc1 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells1_owned; ++c) {
    const Point& xc = mesh1->getCellCentroid(c);
    Kc1->push_back(ana1.TensorDiffusivity(xc, 0.0));
  }

  auto Kc2 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells2_owned; ++c) {
    const Point& xc = mesh2->getCellCentroid(c);
    Kc2->push_back(ana2.TensorDiffusivity(xc, 0.0));
  }

  auto Kc3 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells3_owned; ++c) {
    const Point& xc = mesh3->getCellCentroid(c);
    Kc3->push_back(ana3.TensorDiffusivity(xc, 0.0));
  }

  // populate boundary data
  auto bc1 =
    Teuchos::rcp(new BCs(mesh1, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc1_model = bc1->bc_model();
  std::vector<double>& bc1_value = bc1->bc_value();

  for (int f = 0; f < nfaces1_wghost; ++f) {
    const Point& xf = mesh1->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        (d == 3 && (fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6))) {
      bc1_model[f] = OPERATOR_BC_DIRICHLET;
      bc1_value[f] = ana1.pressure_exact(xf, 0.0);
    }
  }

  auto bc2 =
    Teuchos::rcp(new BCs(mesh2, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc2_model = bc2->bc_model();
  std::vector<double>& bc2_value = bc2->bc_value();

  for (int f = 0; f < nfaces2_wghost; ++f) {
    const Point& xf = mesh2->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[1]) < 1e-6 ||
        (d == 3 && (fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6))) {
      bc2_model[f] = OPERATOR_BC_DIRICHLET;
      bc2_value[f] = ana2.pressure_exact(xf, 0.0);
    }
  }

  auto bc3 =
    Teuchos::rcp(new BCs(mesh3, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc3_model = bc3->bc_model();
  std::vector<double>& bc3_value = bc3->bc_value();

  for (int f = 0; f < nfaces3_wghost; ++f) {
    const Point& xf = mesh3->getFaceCentroid(f);
    if (fabs(xf[1]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        (d == 3 && (fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6))) {
      bc3_model[f] = OPERATOR_BC_DIRICHLET;
      bc3_value[f] = ana3.pressure_exact(xf, 0.0);
    }
  }

  // create state and PDE
  auto S = Teuchos::rcp(new State());
  S->RegisterMesh("mesh1", mesh1);
  S->RegisterMesh("mesh2", mesh2);
  S->RegisterMesh("mesh3", mesh3);

  auto pde = Teuchos::rcp(new PDE_DiffusionMultiMesh(op_list));
  pde->SetVariableTensorCoefficient("mesh1", Kc1);
  pde->SetVariableTensorCoefficient("mesh2", Kc2);
  pde->SetVariableTensorCoefficient("mesh3", Kc3);

  pde->Init(S);

  // create a source terms
  auto op = pde->get_matrix();
  TreeVector rhs(op->DomainMap());
  rhs.SubVector(0)->SetData(op->get_operator_block(0, 0)->rhs());
  rhs.SubVector(1)->SetData(op->get_operator_block(1, 1)->rhs());
  rhs.SubVector(2)->SetData(op->get_operator_block(2, 2)->rhs());

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

  auto& rhs3 = *op->get_operator_block(2, 2)->rhs()->ViewComponent("cell");
  for (int c = 0; c < ncells3_owned; c++) {
    const Point& xc = mesh3->getCellCentroid(c);
    rhs3[0][c] = ana3.source_exact(xc, 0.0) * mesh3->getCellVolume(c);
  }

  // create diffusion problem
  pde->SetBCs("mesh1", bc1, bc1);
  pde->SetBCs("mesh2", bc2, bc2);
  pde->SetBCs("mesh3", bc3, bc3);

  pde->UpdateMatrices(Teuchos::null, Teuchos::null);
  pde->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");

  auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(
    "Hypre AMG", plist.sublist("preconditioners"), "PCG", plist.sublist("solvers"));
  op->set_inverse_parameters(inv_list);
  op->InitializeInverse();
  op->ComputeInverse();
  // std::cout << *op->A() << std::endl; exit(0);

  auto sol = Teuchos::rcp(new TreeVector(op->DomainMap()));
  op->ApplyInverse(rhs, *sol);

  auto flux = Teuchos::rcp(new TreeVector(op->DomainMap()));
  pde->UpdateFlux(sol.ptr(), flux.ptr());

  // error calculation
  Epetra_MultiVector& p1 = *sol->SubVector(0)->Data()->ViewComponent("cell");
  Epetra_MultiVector& p2 = *sol->SubVector(1)->Data()->ViewComponent("cell");
  Epetra_MultiVector& p3 = *sol->SubVector(2)->Data()->ViewComponent("cell");

  double pnorm1, l2_err1, inf_err1, pnorm2, l2_err2, inf_err2, pnorm3, l2_err3, inf_err3;
  ana1.ComputeCellError(p1, 0.0, pnorm1, l2_err1, inf_err1);
  ana2.ComputeCellError(p2, 0.0, pnorm2, l2_err2, inf_err2);
  ana3.ComputeCellError(p3, 0.0, pnorm3, l2_err3, inf_err3);
  CHECK(l2_err1 < tol);
  CHECK(l2_err2 < tol);
  CHECK(l2_err3 < tol);

  double p1min, p1max, p2min, p2max, p3min, p3max;
  p1.MinValue(&p1min);
  p2.MinValue(&p2min);
  p3.MinValue(&p3min);

  p1.MaxValue(&p1max);
  p2.MaxValue(&p2max);
  p3.MaxValue(&p3max);

  if (MyPID == 0) {
    double err = std::sqrt((l2_err1 * l2_err1 + l2_err2 * l2_err2 + l2_err3 * l2_err3) / 
                           (pnorm1 * pnorm1 + pnorm2 * pnorm2 + pnorm3 * pnorm3));
    printf("Total: L2(p) =%10.7f  Inf =%9.6f\n", err, std::max({ inf_err1, inf_err2, inf_err3 }));

    printf("L2(p) =%9.6f  Inf =%9.6f  min/max = %9.6f %9.6f  #cells=%d\n",
           l2_err1 / pnorm1,
           inf_err1,
           p1min,
           p1max,
           ncells1_owned);
    printf("L2(p) =%9.6f  Inf =%9.6f  min/max = %9.6f %9.6f  #cells=%d\n",
           l2_err2 / pnorm2,
           inf_err2,
           p2min,
           p2max,
           ncells2_owned);
    printf("L2(p) =%9.6f  Inf =%9.6f  min/max = %9.6f %9.6f  #cells=%d\n",
           l2_err3 / pnorm3,
           inf_err3,
           p3min,
           p3max,
           ncells3_owned);
  }

  auto& flux1 = *flux->SubVector(0)->Data()->ViewComponent("face", true);
  auto& flux2 = *flux->SubVector(1)->Data()->ViewComponent("face", true);
  auto& flux3 = *flux->SubVector(2)->Data()->ViewComponent("face", true);
  
  double unorm1, unorm2, unorm3;
  ana1.ComputeFaceError(flux1, 0.0, unorm1, l2_err1, inf_err1);
  ana2.ComputeFaceError(flux2, 0.0, unorm2, l2_err2, inf_err2);
  ana3.ComputeFaceError(flux3, 0.0, unorm3, l2_err3, inf_err3);

  if (MyPID == 0) {
    double err = std::sqrt((l2_err1 * l2_err1 + l2_err2 * l2_err2 + l2_err3 * l2_err3) / 
                           (unorm1 * unorm1 + unorm2 * unorm2 + unorm3 * unorm3));
    printf("Total: L2(u) =%10.7f  Inf =%9.6f\n", err, std::max({ inf_err1, inf_err2, inf_err3 }));

    printf("L2(u) =%10.7f  Inf =%9.6f  |u|=%9.6f\n", l2_err1 / unorm1, inf_err1, unorm1);
    printf("L2(u) =%10.7f  Inf =%9.6f  |u|=%9.6f\n", l2_err2 / unorm2, inf_err2, unorm2);
    printf("L2(u) =%10.7f  Inf =%9.6f  |u|=%9.6f\n", l2_err3 / unorm3, inf_err3, unorm3);
  }

  // -- lemma 4.3
  double l2_err(0.0);
  auto interface_weights = pde->get_interface_weights();
  std::vector<Teuchos::RCP<AmanziMesh::Mesh>> meshes({ mesh1, mesh2, mesh3 });

  for (int k = 0; k < 3; ++k) {
    int i1, i2, dir, l1, l2;
    double sum(0.0);
    std::string cutname;
    if (k == 0) {
      i1 = 0;
      i2 = 1;
      cutname = "cut2";
      l1 = 4;
      l2 = 5;
    } else if (k == 1) {
      i1 = 0;
      i2 = 2;
      cutname = "cut1top";
      l1 = 0;
      l2 = 1;
    } else if (k == 2) {
      i1 = 1;
      i2 = 2;
      cutname = "cut1bottom";
      l1 = 2;
      l2 = 3;
    }

    auto& p1_f = *sol->SubVector(i1)->Data()->ViewComponent("face");
    auto& p2_f = *sol->SubVector(i2)->Data()->ViewComponent("face");
    auto& flux1 = *flux->SubVector(i1)->Data()->ViewComponent("face", true);
    auto& flux2 = *flux->SubVector(i2)->Data()->ViewComponent("face", true);

    const auto& block1 = meshes[i1]->getSetEntities(cutname, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    const auto& block2 = meshes[i2]->getSetEntities(cutname, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

    double tot_flux1(0.0), tot_flux2(0.0);
    for (int f : block1) {
      int c = getFaceOnBoundaryInternalCell(*meshes[i1], f);
      meshes[i1]->getFaceNormal(f, c, &dir);
      tot_flux1 += flux1[0][f] * dir;
    }
    ana1.GlobalOp("sum", &tot_flux1, 1);

    for (int f : block2) {
      int c = getFaceOnBoundaryInternalCell(*meshes[i2], f);
      meshes[i2]->getFaceNormal(f, c, &dir);
      tot_flux2 += flux2[0][f] * dir;
    }
    ana1.GlobalOp("sum", &tot_flux2, 1);
    printf("Total interface flux: %9.6f, err =%9.6f\n", tot_flux1, tot_flux1 + tot_flux2);
    CHECK(std::fabs(tot_flux1 + tot_flux2) < 1e-12);

    for (int f : block1) {
      int c = getFaceOnBoundaryInternalCell(*meshes[i1], f);
      meshes[i1]->getFaceNormal(f, c, &dir);
   
      double pavg(p1_f[0][f]);
      for (auto data : interface_weights[l1][f]) {
        int f2 = data.first;
        pavg += data.second * p2_f[0][f2];
      }
      sum += pavg * flux1[0][f] * dir;
    }

    for (int f : block2) {
      int c = getFaceOnBoundaryInternalCell(*meshes[i2], f);
      meshes[i2]->getFaceNormal(f, c, &dir);

      double pavg(p2_f[0][f]);
      for (auto data : interface_weights[l2][f]) {
        int f1 = data.first;
        pavg += data.second * p1_f[0][f1];
      }
      sum += pavg * flux2[0][f] * dir;
    }
    printf("Lemma 4.3, %12.9f >= 0.0\n", sum / 2);
    CHECK(sum >= -1e-12);

    // interface pressure error
    for (int f : block1) {
      const auto& xf = meshes[i1]->getFaceCentroid(f);
      double area = meshes[i1]->getFaceArea(f);
      double tmp = ana1.pressure_exact(xf, 0.0);
      l2_err += std::pow(tmp - p1_f[0][f], 2) * area;
    }
    for (int f : block2) {
      const auto& xf = meshes[i2]->getFaceCentroid(f);
      double area = meshes[i2]->getFaceArea(f);
      double tmp = ana2.pressure_exact(xf, 0.0);
      l2_err += std::pow(tmp - p2_f[0][f], 2) * area;
    }
  }
  printf("Interface absolute error: L2(p)=%12.9f\n", std::sqrt(l2_err));

  // i/o
  Teuchos::ParameterList iolist;
  iolist.set<std::string>("file name base", "plot1");
  OutputXDMF io1(iolist, mesh1, true, false);

  iolist.set<std::string>("file name base", "plot2");
  OutputXDMF io2(iolist, mesh2, true, false);

  iolist.set<std::string>("file name base", "plot3");
  OutputXDMF io3(iolist, mesh3, true, false);

  int loop(0);
  double t(0.0);
  io1.InitializeCycle(t, loop, "");
  io1.WriteVector(*p1(0), "solution1", AmanziMesh::Entity_kind::CELL);
  io1.FinalizeCycle();

  io2.InitializeCycle(t, loop, "");
  io2.WriteVector(*p2(0), "solution2", AmanziMesh::Entity_kind::CELL);
  io2.FinalizeCycle();

  io3.InitializeCycle(t, loop, "");
  io3.WriteVector(*p3(0), "solution3", AmanziMesh::Entity_kind::CELL);
  io3.FinalizeCycle();
}


TEST(OPERATOR_DIFFUSION_TWO_MESH_PROBLEM)
{
  TestDiffusionMultiMesh<Analytic01>(2, 5e-2, 1); // lemma 4.3 gives zero
  TestDiffusionMultiMesh<Analytic01>(2, 5e-2, 1, "test/poly8.exo", "test/median7x8.exo");
  TestDiffusionMultiMesh<Analytic01b>(3, 1e-1, 1);
}
