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
    for (int i = 0; i < d; ++i) {
      xp[i] = scale[i] * xp[i] + shift[i];
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
TestDiffusionMultiMesh(int d, double tol, const std::string& filename = "")
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D multi mesh problem: three meshes\n";

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_multi_mesh_cuts.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  ParameterList region_list = plist.sublist("regions" + to_string(d));
  auto gm = Teuchos::rcp(new GeometricModel(d, region_list, *comm));

  // create meshese
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh1, mesh2, mesh3;
  if (d == 2 && filename == "") {
    mesh1 = meshfactory.create(0.0, 0.5, 0.5, 1.0, 5, 5);
    mesh2 = meshfactory.create(0.0, 0.0, 0.5, 0.5, 5, 5);
    mesh3 = meshfactory.create(0.5, 0.0, 1.0, 1.0, 10, 10);
  } else if (d == 2 && filename != "") {
    mesh1 = meshfactory.create(filename);
    mesh2 = meshfactory.create(filename);
    mesh3 = meshfactory.create(filename);
    ModifyMesh(*mesh1, AmanziGeometry::Point(0.5, 0.5), AmanziGeometry::Point(0.0, 0.5));
    ModifyMesh(*mesh2, AmanziGeometry::Point(0.5, 0.5), AmanziGeometry::Point(0.0, 0.0));
    ModifyMesh(*mesh3, AmanziGeometry::Point(0.5, 1.0), AmanziGeometry::Point(0.5, 0.0));
  } else {
    mesh1 = meshfactory.create(0.0, 0.5, 0.0, 0.5, 1.0, 1.0, 5, 5, 5);
    mesh2 = meshfactory.create(0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 7, 7, 5);
    mesh3 = meshfactory.create(0.5, 0.0, 0.0, 1.0, 1.0, 1.0, 9, 8, 7);
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

  Analytic ana1(mesh1), ana2(mesh2), ana3(mesh3);

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
  TreeVector rhs(op->DomainMap()), sol(op->DomainMap());
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

  op->ApplyInverse(rhs, sol);

  // error calculation
  Epetra_MultiVector& p1 = *sol.SubVector(0)->Data()->ViewComponent("cell");
  Epetra_MultiVector& p2 = *sol.SubVector(1)->Data()->ViewComponent("cell");
  Epetra_MultiVector& p3 = *sol.SubVector(2)->Data()->ViewComponent("cell");

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
    l2_err1 /= pnorm1;
    l2_err2 /= pnorm2;
    l2_err3 /= pnorm3;
    printf("L2 =%9.6f  Inf =%9.6f  min/max = %9.6f %9.6f  #cells=%d\n",
           l2_err1,
           inf_err1,
           p1min,
           p1max,
           ncells1_owned);
    printf("L2 =%9.6f  Inf =%9.6f  min/max = %9.6f %9.6f  #cells=%d\n",
           l2_err2,
           inf_err2,
           p2min,
           p2max,
           ncells2_owned);
    printf("L2 =%9.6f  Inf =%9.6f  min/max = %9.6f %9.6f  #cells=%d\n",
           l2_err3,
           inf_err3,
           p3min,
           p3max,
           ncells3_owned);
  }

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
  TestDiffusionMultiMesh<Analytic01>(2, 5e-2, "test/median7x8.exo");
  TestDiffusionMultiMesh<Analytic01>(2, 5e-2);
  TestDiffusionMultiMesh<Analytic01b>(3, 1e-1);
}
