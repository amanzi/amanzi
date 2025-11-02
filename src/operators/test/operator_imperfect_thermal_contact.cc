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
void
TestImperfectContact(int d,
                     const std::string& filename1 = "",
                     const std::string& filename2 = "",
                     int ibc = 0)
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
  std::string xmlFileName = "test/operator_imperfect_thermal_contact.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  ParameterList region_list = plist.sublist("regions" + std::to_string(d));
  auto gm = Teuchos::rcp(new GeometricModel(d, region_list, *comm));

  // create meshes
  int n1x(10), n1y(4), n1z(6);
  int n2x(10), n2y(5), n2z(6);
  auto mlist = Teuchos::rcp(new Teuchos::ParameterList(plist.sublist("mesh")));
  MeshFactory meshfactory(comm, gm, mlist);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh1, mesh2;
  if (filename1 != "") {
    mesh1 = meshfactory.create(filename1);
    mesh2 = meshfactory.create(filename2);

    ModifyMesh(*mesh1, AmanziGeometry::Point(0.5, 1.0), AmanziGeometry::Point(0.0, 0.0));
    ModifyMesh(*mesh2, AmanziGeometry::Point(0.5, 1.0), AmanziGeometry::Point(0.5, 0.0));
  } else if (d == 2) {
    mesh1 = meshfactory.create(0.0, 0.0, 0.5, 1.0, n1x, n1y);
    mesh2 = meshfactory.create(0.5, 0.0, 1.0, 1.0, n2x, n2y);
  } else {
    mesh1 = meshfactory.create(0.0, 0.0, 0.0, 0.5, 1.0, 1.0, n1x, n1y, n1z);
    mesh2 = meshfactory.create(0.5, 0.0, 0.0, 1.0, 1.0, 1.0, n2x, n2y, n2z);
  }

  int ncells1_owned = mesh1->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces1_wghost = mesh1->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  int ncells2_owned = mesh2->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int nfaces2_wghost = mesh2->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);

  // populate diffusion coefficient
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator");

  double k1(1.0), k2(1.0), L1(0.5), L2(0.5), R;
  R = 1.0 / op_list.sublist("interfaces").sublist("interface 1").sublist("contact conductance").sublist("function-constant").get<double>("value");
  double q = 1.0 / (R + L1 / k1 + L2 / k2);
  // double jump = R * q;

  WhetStone::Tensor K(d, 1);
  auto Kc1 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells1_owned; ++c) {
    const Point& xc = mesh1->getCellCentroid(c);
    K(0, 0) = k1;
    Kc1->push_back(K);
  }

  auto Kc2 = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells2_owned; ++c) {
    const Point& xc = mesh2->getCellCentroid(c);
    K(0, 0) = k2;
    Kc2->push_back(K);
  }

  // populate boundary data
  auto bc1 = Teuchos::rcp(new BCs(mesh1, Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc1_model = bc1->bc_model();
  std::vector<double>& bc1_value = bc1->bc_value();

  for (int f = 0; f < nfaces1_wghost; ++f) {
    const Point& xf = mesh1->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6) {
      bc1_model[f] = OPERATOR_BC_DIRICHLET;
      bc1_value[f] = xf[0];
    }
    else if (ibc == 1 && (fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6)) {
      bc1_model[f] = OPERATOR_BC_DIRICHLET;
      bc1_value[f] = (q / k1) * xf[0];
    }
  }

  auto bc2 = Teuchos::rcp(new BCs(mesh2, Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc2_model = bc2->bc_model();
  std::vector<double>& bc2_value = bc2->bc_value();

  for (int f = 0; f < nfaces2_wghost; ++f) {
    const Point& xf = mesh2->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6) {
      bc2_model[f] = OPERATOR_BC_DIRICHLET;
      bc2_value[f] = xf[0];
    } 
    else if (ibc == 1 && (fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6)) {
      bc2_model[f] = OPERATOR_BC_DIRICHLET;
      bc2_value[f] = 1.0 + (q / k2) * (xf[0] - L1 - L2);
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

  // error analysis
  // -- norms
  auto& p1 = *sol->SubVector(0)->Data()->ViewComponent("cell");
  auto& p2 = *sol->SubVector(1)->Data()->ViewComponent("cell");

  double pnorm1(0.0), l2_err1(0.0), inf_err1(0.0), pexact, err;
  double pnorm2(0.0), l2_err2(0.0), inf_err2(0.0);

  for (int c = 0; c < ncells1_owned; ++c) {
    const auto& xc = mesh1->getCellCentroid(c);
    pexact = (q / k1) * xc[0];

    err = pexact - p1[0][c];
    inf_err1 = std::max(inf_err1, std::fabs(err));
    l2_err1 += err * err * mesh1->getCellVolume(c);
    pnorm1 += pexact * pexact * mesh1->getCellVolume(c);
  }

  for (int c = 0; c < ncells2_owned; ++c) {
    const auto& xc = mesh2->getCellCentroid(c);
    pexact = xc[0];
    pexact = 1.0 + (q / k2) * (xc[0] - L1 - L2);

    err = pexact - p2[0][c];
    inf_err2 = std::max(inf_err2, std::fabs(err));
    l2_err2 += err * err * mesh2->getCellVolume(c);
    pnorm2 += pexact * pexact * mesh2->getCellVolume(c);
  }
  CHECK(l2_err1 < 1e-12 && l2_err2 < 1e-12);

  double p1min, p1max, p2min, p2max;
  p1.MinValue(&p1min);
  p2.MinValue(&p2min);

  p1.MaxValue(&p1max);
  p2.MaxValue(&p2max);

  if (MyPID == 0) {
    err = std::sqrt((l2_err1 * l2_err1 + l2_err2 * l2_err2) / (pnorm1 * pnorm1 + pnorm2 * pnorm2));
    printf("Total: L2(p) =%10.7f  Inf =%9.6f\n", err, std::max(inf_err1, inf_err2));

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

  // -- lemma 4.3
  auto block1 = mesh1->getSetEntities("cut1", AmanziMesh::Entity_kind::FACE, Parallel_kind::OWNED);
  auto block2 = mesh2->getSetEntities("cut1", AmanziMesh::Entity_kind::FACE, Parallel_kind::OWNED);

  auto& flux1 = *flux->SubVector(0)->Data()->ViewComponent("face", true);
  auto& flux2 = *flux->SubVector(1)->Data()->ViewComponent("face", true);

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
  TestImperfectContact(2);
  TestImperfectContact(2, "test/median7x8.exo", "test/poly8.exo", 1);
}
