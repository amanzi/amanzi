/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MultiPhase PK

  Verification of Numerical Jacobian for model 1a with primary
  variables liquid pressure, mole gas fraction, and liquid saturation.
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
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IO.hh"
#include "MeshFactory.hh"
#include "Mesh.hh"
#include "State.hh"
#include "OperatorDefs.hh"

// Multiphase
#include "Multiphase_PK.hh"


/* **************************************************************** */
TEST(MULTIPHASE_MODEL_I)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Multiphase;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: jacobian for model 1a" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/multiphase_jacobian1a.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a MSTK mesh framework
  double X1(200.0);
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, X1, 20.0, 40, 4);

  // create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Multiphase_PK", *plist));

  // create a simple state populate it
  auto state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  // create a solution vector
  ParameterList pk_tree = plist->sublist("PKs").sublist("multiphase");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  auto MPK = Teuchos::rcp(new Multiphase_PK(pk_tree, plist, S, soln));

  // initialize the multiphase process kernel
  MPK->Setup();
  S->Setup();
  S->InitializeFields();

  // set up new primary variables
  int ncells_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  std::string passwd("");
  std::vector<std::string> names({ "pressure_liquid", "mole_fraction_gas", "saturation_liquid" });
  auto& pl = *S->GetW<CompositeVector>(names[0], passwd).ViewComponent("cell");
  auto& xg = *S->GetW<CompositeVector>(names[1], passwd).ViewComponent("cell");
  auto& sl = *S->GetW<CompositeVector>(names[2], passwd).ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    const auto& xc = mesh->getCellCentroid(c);
    pl[0][c] += 1e+6 * std::pow(1.0 - xc[0] / X1, 2);
    sl[0][c] -= 0.5 * std::pow(1.0 - xc[0] / X1, 3);
    xg[0][c] += 0.99 * std::pow(1.0 - xc[0] / X1, 4);
  }

  S->InitializeEvaluators();

  MPK->Initialize();
  S->CheckAllFieldsInitialized();
  WriteStateStatistics(*S, *vo);

  // verify Jacobian
  double dt(1e+7);
  S->set_time(0.0);
  S->set_initial_time(0.0);
  S->set_final_time(dt);
  S->set_last_time(dt);

  // -- create and populate a reference tree vector
  Teuchos::RCP<TreeVector> u1 = Teuchos::rcp(new TreeVector());

  for (const auto& name : names) {
    auto field = Teuchos::rcp(new TreeVector());
    u1->PushBack(field);
    field->SetData(S->GetPtrW<CompositeVector>(name, Tags::DEFAULT, passwd));
  }

  double tol(2e-5);
  auto Jfd = MPK->FiniteDifferenceJacobian(0.0, dt, u1, u1, tol);
  int nJ = Jfd.NumRows();

  // -- create Jacobian matrix
  MPK->UpdatePreconditioner(0.0, u1, dt);
  auto A = MPK->op_tree_pc()->A();

  int num_entries;
  int* indices;
  double* values;

  WhetStone::DenseMatrix Jpk(nJ, nJ);
  Jpk.PutScalar(0.0);

  for (int i = 0; i < nJ; ++i) {
    int ig = i / 3;
    int ik = i % 3;
    int nrow = ik * ncells_owned + ig;

    A->ExtractMyRowView(i, num_entries, values, indices);
    for (int n = 0; n < num_entries; ++n) {
      int jg = indices[n] / 3;
      int jk = indices[n] % 3;
      int ncol = jk * ncells_owned + jg;
      Jpk(nrow, ncol) = values[n];
    }
  }

  int ib, ie, jb, je;
  WhetStone::DenseMatrix norm_fd(3, 3), norm_pk(3, 3), norm_diff(3, 3), rel_diff(3, 3);

  auto diff = Jfd - Jpk;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ib = i * ncells_owned;
      ie = (i + 1) * ncells_owned;

      jb = j * ncells_owned;
      je = (j + 1) * ncells_owned;
      norm_fd(i, j) = Jfd.SubMatrix(ib, ie, jb, je).Norm2();
      norm_pk(i, j) = Jpk.SubMatrix(ib, ie, jb, je).Norm2();
      norm_diff(i, j) = diff.SubMatrix(ib, ie, jb, je).Norm2();

      // if (i == 0 && j == 0) { std::cout << Jpk.SubMatrix(ib, ie, jb, je) << std::endl; exit(0); }
      if (norm_pk(i, j) > 0.0) {
        rel_diff(i, j) = norm_diff(i, j) / norm_pk(i, j);
        CHECK(rel_diff(i, j) < 1e-2);
      }
    }
  }
  // std::cout << Jpk << std::endl; exit(0);

  std::cout << "Finite difference Jacobian, norms of blocks:" << std::endl;
  PrintMatrix(norm_fd, "%12.5g");

  std::cout << "Multiphase Jacobian, norms of blocks:" << std::endl;
  PrintMatrix(norm_pk, "%12.5g");

  std::cout << "Difference of Jacobians:" << std::endl;
  PrintMatrix(norm_diff, "%12.5g");

  std::cout << "Relative Difference:" << std::endl;
  PrintMatrix(rel_diff, "%12.5g");
}
