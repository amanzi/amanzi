/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy PK
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
#include "CompositeVector.hh"
#include "IO.hh"
#include "EnergyPressureEnthalpy_PK.hh"
#include "EnergyPressureTemperature_PK.hh"
#include "MeshFactory.hh"
#include "Operator.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "State.hh"
#include "VerboseObject.hh"
#include "WhetStoneDefs.hh"

/* ****************************************************************
* Analysis of Jacobian for the FEHM EOS.
* ************************************************************** */
template<class Energy>
void Run(const std::string& filename)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  using namespace Amanzi::Energy;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: Jacobian calculation" << std::endl;

  // read parameter list
  Teuchos::ParameterXMLFileReader xmlreader(filename);
  auto plist = Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 5, 5);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("energy");
  auto soln = Teuchos::rcp(new TreeVector());
  auto EPK = Teuchos::rcp(new Energy(pk_tree, plist, S, soln));

  EPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  EPK->Initialize();
  S->CheckAllFieldsInitialized();

  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Energy_PK", *plist));
  WriteStateStatistics(*S, *vo);

  // finite difference Jacobian
  auto u1 = soln;

  auto u0 = Teuchos::rcp(new TreeVector(*u1));
  auto f0 = Teuchos::rcp(new TreeVector(*u1));
  auto f1 = Teuchos::rcp(new TreeVector(*u1));

  int ncells, nfaces, nJ, v;
  double umax, factor, eps(1e-8), t_old(0.0), t_new(10.0);
  std::string kind;

  ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nJ = ncells + nfaces;
  WhetStone::DenseMatrix Jfd(nJ, nJ);

  u0->NormInf(&umax);
  factor = eps * umax;

  Jfd.PutScalar(0.0);
  EPK->FunctionalResidual(t_old, t_new, u0, u1, f0);

  for (int ncol = 0; ncol < nJ; ++ncol) {
    EPK->ChangedSolution();

    if (ncol < nfaces) {
      kind = "face";
      v = ncol;
    } else {
      kind = "cell";
      v = ncol - nfaces;
    }
    auto& u1_v = *u1->Data()->ViewComponent(kind);

    u1_v[0][v] += factor;
    EPK->FunctionalResidual(t_old, t_new, u0, u1, f1);
    u1_v[0][v] -= factor;

    for (int nrow = 0; nrow < nJ; ++nrow) {
      if (nrow < nfaces) {
        kind = "face";
        v = nrow;
      } else {
        kind = "cell";
        v = nrow - nfaces;
      }

      auto& f0_v = *f0->Data()->ViewComponent(kind);
      auto& f1_v = *f1->Data()->ViewComponent(kind);
      Jfd(nrow, ncol) = (f1_v[0][v] - f0_v[0][v]) / factor;
    }
  }

  // numerical Jacobian
  EPK->UpdatePreconditioner(t_old, u1, t_new - t_old);
  auto A = EPK->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW)->A();

  int num_entries;
  int* indices;
  double* values;

  WhetStone::DenseMatrix Jpk(nJ, nJ);
  Jpk.PutScalar(0.0);

  for (int row = 0; row < nJ; ++row) {
    A->ExtractMyRowView(row, num_entries, values, indices);
    for (int n = 0; n < num_entries; ++n) {
      int col = indices[n];
      Jpk(row, col) = values[n];
    }
  }

  // std::cout << Jfd << std::endl;
  // std::cout << Jpk << std::endl;
  std::cout << umax << " " << Jfd.NumRows() << std::endl;
  auto Jdiff = Jfd - Jpk;
  double jdiff = Jdiff.Norm2();
  double jfd = Jfd.Norm2();
  double jpk = Jpk.Norm2();

  std::cout << "|| Jfd - Jpk || = " << jdiff << ",  || Jfd || = " << jfd << ",  || Jpk || = " << jpk
            << std::endl;
  CHECK(jdiff / jfd < 1e-5);
}

TEST(ENERGY_JACOBIAN) {
  ::Run<Amanzi::Energy::EnergyPressureTemperature_PK>("test/energy_jacobian_pt.xml");
  ::Run<Amanzi::Energy::EnergyPressureEnthalpy_PK>("test/energy_jacobian_ph.xml");
}
