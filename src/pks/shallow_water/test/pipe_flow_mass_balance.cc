/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
 Pipe flow PK. Test to check mass conservation.

 */

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IO.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"

// Amanzi::ShallowWater
#include "ShallowWater_PK.hh"
#include "ShallowWater_Helper.hh"

// Amanzi::PipeFlow
#include "PipeFlow_PK.hh"

#include "EvaluatorIndependent.hh"
#include "Evaluator.hh"

using namespace Amanzi;

TEST(PIPE_FLOW_1D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::ShallowWater;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 1D pipe flow mass conservation" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/pipe_flow_mass_balance.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  // create a mesh
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 10.0, 1.0, 25, 1);

  // create a state
  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterMesh("pipe", mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

  Teuchos::ParameterList pf_list = plist->sublist("PKs").sublist("pipe flow");

  // create a pipe flow PK
  PipeFlow_PK PFPK(pf_list, plist, S, soln);
  PFPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);

  PFPK.Initialize();

  S->CheckAllFieldsInitialized();

  S->Get<CompositeVector>("pipe-bathymetry").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("pipe-velocity").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("pipe-total_depth").ScatterMasterToGhosted("cell");

  // create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("PipeFlow", pf_list));
  WriteStateStatistics(*S, *vo);

  // advance in time
  double t_old(0.0), t_new(0.0), dt;

  // initialize io
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  auto io = Teuchos::rcp(new OutputXDMF(iolist, mesh, true, false));

  std::string passwd("state");

  int iter = 0;
  double Tend = 1.0;

  // for mass calculation
  double total_mass_initial_tmp = 0.0, total_mass_final_tmp = 0.0, total_mass_initial,
         total_mass_final;

  int ncells_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                          Amanzi::AmanziMesh::Parallel_kind::OWNED);

  auto& wac_init = *S->GetW<CompositeVector>("pipe-wetted_area", Tags::DEFAULT, "")
                      .ViewComponent("cell"); // wetted area
  auto& htc_init = *S->GetW<CompositeVector>("pipe-total_depth", Tags::DEFAULT, "")
                      .ViewComponent("cell"); // total depth for junction
  auto& Bc = *S->GetW<CompositeVector>("pipe-bathymetry", Tags::DEFAULT, "")
                .ViewComponent("cell"); // bathymetry for junction
  const auto& dirc = *S->Get<CompositeVector>("pipe-direction")
                        .ViewComponent("cell"); // direction to identify junction

  for (int c = 0; c < ncells_owned; ++c) {
    if (std::abs(dirc[0][c] - 1.0) < 1.e-12) {
      total_mass_initial_tmp +=
        wac_init[0][c] * mesh->getCellVolume(c); // length of cell = volume of cell
    } else if (std::abs(dirc[0][c] - 0.0) < 1.e-12) {
      total_mass_initial_tmp += (htc_init[0][c] - Bc[0][c]) * mesh->getCellVolume(c);
    }
  }

  mesh->getComm()->SumAll(&total_mass_initial_tmp, &total_mass_initial, 1);

  while (t_new < Tend) {
    // cycle 1, time t
    dt = PFPK.get_dt();

    t_new = t_old + dt;

    PFPK.AdvanceStep(t_old, t_new);
    PFPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;

    // output data
    if (iter % 200 == 0) {
      io->InitializeCycle(t_new, iter, "");

      const auto& u0 = *S->Get<CompositeVector>("pipe-total_depth").ViewComponent("cell");
      const auto& u1 = *S->Get<CompositeVector>("pipe-velocity").ViewComponent("cell");
      const auto& u2 = *S->Get<CompositeVector>("pipe-bathymetry").ViewComponent("cell");
      const auto& u3 = *S->Get<CompositeVector>("pipe-discharge").ViewComponent("cell");
      const auto& u4 = *S->Get<CompositeVector>("pipe-wetted_area").ViewComponent("cell");
      const auto& u5 = *S->Get<CompositeVector>("pipe-direction").ViewComponent("cell");

      io->WriteVector(*u0(0), "pipe-total_depth", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u1(0), "pipe-velocity x", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u1(1), "pipe-velocity y", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u2(0), "pipe-bathymetry", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u3(0), "pipe-discharge x", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u3(1), "pipe-discharge y", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u4(0), "pipe-wetted_area", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u5(0), "pipe-direction", AmanziMesh::Entity_kind::CELL);
      io->FinalizeCycle();
    }
  }

  auto& wac = *S->GetW<CompositeVector>("pipe-wetted_area", Tags::DEFAULT, "")
                 .ViewComponent("cell"); // wetted area
  auto& htc = *S->GetW<CompositeVector>("pipe-total_depth", Tags::DEFAULT, "")
                 .ViewComponent("cell"); // total depth for junction

  for (int c = 0; c < ncells_owned; ++c) {
    if (std::abs(dirc[0][c] - 1.0) < 1.e-12) {
      total_mass_final_tmp += wac[0][c] * mesh->getCellVolume(c); // length of cell = volume of cell
    } else if (std::abs(dirc[0][c] - 0.0) < 1.e-12) {
      total_mass_final_tmp += (htc[0][c] - Bc[0][c]) * mesh->getCellVolume(c);
    }
  }

  mesh->getComm()->SumAll(&total_mass_final_tmp, &total_mass_final, 1);

  std::cout.precision(14);

  if (MyPID == 0) {
    std::cout << "Calculated total initial mass = " << total_mass_initial << std::endl;
    std::cout << "Calculated total final mass = " << total_mass_final << std::endl;
    std::cout << "Difference = " << total_mass_final - total_mass_initial << std::endl;
    std::cout << "Calculated mass through BCs = " << 0.0 << std::endl;

    CHECK(std::abs((total_mass_final - total_mass_initial) - 0.0) < 1.e-12);
  }

  WriteStateStatistics(*S, *vo);
}
