/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy

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
#include "eos_reg.hh"
#include "evaluators_reg.hh"
#include "FlowEnergyPH_PK.hh"
#include "IO.hh"
#include "MeshFactory.hh"
#include "models_flow_reg.hh"
#include "models_energy_reg.hh"
#include "PK_Factory.hh"
#include "pks_energy_reg.hh"
#include "pks_flow_reg.hh"
#include "pks_mpc_reg.hh"
#include "State.hh"
#include "VerboseObject.hh"
#include "WhetStoneDefs.hh"
#include "FlowEnergyPH_PK.hh"


/* ****************************************************************
* Runs cooling problem
* ************************************************************** */
TEST(MPC_TWO_PHASE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Evaluators;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: flow-energy cooling problem" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/mpc_two_phasePH_heating.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 2.0, 1.0, 2, 1);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  PKFactory factory;
  auto pk_tree = plist->sublist("cycle driver").sublist("time periods").sublist("TP 0").sublist("PK tree");
  auto soln = Teuchos::rcp(new TreeVector());
  auto mpc = factory.CreatePK("transient:mpc1", pk_tree, plist, S, soln);

  mpc->parseParameterList();
  mpc->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  mpc->Initialize();
  S->CheckAllFieldsInitialized();

  WriteStateStatistics(*S);

  // variable timestepping
  const auto& p = *S->Get<CompositeVector>("pressure", Tags::DEFAULT).ViewComponent("cell");
  const auto& T = *S->Get<CompositeVector>("temperature", Tags::DEFAULT).ViewComponent("cell");
  const auto& rho = *S->Get<CompositeVector>("mass_density_liquid", Tags::DEFAULT).ViewComponent("cell");
  const auto& mu = *S->Get<CompositeVector>("viscosity_liquid", Tags::DEFAULT).ViewComponent("cell");
  const auto& state = *S->Get<CompositeVector>("thermodynamic_state", Tags::DEFAULT).ViewComponent("cell");

  CHECK(state[(int)TSPH_t::RGN][0] == 1);
  CHECK(state[(int)TSPH_t::RGN][1] == 1);

  int itrs(0);
  double t(0.0), dt(5.0), t1(50.0e+3);
  while (t < t1) {
    bool fail = mpc->AdvanceStep(t, t + dt, false);
    if (fail) {
      dt /= 2.0;
    } else {
      mpc->CommitStep(t, t + dt, Tags::DEFAULT);
      dt = std::min(50.0, dt * 1.01); 
      // std::cout << "PT: " << t + dt << " " << p[0][0] << " " << state[(int)TSPH_t::X][0] << " " << state[(int)TSPH_t::RGN][0] 
      //                               << " " << p[0][1] << " " << state[(int)TSPH_t::X][1] << " " << state[(int)TSPH_t::RGN][1] << std::endl;

      if (itrs % 20 == 0) WriteStateStatistics(*S);
    }

    t += dt;
    itrs++;
  }

  CHECK(state[(int)TSPH_t::RGN][0] == 2);
  CHECK(state[(int)TSPH_t::RGN][1] == 2);
}
