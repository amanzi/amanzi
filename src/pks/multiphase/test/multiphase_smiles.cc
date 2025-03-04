/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MultiPhase

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
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshExtractedManifold.hh"
#include "State.hh"
#include "OperatorDefs.hh"
#include "OutputXDMF.hh"

// Multiphase
#include "Multiphase_PK.hh"


/* **************************************************************** */
TEST(MULTIPHASE_SMILES)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Multiphase;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: multiphase Smiles diffusion test" << std::endl;

  int d = 2;

  // read parameter list
  auto plist = Teuchos::getParametersFromXmlFile("test/multiphase_smiles.xml");

  // create a MSTK mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(d, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 10.0, 1.0, 51, 1);

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

  // work-around
  Key key("mass_density_gas");
  S->Require<CompositeVector, CompositeVectorSpace>(key, Tags::DEFAULT, key)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  MPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // initialize the multiphase process kernel
  MPK->Initialize();
  S->CheckAllFieldsInitialized();
  WriteStateStatistics(*S, *vo);

  // initialize io
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  auto io = Teuchos::rcp(new OutputXDMF(iolist, mesh, true, false));

  const auto& tccl = *S->Get<CompositeVector>("total_component_concentration_liquid", Tags::DEFAULT)
                        .ViewComponent("cell");
  const auto& tccg = *S->Get<CompositeVector>("total_component_concentration_gas", Tags::DEFAULT)
                        .ViewComponent("cell");

  // loop
  int iloop(0);
  double t(0.0), tend(3.1558e+8), dt(3.1558e+6), dt_max(2 * 3.1558e+6);
  while (t < tend) {
    while (MPK->AdvanceStep(t, t + dt, false)) { dt /= 10; }

    MPK->CommitStep(t, t + dt, Tags::DEFAULT);
    S->advance_cycle();

    S->advance_time(dt);
    t += dt;
    dt = std::min(dt_max, dt * 1.2);
    dt = std::min(tend - t, dt);

    // output solution
    iloop++;
    if (iloop % 5 == 0) WriteStateStatistics(*S, *vo);

    if (iloop % 5000 == 0) {
      io->InitializeCycle(t, iloop, "");
      const auto& u2 = *S->Get<CompositeVector>("molar_density_liquid").ViewComponent("cell");
      const auto& u3 = *S->Get<CompositeVector>("molar_density_gas").ViewComponent("cell");

      io->WriteVector(*u2(0), "liquid tritium", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u3(0), "gas tritium", AmanziMesh::Entity_kind::CELL);
      io->FinalizeCycle();
    }
    printf("AAA: %g %16.10g %16.10g ... %16.10g   %16.10g\n",
           t,
           tccl[0][0],
           tccl[0][1],
           tccl[0][25],
           tccg[0][25]);
  }

  WriteStateStatistics(*S, *vo);

  double sl(0.01);
  double tccl_eff = tccl[0][25] + tccg[0][25] * (1 - sl) / sl;
  printf("Normalized effective liquid concentration: %g \n", tccl_eff * 1e+5);
  CHECK_CLOSE(0.14, tccl_eff / 1e-5, 0.002);
}
