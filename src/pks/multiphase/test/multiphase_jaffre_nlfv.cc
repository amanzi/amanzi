/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Naren Vohra (vohra@lanl.gov)
*/

/*
  MoMas benchmark example: 2 component Hydrogen (H) and water (W) to show gas
  phase appearance/disappearance with assumptions/model as in [Gharbia, Jaffre' 14].
  Primary variables are pressure liquid, saturation liquid, and molar density
  of hydrogen in liquid phase.
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

const std::vector<std::vector<double>> fv_data = {
  { 1.5700e+11, 1.00007e+06, 1.00007e+06 }, { 8.3210e+11, 1.05948e+06, 1.18695e+06 },
  { 1.5072e+12, 1.10941e+06, 1.30413e+06 }, { 2.2922e+12, 1.13316e+06, 1.36986e+06 },
  { 3.0772e+12, 1.14144e+06, 1.40667e+06 }, { 3.8622e+12, 1.13799e+06, 1.42517e+06 },
  { 4.6472e+12, 1.12557e+06, 1.43071e+06 }, { 5.4322e+12, 1.10765e+06, 1.42786e+06 },
  { 6.2172e+12, 1.08837e+06, 1.42119e+06 }, { 7.0022e+12, 1.07047e+06, 1.41372e+06 },
  { 7.7872e+12, 1.05537e+06, 1.40707e+06 }, { 8.5722e+12, 1.04329e+06, 1.40176e+06 },
  { 9.3572e+12, 1.03382e+06, 1.39764e+06 }, { 1.0142e+13, 1.02656e+06, 1.39459e+06 },
  { 1.0927e+13, 1.02092e+06, 1.39226e+06 }, { 1.1712e+13, 1.01651e+06, 1.39046e+06 },
  { 1.2497e+13, 1.01306e+06, 1.38906e+06 }, { 1.3282e+13, 1.01035e+06, 1.38797e+06 },
  { 1.4067e+13, 1.00822e+06, 1.38712e+06 }, { 1.4852e+13, 1.00654e+06, 1.38645e+06 },
  { 1.5637e+13, 1.00522e+06, 1.38593e+06 }, { 1.6422e+13, 8.16296e+05, 1.16264e+06 },
  { 1.7207e+13, 7.69158e+05, 1.08735e+06 }, { 1.7992e+13, 7.72753e+05, 1.05068e+06 },
  { 1.8777e+13, 8.03882e+05, 1.02903e+06 }, { 1.9562e+13, 8.48750e+05, 1.01289e+06 },
  { 2.0347e+13, 9.00026e+05, 9.96174e+05 }, { 2.1132e+13, 9.87602e+05, 9.87602e+05 },
  { 2.1917e+13, 1.00011e+06, 1.00011e+06 }, { 2.2702e+13, 1.00010e+06, 1.00010e+06 },
  { 2.3487e+13, 1.00008e+06, 1.00008e+06 }, { 2.4272e+13, 1.00007e+06, 1.00007e+06 },
  { 2.5057e+13, 1.00006e+06, 1.00006e+06 }, { 2.5842e+13, 1.00005e+06, 1.00005e+06 },
  { 2.6627e+13, 1.00005e+06, 1.00005e+06 }, { 2.7412e+13, 1.00004e+06, 1.00004e+06 },
  { 2.8197e+13, 1.00004e+06, 1.00004e+06 }, { 2.8982e+13, 1.00003e+06, 1.00003e+06 },
  { 2.9767e+13, 1.00003e+06, 1.00003e+06 }, { 3.0552e+13, 1.00002e+06, 1.00002e+06 },
  { 3.1337e+13, 1.00002e+06, 1.00002e+06 },
};

/* **************************************************************** */
void
run_test(const std::string& domain, const std::string& filename)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Multiphase;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: multiphase pk, model Jaffre, domain=" << domain << std::endl;

  int d = (domain == "2D") ? 2 : 3;

  // read parameter list
  auto plist = Teuchos::getParametersFromXmlFile(filename);

  // create a MSTK mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(d, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh;
  mesh = meshfactory.create(0.0, 0.0, 200.0, 20.0, 100, 1);

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

  // Tend = 1000,000 years, dt = 5000 years
  int iloop(0);
  double t(0.0), tend(3.14e+13), dt(1.57e+11), dt_max(1.57e+11);

  // store Newton iterations and timestep size (after successful iteration)
  std::vector<int> newton_iterations_per_step;
  std::vector<double> time_step_size;

  std::vector<std::vector<double>> nlfv_data;

  while (t < tend && iloop < 1000) {
    while (MPK->AdvanceStep(t, t + dt, false)) { dt /= 2.0; }

    MPK->CommitStep(t, t + dt, Tags::DEFAULT);

    // store number of Newton iterations taken (only successful iterations after possible timestep reduction)
    double iter = MPK->bdf1_dae()->number_solver_iterations();
    newton_iterations_per_step.push_back(iter);
    // store timestep size
    time_step_size.push_back(dt);

    S->advance_cycle();
    S->advance_time(dt);
    t += dt;
    dt = std::min(dt_max, dt * 1.6);
    iloop++;

    // output solution
    const auto& u0 = *S->Get<CompositeVector>("pressure_liquid").ViewComponent("cell");
    const auto& u4 = *S->Get<CompositeVector>("pressure_gas").ViewComponent("cell");

    std::vector<double> data({ t, u0[0][0], u4[0][0] });
    nlfv_data.push_back(data);
  }
  WriteStateStatistics(*S, *vo);

  // compute error between two curves
  double w, pl0, pg0, pl1, pg1, pref(1e+6), err_l(0.0), err_g(0.0);
  for (int i = 0; i < nlfv_data.size(); ++i) {
    t = nlfv_data[i][0];
    pl0 = nlfv_data[i][1];
    pg0 = nlfv_data[i][2];

    for (int n = 0; n < fv_data.size() - 1; ++n) {
      if (fv_data[n][0] <= t && t <= fv_data[n + 1][0]) {
        w = (t - fv_data[n][0]) / (fv_data[n + 1][0] - fv_data[n][0]);
        pl1 = fv_data[n + 1][1] * w + fv_data[n][1] * (1 - w);
        pg1 = fv_data[n + 1][2] * w + fv_data[n][2] * (1 - w);
        break;
      }
    }
    err_l += std::fabs(pl1 - pl0) / pref / nlfv_data.size();
    err_g += std::fabs(pg1 - pg0) / pref / nlfv_data.size();
  }

  CHECK(err_l <= 2e-3);
  CHECK(err_g <= 3e-3);
}


TEST(MULTIPHASE_JAFFRE_NLFV)
{
  run_test("2D", "test/multiphase_jaffre_nlfv.xml");
}
