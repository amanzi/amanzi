/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Naren Vohra (vohra@lanl.gov)
*/

/*
2 component Hydrogen (H) and water (W) with assumptions/model as in [Gharbia, Jaffre' 14]. Primary variables are pressure liquid, saturation liquid, and mole fraction of hydrogen in liquid phase.
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
void
run_test(const std::string& filename)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Multiphase;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: multiphase pk, model Jaffre with pl, sl, xg as primary variables"<<std::endl;

  // read parameter list
  auto plist = Teuchos::getParametersFromXmlFile(filename);

  // create a MSTK mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

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
  
  /*
  // work-around
  Key key("mass_density_gas");
  S->Require<CompositeVector, CompositeVectorSpace>(key, Tags::DEFAULT, key)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  */

  MPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // initialize the multiphase process kernel
  MPK->Initialize();
  S->CheckAllFieldsInitialized();
  // S->WriteDependencyGraph();
  WriteStateStatistics(*S, *vo);

  // initialize io
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  auto io = Teuchos::rcp(new OutputXDMF(iolist, mesh, true, false));

  // loop
  int iloop(0);
  double t(0.0), tend(3.14e+12), dt(1.57e+7), dt_max(1.57e+12); // Tend = 1000,000 years, dt = 5000 years
  // store Newton iterations and time step size (after successful iteration)
  std::vector<int> newton_iterations_per_step;
  std::vector<double> time_step_size;
  
  while (t < tend && iloop < 100000) {
    while (MPK->AdvanceStep(t, t + dt, false)) { dt /= 2.0; }

    MPK->CommitStep(t, t + dt, Tags::DEFAULT);

    // store number of Newton iterations taken (only successful iterations after possible time step reduction) 
    double iter = MPK->bdf1_dae()->number_solver_iterations();
    newton_iterations_per_step.push_back(iter);
    // store time step size
    time_step_size.push_back(dt);
    
    S->advance_cycle();
    S->advance_time(dt);
    t += dt;
    dt = std::min(dt_max, dt * 1.6);
    iloop++;

    // output solution
    if (iloop % 2 == 0) {
      io->InitializeCycle(t, iloop, "");
      const auto& u0 = *S->Get<CompositeVector>("pressure_liquid").ViewComponent("cell");
      const auto& u1 = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell");
      const auto& u2 = *S->Get<CompositeVector>("mole_fraction_gas").ViewComponent("cell");
      const auto& u3 = *S->Get<CompositeVector>("pressure_gas").ViewComponent("cell");

      io->WriteVector(*u0(0), "liquid pressure", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u1(0), "saturation_liquid", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u2(0), "mole fraction gas", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u3(0), "gas pressure", AmanziMesh::Entity_kind::CELL);
      io->FinalizeCycle();

      WriteStateStatistics(*S, *vo);
    }
  }

  WriteStateStatistics(*S, *vo);

  // write iteration output to text file
  std::ofstream outFile("iterations_per_time_step.txt");
  for (int i = 0; i < newton_iterations_per_step.size(); ++i) { 
    outFile << newton_iterations_per_step[i] << "," << time_step_size[i] << std::endl;
  }
  outFile.close(); 

  // verification
  double dmin, dmax;
  const auto& sl = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell");
  sl.MinValue(&dmin);
  sl.MaxValue(&dmax);
  CHECK(dmin >= 0.0 && dmax <= 1.0 && dmin <= 1.0);

  const auto& xg = *S->Get<CompositeVector>("mole_fraction_gas").ViewComponent("cell");
  xg.MinValue(&dmin);
  xg.MaxValue(&dmax);
  CHECK(dmin >= 0.0 && dmax >= 0.0 && dmax <= 1.0 && dmin <= 1.0);
}


TEST(MULTIPHASE_JAFFRE_2P2C)
{
  run_test("test/multiphase_jaffre_xg.xml");
}
