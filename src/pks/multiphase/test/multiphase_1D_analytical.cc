/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Naren Vohra (vohra@lanl.gov)
*/

/*
2 component (Hydrogen, Water), 2 phase with assumptions as in [Gharbia, Jaffre' 14]. 1D analytical solution not related to any physical scenario.
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

// Misc
#include "math.h"


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
  if (MyPID == 0) std::cout << "Test: multiphase pk, 1D analytical solution"<<std::endl;

  // read parameter list
  auto plist = Teuchos::getParametersFromXmlFile(filename);

  // create a MSTK mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh;
    
  mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 16, 1);

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
  // S->WriteDependencyGraph();
  WriteStateStatistics(*S, *vo);

  // initialize io
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  auto io = Teuchos::rcp(new OutputXDMF(iolist, mesh, true, false));

  // loop
  int iloop(0);
  double t(0.0), tend(1.0), dt(1.0e-3), dt_max(1.0e-3);
  // store Newton iterations and time step size (after successful iteration)
  std::vector<int> newton_iterations_per_step;
  std::vector<double> time_step_size;
  // error initialize
  double perr_linf_inf = 0.0, perr_linf_l1 = 0.0, perr_linf_l2 = 0.0;  

  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  double pi = M_PI;

  while (t < tend) {

    // output solution
    if (iloop % 1 == 0) {
      io->InitializeCycle(t, iloop, "");
      const auto& u0 = *S->Get<CompositeVector>("pressure_liquid").ViewComponent("cell");
      const auto& u1 = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell");
      const auto& u2 = *S->Get<CompositeVector>("molar_density_liquid").ViewComponent("cell");
      const auto& u3 = *S->Get<CompositeVector>("molar_density_gas").ViewComponent("cell");
      const auto& u4 = *S->Get<CompositeVector>("pressure_gas").ViewComponent("cell");

      io->WriteVector(*u0(0), "liquid pressure", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u1(0), "saturation_liquid", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u2(0), "liquid hydrogen", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u3(0), "gas hydrogen", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u4(0), "gas pressure", AmanziMesh::Entity_kind::CELL);
      io->FinalizeCycle();

      WriteStateStatistics(*S, *vo);
    }

    while (MPK->AdvanceStep(t, t + dt, false)) { dt /= 2.0; }

    MPK->CommitStep(t, t + dt, Tags::DEFAULT);

    // store number of Newton iterations taken (only successful iterations after possible time step reduction) 
    double iter = MPK->bdf1_dae()->number_solver_iterations();
    newton_iterations_per_step.push_back(iter);
    // store time step size
    time_step_size.push_back(dt);
    
    S->advance_cycle();
    S->advance_time(dt);
      
    // calculate error
    auto pl = *S->Get<CompositeVector>("pressure_liquid").ViewComponent("cell");
    double perr_linf_l1_tmp = 0.0, perr_linf_l2_tmp = 0.0;
    for (int c = 0; c < ncells_owned; ++c) { 
      const Amanzi::AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      double x = xc[0], y = xc[1];
      double err = std::abs( pl[0][c] + (1.0/(32*pi*pi))*std::sin(pi*x)*std::exp(-t) );
      if (perr_linf_inf < err) { perr_linf_inf = err; }
      perr_linf_l1_tmp += err * mesh->getCellVolume(c);
      perr_linf_l2_tmp += err * err * mesh->getCellVolume(c);  
    }
    perr_linf_l2_tmp = std::sqrt(perr_linf_l2_tmp);  
    if (perr_linf_l1 < perr_linf_l1_tmp) { perr_linf_l1 = perr_linf_l1_tmp; }
    if (perr_linf_l2 < perr_linf_l2_tmp) { perr_linf_l2 = perr_linf_l2_tmp; }
    // update time
    t += dt;
    dt = std::min(dt_max, dt * 2.0);
    iloop++;

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

  S->Get<CompositeVector>("ncp_fg").NormInf(&dmax);
  CHECK(dmax <= 1.0e-9);

  const auto& xg = *S->Get<CompositeVector>("molar_density_liquid").ViewComponent("cell");
  xg.MinValue(&dmin);
  CHECK(dmin >= 0.0);

  std::cout<<"pl error Linf(Linf) = "<<perr_linf_inf<<std::endl;
  std::cout<<"pl error Linf(L1) = "<<perr_linf_l1<<std::endl;
  std::cout<<"pl error Linf(L2) = "<<perr_linf_l2<<std::endl;
}


TEST(MULTIPHASE_1D_ANALYTICAL)
{
  run_test("2D", "test/multiphase_1D_analytical.xml");
}
