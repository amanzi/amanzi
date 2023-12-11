/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Naren Vohra (vohra@lanl.gov)
*/

/*
2 component (Hydrogen, Water), 2 phase with assumptions as in [Gharbia, Jaffre' 14]. 1D analytical solution not related to any physical scenario.
Manufactured sol:
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
#include "LeastSquare.hh"

// Multiphase
#include "Multiphase_PK.hh"

// Misc
#include "math.h"


/* **************************************************************** */
void
exact_field(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                      Epetra_MultiVector& pl_ex,
                      Epetra_MultiVector& sl_ex,
                      double t)
{
double pi = M_PI;

int ncells_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_kind::OWNED);

for (int c = 0; c < ncells_owned; c++) {
    
    const Amanzi::AmanziGeometry::Point& xc = mesh->getCellCentroid(c);

    double x = xc[0], y = xc[1];

    //pl_ex[0][c] = -(1/(32.0*pi*pi)) * std::sin(pi*x) * std::exp(-1000.0*t);
    //sl_ex[0][c] = (1/8.0) * std::sin(pi*x) * std::exp(-1000.0*t);

    pl_ex[0][c] = x;
    sl_ex[0][c] = 1.0;
  }

}

std::tuple<double, double> 
run_test(int M, double dt0,  const std::string& filename)
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
    
  mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, M, 1);
  //mesh = meshfactory.create("test/struct.exo");

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
  double t(0.0), tend(1.0), dt(dt0), dt_max(dt0);

  // store Newton iterations and time step size (after successful iteration)
  std::vector<int> newton_iterations_per_step;
  std::vector<double> time_step_size;

  // initialize liquid pressure, saturation errors
  double perr_linf_inf = 0.0, perr_linf_l1 = 0.0, perr_linf_l2 = 0.0, perr_l2_l2 = 0.0;
  double slerr_linf_inf = 0.0, slerr_linf_l1 = 0.0, slerr_linf_l2 = 0.0, slerr_l2_l2 = 0.0;
  
  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

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
      const auto& u5 = *S->Get<CompositeVector>("volumetric_flow_rate_liquid").ViewComponent("face");
      
      Epetra_MultiVector pl_ex(u0);
      Epetra_MultiVector sl_ex(u1);

      io->WriteVector(*u0(0), "liquid pressure", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u1(0), "saturation_liquid", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u2(0), "liquid hydrogen", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u3(0), "gas hydrogen", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u4(0), "gas pressure", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*u5(0), "volumetric flow rate liquid", AmanziMesh::Entity_kind::FACE);     

      exact_field(mesh, pl_ex, sl_ex, t);

      io->WriteVector(*pl_ex(0), "exact liquid pressure", AmanziMesh::Entity_kind::CELL);
      io->WriteVector(*sl_ex(0), "exact saturation_liquid", AmanziMesh::Entity_kind::CELL);

      io->FinalizeCycle();

      WriteStateStatistics(*S, *vo);
      //std::cout<<"Time = "<<t<<std::endl;
    }

    while (MPK->AdvanceStep(t, t + dt, false)) { dt /= 2.0; }
    
    /*
    auto& pl_c = *S->Get<CompositeVector>("pressure_liquid").ViewComponent("cell");
    for (int c = 0; c < ncells_owned; ++c) {
      const Amanzi::AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      double x = xc[0];
      pl_c[0][c] = x;
    }
    */

    MPK->CommitStep(t, t + dt, Tags::DEFAULT);

    auto& flux_l = *S->Get<CompositeVector>("volumetric_flow_rate_liquid").ViewComponent("face");
    auto& flux_g = *S->Get<CompositeVector>("volumetric_flow_rate_gas").ViewComponent("face");

    // store number of Newton iterations taken (only successful iterations after possible time step reduction) 
    double iter = MPK->bdf1_dae()->number_solver_iterations();
    newton_iterations_per_step.push_back(iter);
    // store time step size
    time_step_size.push_back(dt);
    
    S->advance_cycle();
    S->advance_time(dt);
      
    // update time
    t += dt;
    // reset time step if reduced
    if (dt < dt0) { std::cout<<"time step reduced to "<<dt<<std::endl; }
    dt = dt0;
    iloop++;

    // calculate error
    auto pl = *S->Get<CompositeVector>("pressure_liquid").ViewComponent("cell");
    auto sl = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell");

    double perr_linf_l1_tmp = 0.0, perr_linf_l2_tmp = 0.0, perr_l2_l2_tmp = 0.0;
    double slerr_linf_l1_tmp = 0.0, slerr_linf_l2_tmp = 0.0, slerr_l2_l2_tmp = 0.0;

    for (int c = 0; c < ncells_owned; ++c) { 
      const Amanzi::AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      double x = xc[0], y = xc[1];
      // exact solution expression
      double err = std::abs( pl[0][c] - x );
      double err_sl = std::abs( sl[0][c] - 1.0 );
      if (perr_linf_inf < err) { perr_linf_inf = err; }
      perr_linf_l1_tmp += err * mesh->getCellVolume(c);
      perr_linf_l2_tmp += err * err * mesh->getCellVolume(c);  
      perr_l2_l2_tmp += err * err * mesh->getCellVolume(c);
  
      if (slerr_linf_inf < err_sl) { slerr_linf_inf = err_sl; }
      slerr_linf_l1_tmp += err_sl * mesh->getCellVolume(c);
      slerr_linf_l2_tmp += err_sl * err_sl * mesh->getCellVolume(c);  
      slerr_l2_l2_tmp += err_sl * err_sl * mesh->getCellVolume(c);

    }
    perr_linf_l2_tmp = std::sqrt(perr_linf_l2_tmp);  
    perr_l2_l2 += perr_l2_l2_tmp * dt;
    if (perr_linf_l1 < perr_linf_l1_tmp) { perr_linf_l1 = perr_linf_l1_tmp; }
    if (perr_linf_l2 < perr_linf_l2_tmp) { perr_linf_l2 = perr_linf_l2_tmp; }


    slerr_linf_l2_tmp = std::sqrt(slerr_linf_l2_tmp);  
    slerr_l2_l2 += slerr_l2_l2_tmp * dt;
    if (slerr_linf_l1 < slerr_linf_l1_tmp) { slerr_linf_l1 = slerr_linf_l1_tmp; }
    if (slerr_linf_l2 < slerr_linf_l2_tmp) { slerr_linf_l2 = slerr_linf_l2_tmp; }


  }
  perr_l2_l2 = std::sqrt(perr_l2_l2);
  slerr_l2_l2 = std::sqrt(slerr_l2_l2);

  WriteStateStatistics(*S, *vo);

  // write iteration output to text file
  std::ofstream outFile("iterations_per_time_step.txt");
  for (int i = 0; i < newton_iterations_per_step.size(); ++i) { 
    outFile << newton_iterations_per_step[i] << "," << time_step_size[i] << std::endl;
  }
  outFile.close(); 

  // output errors
  std::cout<<"Pressure liquid error: "<<std::endl;
  std::cout<<"Linf(Linf) = "<<perr_linf_inf<<std::endl;
  std::cout<<"Linf(L1) = "<<perr_linf_l1<<std::endl;
  std::cout<<"Linf(L2) = "<<perr_linf_l2<<std::endl;
  std::cout<<"L2(L2) = "<<perr_l2_l2<<std::endl;
  std::cout<<" ---------------- "<<std::endl;
  std::cout<<"Saturation liquid error: "<<std::endl;
  std::cout<<"Linf(Linf) = "<<slerr_linf_inf<<std::endl;
  std::cout<<"Linf(L1) = "<<slerr_linf_l1<<std::endl;
  std::cout<<"Linf(L2) = "<<slerr_linf_l2<<std::endl;
  std::cout<<"L2(L2) = "<<slerr_l2_l2<<std::endl;
 
  // return L2(L2) error 
  return std::make_tuple(perr_l2_l2, slerr_l2_l2);
}


TEST(MULTIPHASE_1D_ANALYTICAL)
{
  std::vector<double> h(3), errs_pl(3), errs_sl(3);

  double dt0 = 1.0e-3;
  int i = 0;
  for (int M = 25; M <= 25; M *= 2) {
    h[i] = 1.0/M;
    auto errs = run_test(M, dt0, "test/multiphase_1D_analytical.xml");   
    errs_pl[i] = std::get<0>(errs);
    errs_sl[i] = std::get<1>(errs);
    dt0 /= 2.0;
    i += 1;
  }
  
  // calculate convergence rates for L2(L2) norm
  double rate1 = Amanzi::Utils::bestLSfit(h, errs_pl);
  double rate2 = Amanzi::Utils::bestLSfit(h, errs_sl);
  
  std::cout<<"Pressure: L2(L2) rate = "<<rate1<<std::endl;
  std::cout<<"Saturation: L2(L2) rate = "<<rate2<<std::endl;
  
  // verification
  CHECK(rate1 > 0.9 && rate2 > 0.9);  
}
