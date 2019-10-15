/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "PK.hh"
#include "State.hh"
#include "Teuchos_RCP.hpp"

using namespace Amanzi;
//
// Helper struct to store a run for the below function
// ================================================================================
struct Run {
  Run(const Teuchos::RCP<State>& S_, const Teuchos::RCP<PK>& pk_)
    : S(S_), pk(pk_)
  {}

  Teuchos::RCP<State> S;
  Teuchos::RCP<PK> pk;
};

//
// Helper function that mocks a coordinator to run the test.
// ================================================================================
std::pair<double, double>
run_test(const Teuchos::RCP<State>& S, const Teuchos::RCP<PK>& pk,
         double end_time = 1.0)
{
  pk->Setup();
  S->Setup();

  pk->Initialize();
  S->Initialize();

  int nsteps_good = 0;
  int nsteps_bad = 0;

  double t_final = end_time;
  bool done = false;
  while (!done) {
    double dt = std::min(pk->get_dt(), 1.0);
    S->set_time("next", S->time() + dt);
    S->set_cycle("next", S->cycle() + 1);

    bool fail = pk->AdvanceStep("", "next");

    if (fail) {
      pk->FailStep("", "next");
      nsteps_bad++;
    } else {
      pk->CommitStep("", "next");
      S->set_time("", S->time("next"));
      S->set_cycle("", S->cycle("next"));
      nsteps_good++;
    }

    done = S->time("") >= t_final - 0.0001;
  }
  return std::make_pair(nsteps_good, nsteps_bad);
}

//
// Helper function for creating a run with one MPC and two PKs
// ============================================================================
template <class MPC_t, class PK_A_t, class PK_B_t, class PK_t = Amanzi::PK>
std::unique_ptr<Run>
createRunMPC(const std::string& mpc_name, const std::string& pk_A_name,
             const std::string& pk_B_name,
             const std::string& filename = "test/pks_ode.xml")
{
  std::cout << "Test: " << mpc_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile(filename);
  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  auto comm = getDefaultComm();

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  // make a 1x1 'mesh' for ODEs
  auto mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  S->RegisterDomainMesh(mesh);

  // create the PKs/MPCs
  auto pk_tree_mpc = Teuchos::rcp(new Teuchos::ParameterList(mpc_name));
  auto mpc = Teuchos::rcp(new MPC_t(pk_tree_mpc, global_list, S));

  auto pk_tree_A = Teuchos::rcp(new Teuchos::ParameterList(pk_A_name));
  auto pk_a = Teuchos::rcp(new PK_A_t(pk_tree_A, global_list, S));

  auto pk_tree_B = Teuchos::rcp(new Teuchos::ParameterList(pk_B_name));
  auto pk_b = Teuchos::rcp(new PK_B_t(pk_tree_B, global_list, S));
  mpc->SetChildren(std::vector<Teuchos::RCP<PK_t>>{ pk_a, pk_b });

  return std::make_unique<Run>(S, mpc);
}

//
// Helper function for creating a run with one PK, an ODE
// This makes a 0D mesh (single cell).
// ============================================================================
template <class PK_t>
std::unique_ptr<Run>
createRunODE(const std::string& pk_name,
             const std::string& filename = "test/pks_ode.xml")
{
  std::cout << "Test: " << pk_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile(filename);
  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  auto comm = getDefaultComm();

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  // make a 1x1 'mesh' for ODEs
  auto mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  S->RegisterDomainMesh(mesh);

  // create the PKs/MPCs
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));
  auto pk = Teuchos::rcp(new PK_t(pk_tree, global_list, S));

  return std::make_unique<Run>(S, pk);
}


//
// Helper function for creating a run with PK, a PDE
// Note that this is identical to ODE above but includes a 1D mesh.
// ============================================================================
template <class PK_t>
std::unique_ptr<Run>
createRunPDE(const std::string& pk_name, const std::string& filename)
{
  std::cout << "Test: " << pk_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile(filename);
  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  auto comm = getDefaultComm();

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  // make a 1x1 'mesh' for ODEs
  auto mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, 100, 1);
  S->RegisterDomainMesh(mesh);

  // create the PKs/MPCs
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));
  auto pk = Teuchos::rcp(new PK_t(pk_tree, global_list, S));

  return std::make_unique<Run>(S, pk);
}
