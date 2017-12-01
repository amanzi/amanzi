#include "Teuchos_RCP.hpp"
#include "State.hh"
#include "PK.hh"


//
// Helper struct to store a run for the below function
// ================================================================================
struct Run {
  Run(const Teuchos::RCP<State>& S_,
      const Teuchos::RCP<PK>& pk_)
      : S(S_),
        pk(pk_) {}

  Teuchos::RCP<State> S;
  Teuchos::RCP<PK> pk;
};


//
// Helper function that mocks a coordinator to run the test.
// ================================================================================
std::pair<double,double> run_test(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<PK>& pk) {

  pk->Setup();
  S->Setup();

  pk->Initialize();
  S->Initialize();

  int nsteps_good = 0;
  int nsteps_bad = 0;

  double t_final = 1.0;
  bool done = false;
  while (!done) {
    double dt = std::min(pk->get_dt(), 1.0);
    S->set_time("next", S->time()+dt);
    S->set_cycle("next", S->cycle()+1);

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

    done = S->time("") >= t_final-0.0001;
  }
  return std::make_pair(nsteps_good, nsteps_bad);
}


//
// Helper function for creating a run with one MPC and two PKs
// ============================================================================
template<class MPC_t, class PK_A_t, class PK_B_t, class PK_t=Amanzi::PK>
std::unique_ptr<Run>
createRunMPC(const std::string& mpc_name, const std::string& pk_A_name, const std::string& pk_B_name) {
  std::cout << "Test: " << mpc_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, comm));

  Amanzi::AmanziMesh::FrameworkPreference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::MSTK);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  // make a 1x1 'mesh' for ODEs
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1, gm);
  S->RegisterDomainMesh(mesh);

  // create the PKs/MPCs
  auto pk_tree_mpc = Teuchos::rcp(new Teuchos::ParameterList(mpc_name));
  auto mpc = Teuchos::rcp(new MPC_t(pk_tree_mpc, global_list, S));

  auto pk_tree_A = Teuchos::rcp(new Teuchos::ParameterList(pk_A_name));
  auto pk_a = Teuchos::rcp(new PK_A_t(pk_tree_A, global_list, S));

  auto pk_tree_B = Teuchos::rcp(new Teuchos::ParameterList(pk_B_name));
  auto pk_b = Teuchos::rcp(new PK_B_t(pk_tree_B, global_list, S));
  mpc->SetChildren(std::vector<Teuchos::RCP<PK_t> >{pk_a, pk_b});

  return std::make_unique<Run>(S,mpc);
}






