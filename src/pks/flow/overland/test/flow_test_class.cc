#include "flow_test_class.hh"

FlowTest::FlowTest(Teuchos::ParameterList& plist_,
		   const Teuchos::RCP<AmanziMesh::Mesh>& mesh_,
		   int num_components_) :
  mesh(mesh_), parameter_list(plist_), num_components(num_components_) {

  // create states
  Teuchos::ParameterList state_plist =
    parameter_list.get<Teuchos::ParameterList>("State");
  S0 = Teuchos::rcp(new State(state_plist));
  S0->RegisterDomainMesh(mesh);
  Teuchos::ParameterList flow_plist =
    parameter_list.get<Teuchos::ParameterList>("Flow");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector("solution"));

  // create the PK
  FPK = Teuchos::rcp(new Flow::OverlandFlow(flow_plist, S0, soln));
}

void FlowTest::initialize() {
  // initialize state, including darcy flux, sat, density, poro from parameter list
  S0->Initialize();
  S0->set_time(0.0);

  // initialize flow
  FPK->initialize(S0);
  FPK->commit_state(0.0, S0);
  //  initialize_owned();

  // finish checking state and create the state at the new timestep
  if (!S0->CheckAllInitialized()) {
    std::cout << "DID NOT INITIALIZE THINGS!" << std::endl;
  }
  S1 = Teuchos::rcp(new State(*S0));
  *S1 = *S0;
  FPK->set_states(S0,S0,S1);
}

void FlowTest::commit_step() {
  *S0 = *S1;
}

void FlowTest::initialize_owned() {
  Teuchos::RCP<CompositeVector> pres = S0->GetFieldData("overland_pressure", "overland_flow");

  int c_owned = pres->size("cell");
  for (int c=0; c != c_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*pres)("cell",c) = my_f(xc, 0.0);
  }
  S0->GetRecord("overland_pressure", "overland_flow")->set_initialized();
}


void FlowTest::evaluate_error_pressure(double t, double & L1, double & L2) {
  const Epetra_BlockMap & cmap = mesh->cell_map(true);
  const CompositeVector & pres = *(S1->GetFieldData("overland_pressure"));

  double d;
  L1 = 0.0;
  L2 = 0.0;
  int c_owned = pres.size("cell");
  for (int c=0; c!=c_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    d = pres("cell",c) - my_f(xc, t);

    L1 += fabs(d) * volume;
    L2 += d * d * volume;
  }
  L2 = sqrt(L2);
}

// test problem with constant solution of 1
FlowTestOne::FlowTestOne(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  FlowTest::FlowTest(plist_, mesh_, num_components_) {}

double FlowTestOne::my_f(const AmanziGeometry::Point& x, double t) {
  return 101325.0;
}

// test problem on 2D square with velocity in the <1,1> x-y direction, still virtual
FlowTestTwoDOne::FlowTestTwoDOne(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  FlowTest::FlowTest(plist_, mesh_, num_components_) {}

double FlowTestTwoDOne::my_f(const AmanziGeometry::Point& x, double t) {
  return 101325.0;
}
