#include "transport_test_class.hh"

TransportTest::TransportTest(Teuchos::ParameterList& plist_,
                             const Teuchos::RCP<AmanziMesh::Mesh>& mesh_,
                             int num_components_) :
    mesh(mesh_), parameter_list(plist_), num_components(num_components_) {

  // create states
  Teuchos::ParameterList state_plist =
    parameter_list.get<Teuchos::ParameterList>("State");
  S0 = Teuchos::rcp(new State(state_plist, mesh));
  Teuchos::ParameterList transport_plist =
    parameter_list.get<Teuchos::ParameterList>("Transport");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector("solution"));

  // create the PK
  TPK = Teuchos::rcp(new Transport::PassiveTracer(transport_plist, S0, soln));
}

void TransportTest::initialize() {
  // initialize state, including darcy flux, sat, density, poro from parameter list
  S0->Initialize();

  // initialize transport
  TPK->initialize(S0);
  initialize_tcc();
  initialize_darcy_flux();

  // finish checking state and create the state at the new timestep
  if (!S0->CheckAllInitialized()) {
    std::cout << "DID NOT INITIALIZE THINGS!" << std::endl;
  }
  S1 = Teuchos::rcp(new State(*S0));
  *S1 = *S0;
  TPK->set_states(S0,S0,S1);
}

void TransportTest::commit_step() {
  *S0 = *S1;
}

void TransportTest::initialize_tcc() {
  const Epetra_BlockMap& cmap = mesh->cell_map(true);
  Teuchos::RCP<CompositeVector> tcc =
    S0->GetFieldData("concentration", "transport");

  for (int c=cmap.MinLID(); c<=cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    for (int lcv_comp=0; lcv_comp != num_components; ++lcv_comp) {
      (*tcc)(lcv_comp,c) = my_f(xc, 0.0);
    }
  }
  S0->GetField("concentration", "transport")->set_initialized();
}

void TransportTest::initialize_darcy_flux() {
  const Epetra_BlockMap& fmap = mesh->face_map(true);
  Teuchos::RCP<CompositeVector> darcy_flux =
    S0->GetFieldData("darcy_flux", "state");

  for (int f=fmap.MinLID(); f<=fmap.MaxLID(); f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    const AmanziGeometry::Point& fc = mesh->face_centroid(f);
    (*darcy_flux)(f) = my_u(fc, 0.0) * normal;
  }
  S0->GetField("darcy_flux", "state")->set_initialized();
}

void TransportTest::evaluate_error_tcc(double t, double* L1, double* L2) {
  const Epetra_BlockMap& cmap = mesh->cell_map(true);
  Teuchos::RCP<const CompositeVector> tcc =
    S1->GetFieldData("concentration");

  double d;
  *L1 = 0.0;
  *L2 = 0.0;
  for (int c=cmap.MinLID(); c<=cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    for (int lcv_comp=0; lcv_comp != num_components; ++lcv_comp) {
      d = (*tcc)(lcv_comp,c) - my_f(xc, t);

      *L1 += fabs(d) * volume;
      *L2 += d * d * volume;
    }
  }
  *L2 = sqrt(*L2);
}

AmanziGeometry::Point TransportTest::my_u(const AmanziGeometry::Point& x, double t) {
  return AmanziGeometry::Point(1.0, 0.0, 0.0);
}

// test problem with constant solution of 1
TransportTestOne::TransportTestOne(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  TransportTest::TransportTest(plist_, mesh_, num_components_) {}

double TransportTestOne::my_f(const AmanziGeometry::Point& x, double t) {
  return 1;
}

// test problem with step solution
TransportTestStep::TransportTestStep(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  TransportTest::TransportTest(plist_, mesh_, num_components_) {}

double TransportTestStep::my_f(const AmanziGeometry::Point& x, double t) {
  if (x[0] <= 1 + t) return 1;
  return 0;
}

// test problem with smooth solution
TransportTestSmooth::TransportTestSmooth(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  TransportTest::TransportTest(plist_, mesh_, num_components_) {}

double TransportTestSmooth::my_f(const AmanziGeometry::Point& x, double t) {
  return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
}

// test problem with cubic solution
TransportTestCubic::TransportTestCubic(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  TransportTest::TransportTest(plist_, mesh_, num_components_) {}

double TransportTestCubic::my_f(const AmanziGeometry::Point& x, double t) {
  if( x[0] < 1 + t ) return 1;
  if( x[0] > 3 + t ) return 0;
  double z = (x[0]-1-t) / 2;
  return 2*z*z*z - 3*z*z + 1;
}

// test problem on 2D square with velocity in the <1,1> x-y direction, still virtual
TransportTestTwoDOne::TransportTestTwoDOne(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  TransportTest::TransportTest(plist_, mesh_, num_components_) {}

AmanziGeometry::Point TransportTestTwoDOne::my_u(const AmanziGeometry::Point& x, double t) {
  return AmanziGeometry::Point(1.0, 1.0);
}

double TransportTestTwoDOne::my_f(const AmanziGeometry::Point& x, double t) {
  return 1;
}


