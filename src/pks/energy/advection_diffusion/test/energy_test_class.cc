#include "energy_test_class.hh"

EnergyTest::EnergyTest(Teuchos::RCP<Teuchos::ParameterList> plist_,
                             const Teuchos::RCP<AmanziMesh::Mesh>& mesh_,
                             int num_components_) :
    mesh(mesh_), parameter_list(plist_), num_components(num_components_) {

  // create states
  Teuchos::ParameterList state_plist =
    parameter_list->get<Teuchos::ParameterList>("State");
  S0 = Teuchos::rcp(new State(state_plist));
  S0->RegisterDomainMesh(mesh);

  Teuchos::RCP<Teuchos::ParameterList> energy_plist =
      Teuchos::sublist(parameter_list, "Energy");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

  // create the PK
<<<<<<< HEAD
  EPK = Teuchos::rcp(new Energy::AdvectionDiffusion(S0.ptr(), energy_plist, S0->FEList(), soln));
=======
  EPK = Teuchos::rcp(new Energy::AdvectionDiffusion(energy_plist, S0->FEList(), soln));
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
  EPK->setup(S0.ptr());
  S0->Setup();
}

void EnergyTest::initialize() {
  // initialize energy
  initialize_owned();
  initialize_mass_flux();
  EPK->initialize(S0.ptr());

  // initialize state, including darcy flux, sat, density, poro from parameter list
  S0->Initialize();
  S0->set_time(0.0);


  // finish checking state and create the state at the new timestep
  if (!S0->CheckAllFieldsInitialized()) {
    std::cout << "DID NOT INITIALIZE THINGS!" << std::endl;
  }
  S1 = Teuchos::rcp(new State(*S0));
  *S1 = *S0;
  EPK->set_states(S0,S0,S1);
}

void EnergyTest::commit_step(double dt) {
  EPK->commit_state(dt, S1);
  *S0 = *S1;
}

void EnergyTest::initialize_owned() {
  Teuchos::RCP<CompositeVector> temp = S0->GetFieldData("temperature", "energy");

  int c_owned = temp->size("cell", false);
  Epetra_MultiVector& temp_c = *temp->ViewComponent("cell",false);
  for (int c=0; c != c_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    for (int lcv_comp=0; lcv_comp != num_components; ++lcv_comp) {
      temp_c[lcv_comp][c] = my_f(xc, 0.0);
    }
  }

  Epetra_MultiVector& tc_f = *S0->GetFieldData("thermal_conductivity", "energy")
      ->ViewComponent("face",false);
  tc_f.PutScalar(my_K());

  if (temp->HasComponent("face")) {
    Epetra_MultiVector& temp_f = *temp->ViewComponent("face",false);
    for (int f=0; f != temp_f.MyLength(); ++f) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);
      for (int lcv_comp=0; lcv_comp != num_components; ++lcv_comp) {
        temp_f[lcv_comp][f] = my_f(xf, 0.0);
      }
    }
  }

  S0->GetField("temperature", "energy")->set_initialized();
  S0->GetField("thermal_conductivity", "energy")->set_initialized();
}

void EnergyTest::initialize_mass_flux() {
  const Epetra_BlockMap& fmap = mesh->face_map(true);
  Epetra_MultiVector& mass_flux = *S0->GetFieldData("mass_flux", "state")
      ->ViewComponent("face", true);

  for (int f=fmap.MinLID(); f<=fmap.MaxLID(); f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    const AmanziGeometry::Point& fc = mesh->face_centroid(f);
    mass_flux[0][f] = my_u(fc, 0.0) * normal;
  }
  S0->GetField("mass_flux", "state")->set_initialized();
}

void EnergyTest::evaluate_error_temp(double t, double* L1, double* L2) {
  const Epetra_BlockMap& cmap = mesh->cell_map(true);
  Teuchos::RCP<const CompositeVector> temp = S1->GetFieldData("temperature");

  double d;
  *L1 = 0.0;
  *L2 = 0.0;
  for (int c=cmap.MinLID(); c<=cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    for (int lcv_comp=0; lcv_comp != num_components; ++lcv_comp) {
      d = (*temp)("cell",lcv_comp,c) - my_f(xc, t);

      *L1 += fabs(d) * volume;
      *L2 += d * d * volume;
    }
  }
  *L2 = sqrt(*L2);
}

// // test problem with smooth solution
// EnergyTestSmooth::EnergyTestSmooth(Teuchos::RCP<Teuchos::ParameterList> plist_,
//         const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
//   EnergyTest::EnergyTest(plist_, mesh_, num_components_) {}

// double EnergyTestSmooth::my_f(const AmanziGeometry::Point& x, double t) {
//   return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
// }

// // test problem with cubic solution
// EnergyTestCubic::EnergyTestCubic(Teuchos::RCP<Teuchos::ParameterList> plist_,
//         const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
//   EnergyTest::EnergyTest(plist_, mesh_, num_components_) {}

// double EnergyTestCubic::my_f(const AmanziGeometry::Point& x, double t) {
//   if( x[0] < 1 + t ) return 1;
//   if( x[0] > 3 + t ) return 0;
//   double z = (x[0]-1-t) / 2;
//   return 2*z*z*z - 3*z*z + 1;
// }

// test problem on 2D square with velocity in the <1,1> x-y direction, still virtual
AmanziGeometry::Point EnergyTestTwoDOne::my_u(const AmanziGeometry::Point& x, double t) {
  return AmanziGeometry::Point(1.0, 1.0);
}

double EnergyTestTwoDOne::my_f(const AmanziGeometry::Point& x, double t) {
  return 1.0;
}


