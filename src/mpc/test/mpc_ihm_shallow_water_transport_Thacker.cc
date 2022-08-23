#include "stdlib.h"
#include <iostream>
#include <tuple>
#include "math.h"

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "eos_registration.hh"
#include "IO.hh"
#include "LeastSquare.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "mpc_pks_registration.hh"
#include "numerical_flux_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_transport_registration.hh"
#include "pks_shallow_water_registration.hh"
#include "State.hh"

std::tuple<double, double, double> RunTest(int n, int* ncycles, int* MyPID) {
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  *MyPID = comm->MyPID();
  
  std::string xmlInFileName = "test/mpc_ihm_shallow_water_transport_Thacker.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));
  
  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({Framework::MSTK}));
  auto mesh = factory.create(-3.0, -3.0, 3.0, 3.0, n, n, true, true);
  
  Amanzi::ObservationData obs_data;
  
  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);
  
  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();
  
  // error calculations
  *ncycles = S->Get<int>("cycle", Tags::DEFAULT);
  double err_h, err_v, err_tcc;
  auto h = S->Get<CompositeVector>("ponded_depth", Tags::DEFAULT);
  auto v = S->Get<CompositeVector>("velocity", Tags::DEFAULT);
  auto tcc = S->Get<CompositeVector>("total_component_concentration", Tags::DEFAULT);

  Teuchos::ParameterList sublist;
  CompositeVector h_ex(h), v_ex(v), tcc_ex(tcc);

  // -- height
  std::vector<std::string> subfieldnames = { "cell" };
  sublist = plist->sublist("state").sublist("initial conditions").sublist("ponded_depth");
  sublist.set<double>("time", 1.0);
  Helpers::Initialize(sublist, h_ex, "ponded_depth", &subfieldnames);

  h_ex.Update(-1.0, h, 1.0);
  h_ex.Norm1(&err_h);

  // -- velocity
  sublist = plist->sublist("state").sublist("initial conditions").sublist("velocity");
  sublist.set<double>("time", 1.0);
  Helpers::Initialize(sublist, v_ex, "velocity", &subfieldnames);

  v_ex.Update(-1.0, v, 1.0);
  v_ex.Norm1(&err_v);

  // -- concentration
  sublist = plist->sublist("state").sublist("initial conditions").sublist("total_component_concentration");
  sublist.set<double>("time", 1.0);
  Helpers::Initialize(sublist, tcc_ex, "total_component_conventration", &subfieldnames);

  tcc_ex.Update(-1.0, tcc, 1.0);
  tcc_ex.Norm1(&err_tcc);

  int ncells = h.GlobalLength();
  return std::make_tuple(err_h / ncells, err_v / ncells, err_tcc / ncells);
}


TEST(MPC_DRIVER_IHM_SHALLOW_WATER_TRANSPORT_THACKER) {
  int i(0), ncycles, MyPID;
  std::vector<double> h(3), err_h(3), err_v(3), err_tcc(3);
  for (int n = 16; n < 80; n *= 2, ++i) {
    auto errs = RunTest(n, &ncycles, &MyPID);
    h[i] = 1.0 / n;
    err_h[i] = std::get<0>(errs);
    err_v[i] = std::get<1>(errs);
    err_tcc[i] = std::get<2>(errs);
    if (MyPID == 0) {
      std::cout << "Error: h=" << err_h[i] << " u=" << err_v[i] << " tcc=" << err_tcc[i]
                << ",  dx=" << h[i] << " dt=" << 1.0 / ncycles << std::endl;
    }
  }

  double rate1 = Amanzi::Utils::bestLSfit(h, err_h);
  double rate2 = Amanzi::Utils::bestLSfit(h, err_v);
  double rate3 = Amanzi::Utils::bestLSfit(h, err_tcc);
  if (MyPID == 0) {
    std::cout << "Error convergence rates: " << rate1 << " " << rate2 << " " << rate3 << std::endl;
  }
  CHECK(rate1 > 1.9 && rate2 > 1.9 && rate3 > 1.8);
}

