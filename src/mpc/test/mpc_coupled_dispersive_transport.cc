#include <iostream>
#include "stdlib.h"
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
#include "MeshAudit.hh"
#include "eos_registration.hh"
#include "mdm_transport_registration.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "mpc_pks_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_transport_registration.hh"
#include "State.hh"


double RunTest(int icase, double u_f,
               double mol_diff_f, double mol_diff_m,
               double normal_diff) {
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  
  std::string xmlInFileName = "test/mpc_coupled_dispersive_transport.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // modify the test
  plist->sublist("state").sublist("initial conditions")
        .sublist("fracture-volumetric_flow_rate").sublist("function")
        .sublist("All domain").sublist("function").sublist("dof 1 function")
        .sublist("function-constant").set<double>("value", u_f);

  std::vector<double> tmp_f({ mol_diff_f });
  plist->sublist("PKs").sublist("transport fracture").sublist("molecular diffusion")
        .set<Teuchos::Array<double> >("aqueous values", tmp_f);

  std::vector<double> tmp_m({ mol_diff_m });
  plist->sublist("PKs").sublist("transport matrix").sublist("molecular diffusion")
        .set<Teuchos::Array<double> >("aqueous values", tmp_m);

  plist->sublist("state").sublist("initial conditions")
        .sublist("fracture-normal_diffusion").sublist("function").sublist("All").sublist("function")
        .sublist("function-constant").set<double>("value", normal_diff);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({Framework::MSTK}));
  auto mesh = factory.create(0.0, 0.0, 0.0, 12.0, 1.0, 6.0, 96, 1, 8);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;    
  
  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  auto mesh_fracture = factory.create(mesh, names, AmanziMesh::FACE);

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  const auto& tcc_m = *S->Get<CompositeVector>("total_component_concentration").ViewComponent("cell");
  double cmin(1e+99), cmax(-1e+99);
  for (int c = 0; c < tcc_m.MyLength(); ++c) {
    const auto& xc = mesh->cell_centroid(c);
    if (std::fabs(xc[0] - 4.0625) < 1e-3) {
      cmin = std::min(cmin, tcc_m[0][c]);
      cmax = std::max(cmax, tcc_m[0][c]);
    }
    CHECK(tcc_m[0][c] < 1.0);
  }

  const auto& tcc_f = *S->Get<CompositeVector>("fracture-total_component_concentration").ViewComponent("cell");
  double fmin(1e+99), fmax(-1e+99), err(0.0);

  double b = (*S->Get<CompositeVector>("fracture-aperture").ViewComponent("cell"))[0][0];
  double kn = (*S->Get<CompositeVector>("fracture-normal_diffusion").ViewComponent("cell"))[0][0];

  double t(1e+5), Df(mol_diff_f / b), Dm(mol_diff_m);
  CHECK_CLOSE(kn, Dm / 2, 1e-12);

  for (int c = 0; c < tcc_f.MyLength(); ++c) {
    const auto& xc = mesh_fracture->cell_centroid(c);
    if (std::fabs(xc[0] - 4.0625) < 1e-3) {
      fmin = std::min(fmin, tcc_f[0][c]);
      fmax = std::max(fmax, tcc_f[0][c]);
    }
    if (icase == 1) {
      err += fabs(tcc_f[0][c] - std::erfc(xc[0] / std::sqrt(Df * t) / 2));
    } else if (icase == 2) {
      double Tm = u_f * t / b - xc[0];
      if (Tm > 0.0) err += fabs(tcc_f[0][c] - std::erfc(xc[0] * kn * std::sqrt(b / (u_f * Dm * Tm)) / 2));
    } else{
      std::cout << xc[0] << " " << tcc_f[0][c] << std::endl;
    }
  }
  double err_tmp(err);
  mesh->get_comm()->SumAll(&err_tmp, &err, 1);
  err /= tcc_f.GlobalLength();
  std::cout << "Mean error in fracture: " << err << std::endl;

  if (fmin < 1.0e+98) CHECK_CLOSE(fmin, fmax, 1e-12);

  return err;
}


// double Integrate(const WhetStone::Function& fun) {
// }


TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_0) {
  // no coupling, only diffusion
  double err = RunTest(1, 0.0, 1e-6, 0.0, 0.0);
  CHECK(err < 3.0e-3);
}

TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_1) {
  // d(a C)/dt + d(u C)/dx - (Dm/2) (Cm - C) = 0
  // kn = Dm/2 for the analytic solution
  double err = RunTest(2, 1.0e-5, 0.0, 1.0e-7, 1.0e-7/2);
  CHECK(err < 2.0e-2);
}

/*
TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_2) {
  // d(a C)/dt + d(u C)/dx - a Df d^2(C)/dx^2 - (Dm/2) (Cm - C) = 0
  double err = RunTest(3, 1.0e-5, 1.0e-5, 1.0e-6, 1.0e-6/2);
  CHECK(err < 100);
}
*/


