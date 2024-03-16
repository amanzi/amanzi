/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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
#include "bilinear_form_reg.hh"
#include "CycleDriver.hh"
#include "MeshAudit.hh"
#include "eos_reg.hh"
#include "evaluators_flow_reg.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "models_flow_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_reg.hh"
#include "pks_mechanics_reg.hh"
#include "pks_mpc_reg.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

void
RunTest(const std::string xmlInFileName)
{
  Comm_ptr_type comm = Amanzi::getDefaultComm();

  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  int nx(40), ny(6);
  double L(10.0);
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  mesh_list->set<bool>("request faces", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, L, 3.0, nx, ny);

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  Amanzi::ObservationData obs_data;
  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  S = cycle_driver.Go();

  // compute analytic solution
  double g(10.0), p0(1e+7), rho(1e+3), mu(1e-3), k(6e-18), phi_cp(1e-11);
  double E(3e+10), nu(0.0), t(36000.0);
  double Ss = rho * g * (phi_cp + 1.0 / E);
  double cv = rho * g / mu * k / Ss;  // effective consolidation coefficient

  std::vector<double> pex(nx, 0.0), eex(nx, 0.0);
  for (int n = 0; n < 100; ++n) {
    int k = 2 * n + 1;
    double factor = M_PI / 2 / L;
    for (int c = 0; c < nx; ++c) {
      double x = (mesh->getCellCentroid(ny * c))[0];
      pex[c] += (1.0 / k) * std::sin(k * factor * x) * std::exp(-k * k * factor * factor * cv * t);
    }
  }
  for (int c = 0; c < nx; ++c) {
    pex[c] *= 4 * p0 / M_PI;
    eex[c] = (pex[c] - p0) * (1 + nu) / E;
  }

  // compute error
  double perr(0.0), eerr(0.0), pnorm, enorm; 
  const auto& p = *S->Get<CompositeVector>("pressure", Tags::DEFAULT).ViewComponent("cell");
  const auto& e = *S->Get<CompositeVector>("volumetric_strain", Tags::DEFAULT).ViewComponent("cell");
  p.Norm2(&pnorm);
  e.Norm2(&enorm);

  for (int c = 0; c < nx; ++c) {
    perr += std::pow(pex[c] - p[0][ny * c], 2);
    eerr += std::pow(eex[c] - e[0][ny * c], 2);
    // std::cout << c << " " << pex[c] << " " << p[0][ny * c] << " " << eex[c] << " " << e[0][ny * c] << std::endl;
  }
  perr = std::sqrt(perr) / pnorm;
  eerr = std::sqrt(eerr) / enorm;
  CHECK(perr < 0.0025);
  CHECK(eerr < 0.005);
}


TEST(MPC_DRIVER_BIOT_CONSOLIDATION)
{
  RunTest("test/mpc_biot_consolidation.xml");
}

TEST(MPC_DRIVER_BIOT_CONSOLIDATION_DARCY)
{
  RunTest("test/mpc_biot_consolidation_darcy.xml");
}
