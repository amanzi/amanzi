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

#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "models_energy_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_energy_reg.hh"
#include "pks_flow_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"
#include "models_flow_reg.hh"


std::tuple<double, Amanzi::CompositeVector, Amanzi::CompositeVector>
RunTest(const std::string& model)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  auto plist = Teuchos::getParametersFromXmlFile("test/mpc_thermal_richards_highPT.xml");
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");

  if (model == "FEHM") {
    plist->sublist("state").sublist("evaluators").sublist("internal_energy_liquid")
      .sublist("IEM parameters").sublist("EntireDomain").sublist("IEM parameters")
      .set<std::string>("iem type", "linear");

    plist->sublist("state").sublist("evaluators").sublist("molar_density_liquid")
      .sublist("EOS parameters").set<std::string>("eos type", "liquid water FEHM");

    plist->sublist("state").sublist("evaluators").sublist("viscosity_liquid")
      .sublist("EOS parameters").set<std::string>("eos type", "liquid water FEHM");

    plist->sublist("PKs").sublist("transient:energy").sublist("thermal conductivity evaluator")
      .sublist("TCM_0").sublist("liquid phase")
      .set<std::string>("eos type", "constant");
  }

  auto comm = Amanzi::getDefaultComm();
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 100.0, 20.0, 50, 4);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  double time = S->Get<double>("time");
  return {time, S->Get<CompositeVector>("pressure"), S->Get<CompositeVector>("temperature")};
}


TEST(MPC_DRIVER_THERMAL_RICHARDS_HIGH_PT)
{
  auto [time1, p1, t1] = RunTest("FEHM");
  auto [time2, p2, t2] = RunTest("IAPW97");

  double perr, terr, pnorm, tnorm;
  p2.Update(-1.0, p1, 1.0);
  t2.Update(-1.0, t1, 1.0);

  CHECK(std::abs(time1 - time2) < 1.0e-12 * time1);
  
  for (auto comp = p1.begin(); comp != p1.end(); ++comp) {
    auto& p1v = *p1.ViewComponent(*comp);
    auto& p2v = *p2.ViewComponent(*comp);
    p1v.Norm2(&pnorm);
    p2v.Norm2(&perr);
    perr /= pnorm;

    auto& t1v = *t1.ViewComponent(*comp);
    auto& t2v = *t2.ViewComponent(*comp);
    t1v.Norm2(&tnorm);
    t2v.Norm2(&terr);
    terr /= tnorm;

    printf("\nP-error =%10.5f  norm =%12.5g  comp=%s\n", perr, pnorm, comp->c_str());
    printf("T-error =%10.5f  norm =%12.5g\n", terr, tnorm);
    CHECK(perr < 0.005);
    CHECK(terr < 0.002);
  }
}
