/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy PK
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
#include "UnitTest++.h"

// Amanzi
#include "CompositeVector.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "IAPWS97.hh"
#include "MeshFactory.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "ThermodynamicStateEvaluators.hh"
#include "TotalEnergyEvaluatorPH.hh"
#include "VerboseObject.hh"

TEST(EVALUATORS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Energy;
  using namespace Amanzi::Evaluators;

  using CV_t = CompositeVector;
  using CVS_t = CompositeVectorSpace;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  std::string domain("domain");

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 17, 17);

  Teuchos::ParameterList slist;
  slist.sublist("verbose object").set<std::string>("verbosity level", "extreme");
  Teuchos::RCP<State> S = Teuchos::rcp(new Amanzi::State(slist));
  S->RegisterDomainMesh(mesh);

  Key state_key = Keys::getKey(domain, "thermodynamic_state");
  Key pressure_key = Keys::getKey(domain, "pressure");
  Key enthalpy_key = Keys::getKey(domain, "enthalpy");
  Key mass_density_key = Keys::getKey(domain, "mass_density_liquid");
  Key molar_density_key = Keys::getKey(domain, "molar_density_liquid");
  Key temperature_key = Keys::getKey(domain, "temperature");
  Key conductivity_key = Keys::getKey(domain, "thermal_conductivity");

  // thermodynamics
  S->Require<CV_t, CVS_t>(state_key, Tags::DEFAULT, state_key)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, TS_t_size);

  Teuchos::ParameterList elist1(state_key);
  elist1.set<std::string>("tag", "");
  elist1.sublist("verbose object").set<std::string>("verbosity level", "extreme");

  auto eval1 = Teuchos::rcp(new ThermodynamicStateEvaluator(elist1));
  S->SetEvaluator(state_key, Tags::DEFAULT, eval1);

  // primary variables
  S->Require<CV_t, CVS_t>(pressure_key, Tags::DEFAULT)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  AddDefaultPrimaryEvaluator(S, pressure_key);
 
  S->Require<CV_t, CVS_t>(enthalpy_key, Tags::DEFAULT)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  AddDefaultPrimaryEvaluator(S, enthalpy_key);

  // densities
  S->Require<CV_t, CVS_t>(mass_density_key, Tags::DEFAULT, mass_density_key)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S->Require<CV_t, CVS_t>(molar_density_key, Tags::DEFAULT, molar_density_key)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  Teuchos::ParameterList elist2(mass_density_key);
  elist2.set<std::string>("tag", "");
  elist2.sublist("verbose object").set<std::string>("verbosity level", "extreme");

  auto eval2 = Teuchos::rcp(new DensityEvaluator(elist2));
  S->SetEvaluator(mass_density_key, Tags::DEFAULT, eval2);

  // auto eval2b = Teuchos::rcp(new DensityEvaluator(elist2));
  // S->SetEvaluator(molar_density_key, Tags::DEFAULT, eval2b);

  S->RequireDerivative<CV_t, CVS_t>(
    mass_density_key, Tags::DEFAULT, pressure_key, Tags::DEFAULT, mass_density_key)
    .SetGhosted();
  S->RequireDerivative<CV_t, CVS_t>(
    mass_density_key, Tags::DEFAULT, enthalpy_key, Tags::DEFAULT, mass_density_key)
    .SetGhosted();

  // temperature
  S->Require<CV_t, CVS_t>(temperature_key, Tags::DEFAULT, temperature_key)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S->RequireDerivative<CV_t, CVS_t>(
    temperature_key, Tags::DEFAULT, pressure_key, Tags::DEFAULT, temperature_key)
    .SetGhosted();
  S->RequireDerivative<CV_t, CVS_t>(
    temperature_key, Tags::DEFAULT, enthalpy_key, Tags::DEFAULT, temperature_key)
    .SetGhosted();

  Teuchos::ParameterList elist3(temperature_key);
  elist3.set<std::string>("tag", "");
  elist3.sublist("verbose object").set<std::string>("verbosity level", "extreme");

  auto eval3 = Teuchos::rcp(new TemperatureEvaluator(elist3));
  S->SetEvaluator(temperature_key, Tags::DEFAULT, eval3);

  // thermal conductivity
  S->Require<CV_t, CVS_t>(conductivity_key, Tags::DEFAULT, conductivity_key)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S->RequireDerivative<CV_t, CVS_t>(
    conductivity_key, Tags::DEFAULT, pressure_key, Tags::DEFAULT, conductivity_key)
    .SetGhosted();
  S->RequireDerivative<CV_t, CVS_t>(
    conductivity_key, Tags::DEFAULT, enthalpy_key, Tags::DEFAULT, conductivity_key)
    .SetGhosted();

  Teuchos::ParameterList elist4(conductivity_key);
  elist4.set<std::string>("tag", "");
  elist4.sublist("verbose object").set<std::string>("verbosity level", "extreme");

  auto eval4 = Teuchos::rcp(new ThermalConductivityEvaluator(elist4));
  S->SetEvaluator(conductivity_key, Tags::DEFAULT, eval4);

  S->Setup();

  // initial value, replacement for PK initialization
  std::string passwd("");
  auto& p_c = *S->GetW<CompositeVector>(pressure_key, passwd).ViewComponent("cell");
  auto& h_c = *S->GetW<CompositeVector>(enthalpy_key, passwd).ViewComponent("cell");

  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  for (int c = 0; c < ncells; ++c) {
    p_c[0][c] = 31.0e+6 * (1.0 - 0.5 * c / ncells);
    h_c[0][c] = 1000.0 * 18.015 * (1.0 + 2.0 * c / ncells);
  }

  S->InitializeFields();
  S->InitializeEvaluators();

  Tag tag = Tags::DEFAULT;
  auto eval_p = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(S->GetEvaluatorPtr(pressure_key, tag));
  auto eval_h = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(S->GetEvaluatorPtr(enthalpy_key, tag));

  // test 1
  std::cout << "\n========== Test 1a (state + variable update) ==========\n ";
  eval_p->SetChanged();
  S->GetEvaluator(mass_density_key).Update(*S, "test");

  std::cout << "\n========== Test 1b (no update) ==========\n ";
  S->GetEvaluator(temperature_key).Update(*S, "test");

  // test 2
  std::cout << "\n========== Test 2a (state + variable update) ==========\n ";
  eval_h->SetChanged();
  S->GetEvaluator(mass_density_key).Update(*S, "test");

  std::cout << "\n========== Test 2b (no state update) ==========\n ";
  S->GetEvaluator(temperature_key).Update(*S, "test");

  // test 3
  std::cout << "\n========== Test 3a (state + variable + derivative update) ==========\n ";
  eval_p->SetChanged();
  S->GetEvaluator(mass_density_key).UpdateDerivative(*S, "test", pressure_key, Tags::DEFAULT);

  std::cout << "\n========== Test 3b (state + variable + derivative update) ==========\n ";
  eval_h->SetChanged();
  S->GetEvaluator(molar_density_key).UpdateDerivative(*S, "test", enthalpy_key, Tags::DEFAULT);

  std::cout << "\n========== Test 3c (variable + derivative update) ==========\n ";
  S->GetEvaluator(temperature_key).UpdateDerivative(*S, "test", pressure_key, Tags::DEFAULT);

  std::cout << "\n========== Test 3d (derivative update) ==========\n ";
  S->GetEvaluator(temperature_key).UpdateDerivative(*S, "test", enthalpy_key, Tags::DEFAULT);

  // test 4
  std::cout << "\n========== Test 4 (state + variable + derivative update) ==========\n ";
  eval_h->SetChanged();
  S->GetEvaluator(conductivity_key).UpdateDerivative(*S, "test", pressure_key, Tags::DEFAULT);

  // test for derivatives
  // -- wrt pressure
  std::cout << "\n========== Test 5 (Dfield/Dpressure) ==========\n ";
  eval_p->SetChanged();
  S->GetEvaluator(mass_density_key).Update(*S, "test");
  auto rho1_c = *S->Get<CV_t>(mass_density_key, tag).ViewComponent("cell");

  S->GetEvaluator(mass_density_key).UpdateDerivative(*S, "test", pressure_key, tag);
  auto drho_c = *S->GetDerivative<CV_t>(mass_density_key, tag, pressure_key, tag).ViewComponent("cell");

  S->GetEvaluator(temperature_key).Update(*S, "test");
  auto T1_c = *S->Get<CompositeVector>(temperature_key, tag).ViewComponent("cell");

  S->GetEvaluator(temperature_key).UpdateDerivative(*S, "test", pressure_key, Tags::DEFAULT);
  auto dT_c = *S->GetDerivative<CV_t>(temperature_key, tag, pressure_key, tag).ViewComponent("cell");

  double eps(1e-8);
  for (int c = 0; c < ncells; ++c) p_c[0][c] *= (1.0 + eps);

  eval_p->SetChanged();
  S->GetEvaluator(mass_density_key).Update(*S, "test");
  auto rho2_c = *S->Get<CV_t>(mass_density_key, tag).ViewComponent("cell");

  S->GetEvaluator(temperature_key).Update(*S, "test");
  auto T2_c = *S->Get<CV_t>(temperature_key, tag).ViewComponent("cell");

  double PC(22.064e+6), HC(3.7585e+4);
  double der_fd, dp, dh, tol(1e-3);
  for (int c = 0; c < ncells; ++c) {
    dp = std::fabs(p_c[0][c] - PC) / PC;
    dh = std::fabs(h_c[0][c] - HC) / HC;

    der_fd = (rho2_c[0][c] - rho1_c[0][c]) / (eps * p_c[0][c]);
    if (dp + dh > 0.25) CHECK_CLOSE(der_fd, drho_c[0][c], tol * std::fabs(der_fd));

    der_fd = (T2_c[0][c] - T1_c[0][c]) / (eps * p_c[0][c]);
    CHECK_CLOSE(der_fd, dT_c[0][c], 5 * tol * std::fabs(der_fd));
  }

  // -- wrt enthalpy
  std::cout << "\n========== Test 6 (Dfield/Denthalpy) ==========\n ";
  eval_p->SetChanged();
  eval_h->SetChanged();
  S->GetEvaluator(mass_density_key).Update(*S, "test");
  rho1_c = *S->Get<CV_t>(mass_density_key, tag).ViewComponent("cell");

  S->GetEvaluator(mass_density_key).UpdateDerivative(*S, "test", enthalpy_key, tag);
  drho_c = *S->GetDerivative<CV_t>(mass_density_key, tag, enthalpy_key, tag).ViewComponent("cell");

  S->GetEvaluator(temperature_key).Update(*S, "test");
  T1_c = *S->Get<CV_t>(temperature_key, tag).ViewComponent("cell");

  S->GetEvaluator(temperature_key).UpdateDerivative(*S, "test", enthalpy_key, Tags::DEFAULT);
  dT_c = *S->GetDerivative<CV_t>(temperature_key, tag, enthalpy_key, tag).ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) h_c[0][c] *= (1.0 + eps);

  eval_h->SetChanged();
  S->GetEvaluator(mass_density_key).Update(*S, "test");
  rho2_c = *S->Get<CV_t>(mass_density_key, tag).ViewComponent("cell");

  S->GetEvaluator(temperature_key).Update(*S, "test");
  T2_c = *S->Get<CV_t>(temperature_key, tag).ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) {
    der_fd = (rho2_c[0][c] - rho1_c[0][c]) / (eps * h_c[0][c]);
    CHECK_CLOSE(der_fd, drho_c[0][c], tol * std::fabs(der_fd));

    der_fd = (T2_c[0][c] - T1_c[0][c]) / (eps * h_c[0][c]);
    CHECK_CLOSE(der_fd, dT_c[0][c], 1e-3 * std::fabs(der_fd));
  }

  // -- wrt pressure with chain rule
  std::cout << "\n========== Test 7 (Dfield/Dpressure) ==========\n ";
  eval_p->SetChanged();

  S->GetEvaluator(conductivity_key).Update(*S, "test");
  auto tc1_c = *S->Get<CV_t>(conductivity_key, tag).ViewComponent("cell");

  S->GetEvaluator(conductivity_key).UpdateDerivative(*S, "test", pressure_key, tag);
  auto dtc_c = *S->GetDerivative<CV_t>(conductivity_key, tag, pressure_key, tag).ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) p_c[0][c] *= (1.0 + eps);

  eval_p->SetChanged();
  S->GetEvaluator(conductivity_key).Update(*S, "test");
  auto tc2_c = *S->Get<CV_t>(conductivity_key, tag).ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) {
    if (167 < c && c < 203) continue;
    der_fd = (tc2_c[0][c] - tc1_c[0][c]) / (eps * p_c[0][c]);
    CHECK_CLOSE(der_fd, dtc_c[0][c], tol * std::fabs(der_fd));
  }

  // -- wrt enthalpy with chain rule
  std::cout << "\n========== Test 7 (Dfield/Denthalpy) ==========\n ";
  eval_p->SetChanged();

  S->GetEvaluator(conductivity_key).Update(*S, "test");
  tc1_c = *S->Get<CV_t>(conductivity_key, tag).ViewComponent("cell");

  S->GetEvaluator(conductivity_key).UpdateDerivative(*S, "test", enthalpy_key, tag);
  dtc_c = *S->GetDerivative<CV_t>(conductivity_key, tag, enthalpy_key, tag).ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) h_c[0][c] *= (1.0 + eps);

  eval_h->SetChanged();
  S->GetEvaluator(conductivity_key).Update(*S, "test");
  tc2_c = *S->Get<CV_t>(conductivity_key, tag).ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) {
    der_fd = (tc2_c[0][c] - tc1_c[0][c]) / (eps * h_c[0][c]);
    CHECK_CLOSE(der_fd, dtc_c[0][c], tol * std::fabs(der_fd));
  }
}

