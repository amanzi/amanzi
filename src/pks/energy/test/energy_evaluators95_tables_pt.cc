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
#include <chrono>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CommonDefs.hh"
#include "CompositeVector.hh"
#include "EnergyPressureTemperature_PK.hh"
#include "evaluators_reg.hh"
#include "IAPWS95.hh"
#include "IAPWS95_RaggedSplineRhoT.hh"
#include "MeshFactory.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "IAPWS95_StateEvaluators.hh"
#include "VerboseObject.hh"

double RunTest(int icase, const std::string& iapws95_model)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Energy;
  using namespace Amanzi::Evaluators;

  using CV_t = CompositeVector;
  using CVS_t = CompositeVectorSpace;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: derivative tables: spline=" << (icase > 0) << std::endl;

  // read parameter list
  std::string xmlFileName = "test/energy_iapws95.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  plist->sublist("state").sublist("initial conditions").sublist("pressure")
    .sublist("function").sublist("domain").sublist("function")
    .sublist("function-constant").set<double>("value", 39.0e+6);

  AmanziEOS::IAPWS95 eos(*plist);

  plist->sublist("PKs").sublist("energy").sublist("thermal conductivity evaluator")
    .sublist("All").sublist("liquid phase")
    .set<bool>(iapws95_model, true);
  plist->sublist("state").sublist("evaluators").sublist("thermodynamic_state")
    .set<bool>(iapws95_model, true);
  plist->sublist("state").sublist("evaluators").sublist("viscosity_liquid")
    .set<bool>(iapws95_model, true);

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  int n = 100;
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(1.0, 0.0, 2.0, 0.2, n, n);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("energy");
  auto soln = Teuchos::rcp(new TreeVector());
  auto EPK = Teuchos::rcp(new EnergyPressureTemperature_PK(pk_tree, plist, S, soln));

  EPK->Setup();

  std::string passwd("");
  Key pressure_key = Keys::getKey("", "pressure");
  Key temperature_key = Keys::getKey("", "temperature");
  Key enthalpy_key = Keys::getKey("", "enthalpy");
  Key state_key = Keys::getKey("", "thermodynamic_state");
  Key density_key = Keys::getKey("", "mass_density_liquid");
  Key viscosity_key = Keys::getKey("", "viscosity_liquid");
  Key conductivity_key = Keys::getKey("", "thermal_conductivity");

  // add viscosity to state
  S->Require<CV_t, CVS_t>(viscosity_key, Tags::DEFAULT, viscosity_key)
    .SetMesh(mesh)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
  S->RequireEvaluator(viscosity_key, Tags::DEFAULT);

  S->RequireDerivative<CV_t, CVS_t>(viscosity_key, Tags::DEFAULT, pressure_key, Tags::DEFAULT,
                                    viscosity_key).SetGhosted();
  S->RequireDerivative<CV_t, CVS_t>(conductivity_key, Tags::DEFAULT, pressure_key, Tags::DEFAULT,
                                    conductivity_key).SetGhosted();

  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  EPK->Initialize();
  S->CheckAllFieldsInitialized();

  // populate (p,T) table
  double p_min = 0.2;
  double p_max = 50.0;
  double T_min = 280.0;
  double T_max = 950.0;

  auto& p_c = *S->GetW<CompositeVector>(pressure_key, pressure_key).ViewComponent("cell");
  auto& T_c = *S->GetW<CompositeVector>(temperature_key, passwd).ViewComponent("cell");

  int c = 0;
  for (double i = 0; i < n; i++) {
    for (double j = 0; j < n; j++) {
      p_c[0][c] = (p_min + (p_max - p_min) * i / double(n)) * 1e+6;
      T_c[0][c] = (T_min + (T_max - T_min) * j / double(n));
      c++;
    }
  }

  auto start = std::chrono::steady_clock::now();

  Tag tag = Tags::DEFAULT;
  auto eval_p = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(S->GetEvaluatorPtr(pressure_key, tag));
  auto eval_T = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(S->GetEvaluatorPtr(temperature_key, tag));
  eval_p->SetChanged();
  eval_T->SetChanged();

  // compute a selective derivative
  S->GetEvaluator(density_key).UpdateDerivative(*S, "test", pressure_key, Tags::DEFAULT);
  auto& drhodp = *S->GetDerivative<CV_t>(density_key, tag, pressure_key, tag).ViewComponent("cell");

  S->GetEvaluator(density_key).UpdateDerivative(*S, "test", temperature_key, Tags::DEFAULT);
  auto& drhodT = *S->GetDerivative<CV_t>(density_key, tag, temperature_key, tag).ViewComponent("cell");

  S->GetEvaluator(enthalpy_key).UpdateDerivative(*S, "test", temperature_key, Tags::DEFAULT);
  auto& dhdT = *S->GetDerivative<CV_t>(enthalpy_key, tag, temperature_key, tag).ViewComponent("cell");

  auto& state_c = *S->Get<CV_t>(state_key, tag).ViewComponent("cell");

  auto end = std::chrono::steady_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Elapsed time: " << elapsed.count() << " ms\n";

  c = 0;
  std::ofstream out("field.dat");
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      CHECK(drhodp[0][c] > 0.0);
      CHECK(drhodT[0][c] < 0.0);
      CHECK(dhdT[0][c] > 0.0);
      out << T_c[0][c] << " " << p_c[0][c] * 1e-6 << " " << drhodT[0][c] << std::endl;
      // out << T_c[0][c] << " " << p_c[0][c] * 1e-6 << " " << state_c[(int)TSPH_t::RHO][c] << std::endl;
      c++;
    }
  }
  out.close();

  return elapsed.count();
}

TEST(EVALUATOR_DERIVATIVE95_TABLES_PT)
{
  double t0 = RunTest(0, "use iapws95");
  double t2 = RunTest(2, "use iapws95 spline rho/T");
  CHECK(t0 > 2 * t2);
}
