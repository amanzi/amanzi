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
#include "CommonDefs.hh"
#include "CompositeVector.hh"
#include "EnergyPressureEnthalpy_PK.hh"
#include "evaluators_reg.hh"
#include "IAPWS97.hh"
#include "MeshFactory.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "ThermodynamicStateEvaluators.hh"
#include "VerboseObject.hh"

TEST(EVALUATOR_DERIVATIVE_TABLES_PH)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Energy;

  using CV_t = CompositeVector;
  using CVS_t = CompositeVectorSpace;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: derivative tables" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/energy_pressure_enthalpy.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  int n = 100;
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(1.0, 0.0, 2.0, 0.2, 2 * n, 2 * n);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("energy");
  auto soln = Teuchos::rcp(new TreeVector());
  auto EPK = Teuchos::rcp(new EnergyPressureEnthalpy_PK(pk_tree, plist, S, soln));

  EPK->Setup();

  // add viscosity to state
  std::string passwd("");
  Key pressure_key = Keys::getKey("", "pressure");
  Key enthalpy_key = Keys::getKey("", "enthalpy");
  Key density_key = Keys::getKey("", "mass_density_liquid");
  Key temperature_key = Keys::getKey("", "temperature");
  Key viscosity_key = Keys::getKey("", "viscosity_liquid");
  Key conductivity_key = Keys::getKey("", "thermal_conductivity");

  S->Require<CV_t, CVS_t>(viscosity_key, Tags::DEFAULT, viscosity_key)
    .SetMesh(mesh)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

  Teuchos::ParameterList elist(viscosity_key);
  elist.set<std::string>("tag", "");
  auto eval = Teuchos::rcp(new Evaluators::ViscosityEvaluator(elist));
  S->SetEvaluator(viscosity_key, Tags::DEFAULT, eval);

  S->RequireDerivative<CV_t, CVS_t>(viscosity_key, Tags::DEFAULT, pressure_key, Tags::DEFAULT,
                                    viscosity_key).SetGhosted();
  S->RequireDerivative<CV_t, CVS_t>(viscosity_key, Tags::DEFAULT, enthalpy_key, Tags::DEFAULT,
                                    viscosity_key).SetGhosted();

  // add derivative of thermal conductivity to state
  S->RequireDerivative<CV_t, CVS_t>(conductivity_key, Tags::DEFAULT, pressure_key, Tags::DEFAULT,
                                    conductivity_key).SetGhosted();
  S->RequireDerivative<CV_t, CVS_t>(conductivity_key, Tags::DEFAULT, enthalpy_key, Tags::DEFAULT,
                                    conductivity_key).SetGhosted();

  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  EPK->Initialize();
  S->CheckAllFieldsInitialized();

  // populate (p,h) table
  auto& p_c = *S->GetW<CompositeVector>(pressure_key, pressure_key).ViewComponent("cell");
  auto& h_c = *S->GetW<CompositeVector>(enthalpy_key, passwd).ViewComponent("cell");

  AmanziEOS::IAPWS97 eos(*plist);

  int c(0);
  double scale(50.0 / n);
  double dp(5.0e-2 * scale), dh(15.0 * scale), p, h; 

  for (int i = -n; i < n; ++i) {
    for (int j = -n; j < n; ++j) {
      p = eos.PC + i * dp; // Mpa
      h = eos.HC + j * dh;

      p_c[0][c] = p * 1.0e+6;
      h_c[0][c] = h * CommonDefs::ENTHALPY_FACTOR;
      c++;
    }
  }

  Tag tag = Tags::DEFAULT;
  auto eval_p = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(S->GetEvaluatorPtr(pressure_key, tag));
  auto eval_h = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(S->GetEvaluatorPtr(enthalpy_key, tag));
  eval_p->SetChanged();
  eval_h->SetChanged();

  // compare Maxwell relations
  S->GetEvaluator("thermodynamic_state").Update(*S, "test");
  auto& state_c = *S->Get<CV_t>("thermodynamic_state", tag).ViewComponent("cell");

  c = 0;
  for (int i = -n; i < n; ++i) {
    for (int j = -n; j < n; ++j) {
      int rgn = state_c[(int)Evaluators::TS_t::RGN][c];
      double cp = state_c[(int)Evaluators::TS_t::CP][c];
      if (rgn != 4) CHECK(cp > 0.0);
      double value = (rgn != 4) ? std::log(cp) : 0.0;

      double kt = state_c[(int)Evaluators::TS_t::KT][c];
      if (rgn != 4) CHECK(kt > 0.0);
      value = (rgn != 4) ? std::log(kt) : 0.0;

      double av = state_c[(int)Evaluators::TS_t::AV][c];
      if (rgn != 4) CHECK(av > 0.0);
      value = (rgn != 4) ? std::log(av) : 0.0;

      double T = state_c[(int)Evaluators::TS_t::T][c];
      double rho = state_c[(int)Evaluators::TS_t::RHO][c];
      value = cp * kt - T * av * av / rho;
      if (rgn != 4) CHECK(value > 0.0);
      value = (rgn != 4) ? std::log(value) : 0.0;

      // std::cout << p_c[0][c] * 1e-6 << " " << h_c[0][c] / ENTHALPY_FACTOR << " " << value << std::endl;
      c++;
    }
  }

  // compute a selective derivative
  Key field = density_key;
  Key wrt = enthalpy_key;
  S->GetEvaluator(field).UpdateDerivative(*S, "test", wrt, Tags::DEFAULT);
  auto& der_c = *S->GetDerivative<CV_t>(field, tag, wrt, tag).ViewComponent("cell");
  auto& field_c = *S->Get<CV_t>(field, tag).ViewComponent("cell");

  c = 0;
  for (int i = -n; i < n; ++i) {
    for (int j = -n; j < n; ++j) {
      // std::cout << p_c[0][c] * 1e-6 << " " << h_c[0][c] / ENTHALPY_FACTOR << " " << field_c[0][c] << std::endl;
      // std::cout << p_c[0][c] * 1e-6 << " " << h_c[0][c] / ENTHALPY_FACTOR << " " << der_c[0][c] << std::endl;
      c++;
    }
  }
}

