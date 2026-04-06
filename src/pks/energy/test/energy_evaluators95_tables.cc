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
#include "EnergyPressureTemperature_PK.hh"
#include "evaluators_reg.hh"
#include "IAPWS95.hh"
#include "MeshFactory.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "IAPWS95_StateEvaluators.hh"
#include "VerboseObject.hh"

TEST(EVALUATOR_DERIVATIVE95_TABLES_PT)
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
  std::string xmlFileName = "test/energy_iapws95.xml";
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
  auto EPK = Teuchos::rcp(new EnergyPressureTemperature_PK(pk_tree, plist, S, soln));

  EPK->Setup();

  std::string passwd("");
  Key pressure_key = Keys::getKey("", "pressure");
  Key temperature_key = Keys::getKey("", "temperature");
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
  auto& p_c = *S->GetW<CompositeVector>(pressure_key, pressure_key).ViewComponent("cell");
  auto& T_c = *S->GetW<CompositeVector>(temperature_key, passwd).ViewComponent("cell");

  AmanziEOS::IAPWS95 eos(*plist);

  int c(0);
  double scale(50.0 / n);
  double dp(5.0e-2 * scale), dT(1.0 * scale), p, T; 

  for (int i = -n; i < n; ++i) {
    for (int j = -n; j < n; ++j) {
      p = eos.PC + i * dp; // Mpa
      T = eos.TC + j * dT + 0.11;

      p_c[0][c] = p * 1.0e+6;
      T_c[0][c] = T;
      c++;
    }
  }

  Tag tag = Tags::DEFAULT;
  auto eval_p = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(S->GetEvaluatorPtr(pressure_key, tag));
  auto eval_T = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(S->GetEvaluatorPtr(temperature_key, tag));
  eval_p->SetChanged();
  eval_T->SetChanged();

  // compute a selective derivative
  Key field = viscosity_key;
  Key wrt = pressure_key;
  S->GetEvaluator(field).UpdateDerivative(*S, "test", wrt, Tags::DEFAULT);
  auto& der_c = *S->GetDerivative<CV_t>(field, tag, wrt, tag).ViewComponent("cell");
  auto& field_c = *S->Get<CV_t>(field, tag).ViewComponent("cell");

  c = 0;
  for (int i = -n; i < n; ++i) {
    for (int j = -n; j < n; ++j) {
      // std::cout << p_c[0][c] * 1e-6 << " " << T_c[0][c] << " " << field_c[0][c] << std::endl;
      // std::cout << p_c[0][c] * 1e-6 << " " << T_c[0][c] << " " << der_c[0][c] << std::endl;
      c++;
    }
  }
}

