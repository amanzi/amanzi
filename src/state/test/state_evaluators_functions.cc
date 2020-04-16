/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Tests the math functionality in some basic evaluators.

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "evaluator/EvaluatorIndependentFunction.hh"
#include "evaluator/EvaluatorIndependentFromFile.hh"
#include "evaluator/EvaluatorSecondaryMonotypeAdditive.hh"
#include "evaluator/EvaluatorSecondaryMonotypeMultiplicative.hh"
#include "evaluator/EvaluatorSecondaryMonotypeFromFunction.hh"

using namespace Amanzi;

struct tester {
  tester()
  {
    std::string xmlFileName = "test/state_evaluators_functions.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    Teuchos::ParameterList plist = xmlreader.getParameters();

    S = std::unique_ptr<State>(new State(plist.sublist("state")));

    auto comm = getDefaultComm();
    auto gm = Teuchos::rcp(
      new AmanziGeometry::GeometricModel(2, plist.sublist("regions"), *comm));
    AmanziMesh::MeshFactory meshfactory(comm, gm);
    auto mesh = meshfactory.create(0., -1., 4., 1., 2, 1);
    S->RegisterDomainMesh(mesh);

    times = { 0., 1., 2., 3., 4. };
  }

  std::unique_ptr<State> S;
  std::vector<double> times;

  void require(const std::string& field)
  {
    S->Require<CompositeVector, CompositeVectorSpace>(field)
      .SetMesh(S->GetMesh())
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireEvaluator(field);
  }

  void setup()
  {
    S->Setup();
    S->Initialize();
  }

  void run_test(const std::string& fname,
                const std::vector<std::array<double, 2>>& results)
  {
    std::cout << std::endl
              << std::endl
              << "Running test: " << fname << std::endl
              << "-----------------------------------" << std::endl;
    require(fname);
    setup();

    for (int i = 0; i != 5; ++i) {
      S->set_time(times[i]);
      S->GetEvaluator(fname).Update(*S, "test");

      auto fv = S->Get<CompositeVector>(fname).ViewComponent<DefaultHost>(
        "cell", false);
      CHECK_CLOSE(results[i][0], fv(0, 0), 1.e-10);
      CHECK_CLOSE(results[i][1], fv(1, 0), 1.e-10);
    }
  }
};


SUITE(STATE_EVALUATORS_FUNCTIONS)
{
  TEST_FIXTURE(tester, TEST_THE_TEST)
  {
    CHECK_CLOSE(1.0, S->GetMesh()->cell_centroid(0)[0], 1.e-10);
    CHECK_CLOSE(3.0, S->GetMesh()->cell_centroid(1)[0], 1.e-10);
  }

  TEST_FIXTURE(tester, INDEPENDENT_FROM_FUNCTION)
  {
    std::vector<std::array<double, 2>> result = {
      { 1, 3 }, { 1, 3 }, { 1, 3 }, { 1, 3 }, { 1, 3 }
    };

    run_test("field1", result);
  }


  TEST_FIXTURE(tester, INDEPENDENT_FROM_FILE)
  {
    std::vector<std::array<double, 2>> result = {
      { 2, 4 }, { 2, 4 }, { 3, 5 }, { 4, 6 }, { 4, 6 }
    };

    run_test("field2", result);
  }


  TEST_FIXTURE(tester, INDEPENDENT_FROM_FILE_WITH_TIME_FUNC)
  {
    std::vector<std::array<double, 2>> result = {
      { 2, 4 }, { 3, 5 }, { 4, 6 }, { 4, 6 }, { 4, 6 }
    };

    run_test("field3", result);
  }


  TEST_FIXTURE(tester, SECONDARY_ADDITIVE)
  {
    std::vector<std::array<double, 2>> result = {
      { 4, 10 }, { 4, 10 }, { 5, 11 }, { 6, 12 }, { 6, 12 }
    };
    run_test("field4", result);
  }

  TEST_FIXTURE(tester, SECONDARY_MULTIPLICATIVE)
  {
    std::vector<std::array<double, 2>> result = {
      { 2, 12 }, { 2, 12 }, { 3, 15 }, { 4, 18 }, { 4, 18 }
    };
    run_test("field5", result);
  }

  TEST_FIXTURE(tester, SECONDARY_FUNCTION)
  {
    std::vector<std::array<double, 2>> result = {
      { 8, 18 }, { 8, 18 }, { 11, 21 }, { 14, 24 }, { 14, 24 }
    };
    run_test("field6", result);
  }
}
