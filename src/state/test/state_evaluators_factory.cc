/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Tests the functionality of a full example, working from input file.

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "evaluator/Evaluator.hh"


using namespace Amanzi;


TEST(EVALUATOR_FACTORY) {    // create a mesh
  auto comm = Amanzi::getDefaultComm();

  std::string xmlFileName = "test/state_evaluators_factory.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  auto gm = Teuchos::rcp(
      new AmanziGeometry::GeometricModel(2, plist.sublist("regions"), *comm));
  AmanziMesh::MeshFactory meshfac(comm, gm);
  auto mesh = meshfac.create(plist.sublist("mesh"));

  // create a state
  State S(plist.sublist("state"));
  S.RegisterDomainMesh(mesh);

  // require the time
  S.Require<double>("time", "", "time");

  // require water content and evaluator
  S.Require<CompositeVector, CompositeVectorSpace>("water_content", "")
      .SetMesh(S.GetMesh())
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S.RequireDerivative<CompositeVector,CompositeVectorSpace>("water_content", "", "pressure", "");
  S.RequireEvaluator("water_content", "");

  // setup and initialize
  S.Setup();
  S.set_time(0.);
  S.Initialize();

  // evaluate
  S.GetEvaluator("water_content").Update(S, "test_request");
  {
    auto cvv =
        S.Get<CompositeVector>("water_content", "").ViewComponent<HostSpaceSpecial>(
            "cell", 0, false);
    CHECK_CLOSE(1212.946, cvv(0), 1.e-3);
    CHECK_CLOSE(9252.804, cvv(1), 1.e-3);
  }

  CHECK(S.GetEvaluator("water_content", "").IsDifferentiableWRT(S, "pressure", ""));
  S.GetEvaluator("water_content").UpdateDerivative(S, "test_request", "pressure", "");
  {
    auto cvv = S.GetDerivative<CompositeVector>("water_content", "", "pressure", "")
               .ViewComponent<HostSpaceSpecial>("cell", 0, false);
    CHECK_CLOSE(0.10689, cvv(0), 1.e-5);
    CHECK_CLOSE(0.292203, cvv(1), 1.e-5);
  }
}

  
