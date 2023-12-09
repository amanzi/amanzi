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
#include "xercesc/util/PlatformUtils.hpp"
#include "xercesc/parsers/XercesDOMParser.hpp"

// Amanzi
#include "AmanziUnstructuredGridSimulationDriver.hh"
#include "IO.hh"
#include "InputConverterU.hh"
#include "CycleDriver.hh"
#include "MeshAudit.hh"
#include "eos_reg.hh"
#include "evaluators_flow_reg.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "models_flow_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"

#include "test/EvaluatorAperture.hh"
#include "VerifyDriver.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

Utils::RegisteredFactory<Evaluator, EvaluatorAperture> EvaluatorAperture::reg_("aperture");

void
RunTest(const std::string& xmlInFileName)
{
  xercesc_3_2::XercesDOMParser* parser = Amanzi::AmanziInput::CreateXMLParser();
  xercesc_3_2::DOMDocument* doc = Amanzi::AmanziInput::OpenXMLInput(parser, xmlInFileName);
  AmanziUnstructuredGridSimulationDriver simulator(xmlInFileName, doc, "");

  auto comm = Amanzi::getDefaultComm();
  simulator.set_comm(comm);

  std::string domain;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> submesh;

  auto gm = simulator.InitGeometricModel();
  simulator.InitMesh(gm, mesh, domain, submesh);

  Key key("fracture-aperture");
  auto plist = simulator.get_plist();

  plist->sublist("state")
    .sublist("evaluators")
    .sublist(key)
    .set<std::string>("evaluator type", "aperture")
    .set<std::string>("pressure key", "fracture-pressure");

  Teuchos::ParameterList state_plist = plist->sublist("state");
  auto S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);
  if (submesh.get()) S->RegisterMesh(domain, submesh);

  Amanzi::ObservationData obs_data;
  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);

  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(submesh)->SetGhosted(true);
  cvs->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  double dt0(-1.0), dtg(1.0), t0, t1(0.0);
  while (t1 < 1000) {
    t0 = t1;
    t1 = t0 + dtg;
    std::cout << "========================================\n"
              << "New loop from " << t0 << " to " << t1
              << "\n========================================\n";

    cycle_driver.Go(t0, t1, &dt0);
    dtg *= 1.5;

    // resetting aperture to mimic coupling
    Teuchos::ParameterList elist(key);
    elist.set<std::string>("my key", key)
      .set<std::string>("evaluator type", "aperture")
      .set<std::string>("pressure key", "fracture-pressure")
      .set<double>("normal stiffness", 1.0e+11)
      .set<std::string>("tag", "");

    auto eval = Teuchos::rcp(new EvaluatorAperture(elist));
    S->SetEvaluator(key, Tags::DEFAULT, eval);
  }

  // test aperture opening (5% error)
  VerifyDriverAperture(*S, 0.05);
}


TEST(DRIVER_VIA_MEMORY_LOOP)
{
  RunTest("test/driver_via_memory.xml");
}
