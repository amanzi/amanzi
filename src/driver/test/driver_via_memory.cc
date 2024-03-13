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

  auto plist = simulator.get_plist();
  plist->sublist("state")
    .sublist("evaluators")
    .sublist("fracture-aperture")
    .set<std::string>("evaluator type", "aperture")
    .set<std::string>("pressure key", "fracture-pressure");

  Teuchos::ParameterList state_plist = plist->sublist("state");
  auto S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);
  if (submesh.get()) S->RegisterMesh(domain, submesh);

  Amanzi::ObservationData obs_data;
  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  // test aperture opening (4% error)
  VerifyDriverAperture(*S, 0.04);
}


TEST(DRIVER_VIA_MEMORY_EVALUATOR)
{
  RunTest("test/driver_via_memory.xml");
}
