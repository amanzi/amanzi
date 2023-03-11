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
#include "CycleDriver.hh"
#include "MeshAudit.hh"
#include "energy_tcm_registration.hh"
#include "energy_iem_registration.hh"
#include "eos_registration.hh"
#include "evaluators_flow_reg.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "models_flow_reg.hh"
#include "mpc_pks_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_energy_registration.hh"
#include "pks_flow_reg.hh"
#include "pks_transport_registration.hh"
#include "State.hh"

// MPC
#include "EvaluatorAperture.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

Utils::RegisteredFactory<Evaluator, EvaluatorAperture> EvaluatorAperture::reg_("aperture");

void
RunTest(const std::string xmlInFileName)
{
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  // setup a piecewice linear solution with a jump
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(-12.5, -0.5, -0.5, 12.5, 0.5, 0.5, 125, 1, 4, true, true);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  // auto mesh_fracture = factory.create(mesh, names, AmanziMesh::FACE);
  auto mesh_fracture = Teuchos::rcp(new MeshExtractedManifold(
    mesh, "fracture", AmanziMesh::FACE, comm, gm, mesh_list, true, false));

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();


  // test aperture opening (5% error)
  double t = 1000.0;
  const auto& a_c = *S->Get<CompositeVector>("fracture-aperture").ViewComponent("cell");

  for (int c = 0; c < a_c.MyLength(); ++c) {
    const auto& xc = mesh_fracture->cell_centroid(c);
    double w = (xc[0] + 12.5) * std::sqrt(2000.0 / (t + 1e-10)) - 12.5;
    double tmp = (w > 7.4) ? 1.0e-5 :
                             1.257136447580295e-05 - 5.307442110525289e-07 * w +
                               1.809065105747629e-08 * w * w + 1.184392778715966e-09 * w * w * w -
                               3.620526561612838e-11 * w * w * w * w;
    if (c < 10 && MyPID == 0)
      std::cout << "aperture: " << a_c[0][c] << ",  exact: " << tmp << std::endl;
    CHECK(std::fabs(a_c[0][c] - tmp) < 2e-2 * tmp);
  }
}


TEST(MPC_DRIVER_COUPLED_FLOW_APERTURE_DARCY)
{
  RunTest("test/mpc_coupled_flow_aperture_darcy.xml");
}

TEST(MPC_DRIVER_COUPLED_FLOW_APERTURE_RICHARDS)
{
  RunTest("test/mpc_coupled_flow_aperture.xml");
}
