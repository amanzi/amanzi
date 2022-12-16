/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
   Flow PK

*/

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "AmanziComm.hh"
#include "Epetra_SerialComm.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

// Flow
#include "Richards_PK.hh"

#include "Analytic01.hh"
#include "Richards_SteadyState.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Flow;

void
RunTestConvergence(std::string input_xml)
{
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0)
    std::cout << "\nConvergence analysis on random meshes: " << input_xml << std::endl;

  std::string xmlFileName = input_xml;
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // convergence estimate
  int nmeshes = plist->get<int>("number of meshes", 1);
  std::vector<double> h, p_error, v_error;

  for (int n = 0; n < nmeshes; n++) { // Use "n < 3" for the full test
    Teuchos::ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

    Preference pref;
    pref.clear();
    pref.push_back(Framework::MSTK);
    pref.push_back(Framework::STK);

    MeshFactory meshfactory(comm, gm);
    meshfactory.set_preference(pref);
    Teuchos::RCP<const Mesh> mesh;
    if (n == 0) {
      //mesh = meshfactory.create("test/test_nice.exo");
      mesh = meshfactory.create("test/random_mesh1.exo");
    } else if (n == 1) {
      mesh = meshfactory.create("test/random_mesh2.exo");
    } else if (n == 2) {
      mesh = meshfactory.create("test/random_mesh3.exo");
    }

    /* create a simple state and populate it */
    Amanzi::VerboseObject::global_hide_line_prefix = false;

    Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
    Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
    Richards_PK* RPK = new Richards_PK(plist, "flow", S, soln);

    RPK->Setup();
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();

    // create Richards process kernel
    RPK->Initialize();
    S->CheckAllFieldsInitialized();

    // solver the problem
    TI_Specs ti_specs;
    ti_specs.T0 = 0.0;
    ti_specs.dT0 = 1.0;
    ti_specs.T1 = 1.0e+5;
    ti_specs.max_itrs = 2000;

    AdvanceToSteadyState(S, *RPK, ti_specs, soln);
    RPK->CommitStep(0.0, 1.0, Tags::DEFAULT);

    S->Get<CompositeVector>("volumetric_flow_rate").ScatterMasterToGhosted("face");
    const auto& p = *S->Get<CompositeVector>("pressure").ViewComponent("cell");
    const auto& flux = *S->Get<CompositeVector>("volumetric_flow_rate").ViewComponent("face", true);

    double pressure_err, flux_err, div_err; // error checks
    pressure_err = CalculatePressureCellError(mesh, p);
    flux_err = CalculateDarcyFluxError(mesh, flux);
    div_err = CalculateDarcyDivergenceError(mesh, flux);

    int num_bdf1_steps = ti_specs.num_itrs;
    printf("mesh=%d bdf1_steps=%d  L2_pressure_err=%7.3e  l2_flux_err=%7.3e  L2_div_err=%7.3e\n",
           n,
           num_bdf1_steps,
           pressure_err,
           flux_err,
           div_err);

    CHECK(pressure_err < 2.2e-1 && flux_err < 2e-1 && div_err < 2e-1);

    GMV::open_data_file(*mesh, (std::string) "flow_richards.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();

    delete RPK;
  }
}


/* *****************************************************************
* Run with various discretization methods
* **************************************************************** */
TEST(FLOW_RICHARDS_CONVERGENCE_NLFV)
{
  RunTestConvergence("test/flow_richards_random_nlfv.xml");
}

TEST(FLOW_RICHARDS_CONVERGENCE_FV)
{
  RunTestConvergence("test/flow_richards_random_fv.xml");
}

TEST(FLOW_RICHARDS_CONVERGENCE_MFD)
{
  RunTestConvergence("test/flow_richards_random.xml");
}
