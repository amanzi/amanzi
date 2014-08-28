/*
  This is the flow test of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Richards_PK.hh"


/* **************************************************************** */
TEST(FLOW_2D_RICHARDS_SEEPAGE) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: 2D Richards, seepage boundary condition" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/flow_richards_seepage.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 100.0, 50.0, 100, 50, gm); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  ParameterList state_list;
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Richards_PK* RPK = new Richards_PK(plist, S);
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  RPK->InitializeFields();
  S->CheckAllFieldsInitialized();

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell");
  
  for (int c = 0; c != K.MyLength(); ++c) {
    K[0][c] = 5.0e-13;
    K[1][c] = 5.0e-14;
  }

  double rho = *S->GetScalarData("fluid_density", passwd) = 998.0;
  *S->GetScalarData("fluid_viscosity", passwd) = 0.00089;
  Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", passwd);
  double g = gravity[1] = -9.81;

  /* create the initial pressure function */
  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell");
  Epetra_MultiVector& lambda = *S->GetFieldData("pressure", passwd)->ViewComponent("face");

  double p0(101325.0), z0(30.0);
  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[0][c] = p0 + rho * g * (xc[1] - z0);
  }
  RPK->DeriveFaceValuesFromCellValues(p, lambda); 

  /* create Richards process kernel */
  RPK->Initialize(S.ptr());
  RPK->ti_specs_sss().T1 = 1e+10;
  RPK->ti_specs_sss().dTmax = 1e+8;
  RPK->ti_specs_sss().residual_tol = 1e-9;
  RPK->ti_specs_sss().max_itrs = 30;

  RPK->InitSteadyState(0.0, 0.01);

  /* solve the steady-state problem */
  RPK->AdvanceToSteadyState(0.0, 0.01);
  RPK->CommitState(0.0, S.ptr());

  const Epetra_MultiVector& ws = *S->GetFieldData("water_saturation")->ViewComponent("cell");
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::write_cell_data(ws, 0, "saturation");
    GMV::close_data_file();
  }

  delete RPK;
}
