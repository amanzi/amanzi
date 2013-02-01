/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

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

#include "BDF2_TI.hh"
#include "MeshFactory.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Flow_State.hh"
#include "Richards_PK.hh"


/* **************************************************************** */
TEST(FLOW_2D_RICHARDS_SEEPAGE) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) cout << "Test: 2D Richards, seepage boundary condition" << endl;

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_seepage.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<Mesh> mesh = meshfactory(0.0, 0.0, 100.0, 50.0, 100, 50, gm); 

  // create and populate flow state
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(mesh));
  FS->set_permeability(5.0e-13, 5.0e-14, "Material 1");
  FS->set_porosity(0.2);
  FS->set_fluid_viscosity(0.00089);
  FS->set_fluid_density(998);
  FS->set_gravity(-9.81);

  // create Richards process kernel
  Richards_PK* RPK = new Richards_PK(parameter_list, FS);
  RPK->InitPK();
  RPK->InitSteadyState(0.0, 0.01);

  // solve the steady-state problem
  RPK->AdvanceToSteadyState(0.0, 0.01);
  RPK->CommitState(FS);

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(FS->ref_pressure(), 0, "pressure");
    GMV::write_cell_data(FS->ref_water_saturation(), 0, "saturation");
    GMV::close_data_file();
  }

  delete RPK;
}
