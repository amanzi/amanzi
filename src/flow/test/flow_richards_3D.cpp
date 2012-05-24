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
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh.hh"
#include "Mesh_STK.hh"
#include "gmv_mesh.hh"

#include "State.hpp"
#include "Flow_State.hpp"
#include "Richards_PK.hpp"
#include "BDF2_Dae.hpp"


/* **************************************************************** */
TEST(FLOW_3D_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) cout << "Test: 3D parallel Richards" << endl;

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_3D.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create a mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);
  ParameterList unstructured_list = parameter_list.get<ParameterList>("Mesh").get<ParameterList>("Unstructured");
  ParameterList mesh_list = unstructured_list.get<ParameterList>("Generate Mesh");
  RCP<Mesh> mesh = rcp(new Mesh_STK(mesh_list, &comm, gm));

  // create flow state
  ParameterList state_list = parameter_list.get<ParameterList>("State");
  State* S = new State(state_list, mesh);
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(*S));

  // create Richards process kernel
  Richards_PK* RPK = new Richards_PK(parameter_list, FS);
  RPK->set_standalone_mode(true);
  RPK->InitPK();
  RPK->InitSteadyState(0.0, 1e-8);
  RPK->PrintStatistics();

  // solve the problem
  S->set_time(0.0);
  RPK->AdvanceToSteadyState();
  RPK->CommitStateForTransport(FS);

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(FS->ref_pressure(), 0, "pressure");
    GMV::write_cell_data(*(S->get_vertical_permeability()), 0, "vert_permeability");
    GMV::close_data_file();
  }

  delete gm;
  delete RPK;
  delete S;
}
