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

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Flow_State.hh"
#include "Richards_PK.hh"


/* **************************************************************** */
TEST(FLOW_3D_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) cout << "Test: 3D Richards, crib model" << endl;

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_bc_cribs.xml";
  
  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  // create a mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);

  MeshFactory factory(&comm);
  ParameterList mesh_list = parameter_list.get<ParameterList>("Mesh").get<ParameterList>("Unstructured");
  ParameterList factory_list = mesh_list.get<ParameterList>("Generate Mesh");
  Teuchos::RCP<Mesh> mesh(factory(factory_list, gm));

  // create flow state
  ParameterList state_list = parameter_list.get<ParameterList>("State");
  State S(state_list, mesh);
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(S));

  // create Richards process kernel
  Richards_PK* RPK = new Richards_PK(parameter_list, FS);
  RPK->InitPK();
  RPK->InitSteadyState(0.0, 1e-7);  // dT0 is not used
  RPK->PrintStatistics();

  // solve the problem
  S.set_time(0.0);
  RPK->AdvanceToSteadyState(0.0, 1e-7);
  RPK->CommitState(FS);

  // derive dependent variable
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector  saturation(pressure);
  RPK->DeriveSaturationFromPressure(pressure, saturation);

  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(pressure, 0, "pressure");
  GMV::write_cell_data(saturation, 0, "saturation");
  GMV::write_cell_data(*(S.get_vertical_permeability()), 0, "vert_permeability");
  GMV::close_data_file();

  // check the pressure profile
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  //for (int c = 0; c < ncells; c++) cout << (mesh->cell_centroid(c))[2] << " " << pressure[c] << endl;
  for (int c = 0; c < ncells; c++) CHECK(pressure[c] > 4500.0 && pressure[c] < 101325.0);

  delete RPK;
}
