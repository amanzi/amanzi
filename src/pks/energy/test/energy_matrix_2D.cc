/*
  This is the energy component of the Amanzi code. 

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
#include "EnergyTwoPhase_PK.hh"

/* **************************************************************** */
TEST(ENERGY_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Energy;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: 2D Thermal, homogeneous medium" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/energy_matrix_2D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::RCP<const Teuchos::ParameterList> plist = 
      Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  /* create a mesh framework */
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 18, 18, gm);

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  Teuchos::ParameterList state_list;
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Energy_PK* EPK = new EnergyTwoPhase_PK(plist, S);
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  EPK->InitializeFields();
  S->CheckAllFieldsInitialized();

  /* modify the default state for the problem at hand */
  /* create the initial temperature function */
  std::string passwd("state");
  Epetra_MultiVector& temperature = *S->GetFieldData("temperature", passwd)->ViewComponent("cell");

  /* initialize the Energy process kernel */

  /* solve the problem */

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"energy.gmv");
    GMV::start_data();
    GMV::write_cell_data(temperature, 0, "temperature");
    GMV::close_data_file();
  }

  delete EPK;
}
