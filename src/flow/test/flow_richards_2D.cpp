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

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "gmv_mesh.hh"

#include "State.hpp"
#include "Flow_State.hpp"
#include "Richards_PK.hpp"
#include "BDF2_Dae.hpp"


/* **************************************************************** */
TEST(FLOW_2D_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) cout << "Test: 2D Richards, 2-layer model" << endl;

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_2D.xml";
  // DEPRECATED  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  // create an MSTK mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<Mesh> mesh = meshfactory(0.0, -2.0, 1.0, 0.0, 18, 18, gm);

  // create and populate flow state
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(mesh));
  FS->set_permeability(0.1, 2.0, "Material 1");
  FS->set_permeability(0.5, 0.5, "Material 2");
  FS->set_porosity(0.2);
  FS->set_fluid_viscosity(1.0);
  FS->set_fluid_density(1.0);
  FS->set_gravity(-1.0);

  Epetra_Vector& p = FS->ref_pressure();
  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[c] = xc[1] * (xc[1] + 2.0);
  }

  // create Richards process kernel
  Richards_PK* RPK = new Richards_PK(parameter_list, FS);
  RPK->InitPK();
  RPK->InitializeAuxiliaryData();
  RPK->InitSteadyState(0.0, 1e-8);
  RPK->ResetErrorControl(AmanziFlow::FLOW_TI_ERROR_CONTROL_PRESSURE);

  // solve the problem
  RPK->AdvanceToSteadyState(0.0, 0.1);
  RPK->CommitState(FS);

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(FS->ref_pressure(), 0, "pressure");
    GMV::close_data_file();
  }


  // check the pressure 
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) CHECK(p[c] > 0.0 && p[c] < 2.0);


  delete RPK;
}
