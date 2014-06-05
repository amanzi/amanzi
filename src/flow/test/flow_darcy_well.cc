/*
  This is the flow component of the Amanzi code. 

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
#include "MeshAudit.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Darcy_PK.hh"


/* **************************************************************** */
TEST(FLOW_2D_DARCY_WELL) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: 2D specific storage Darcy, homogeneous medium" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/flow_darcy_well.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* create an MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(-10.0, -5.0, 10.0, 0.0, 200, 50, gm);

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Darcy_PK* DPK = new Darcy_PK(plist, S);
  S->Setup();
  S->InitializeFields();

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell", false);
  
  for (int c = 0; c < K.MyLength(); c++) {
    K[0][c] = 0.1;
    K[1][c] = 2.0;
  }

  *S->GetScalarData("fluid_density", passwd) = 1.0;
  *S->GetScalarData("fluid_viscosity", passwd) = 1.0;
  Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", passwd);
  gravity[1] = -1.0;
  S->GetFieldData("porosity", passwd)->PutScalar(0.2);
  S->GetFieldData("specific_storage", passwd)->PutScalar(0.1);

  /* initialize the Darcy process kernel */
  DPK->InitPK();
  DPK->InitTransient(0.0, 1e-8);

  /* transient solution */
  double dT = 0.5;
  for (int n = 0; n < 10; n++) {
    double dT_actual(dT);
    DPK->Advance(dT, dT_actual);
    DPK->CommitState(S);

    if (MyPID == 0) {
      const Epetra_MultiVector& p = *S->GetFieldData("pressure")->ViewComponent("cell");
      GMV::open_data_file(*mesh, (std::string)"flow.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "pressure");
      GMV::close_data_file();
    }
  }

  delete DPK;
}
