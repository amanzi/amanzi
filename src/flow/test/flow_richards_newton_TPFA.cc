/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Daniil Svyatskiy (dasvyat@lanl.gov)
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

#include "BDF1_TI.hh"

/* ******************************** */
TEST(NEWTON_RICHARD_STEADY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) cout << "Test: orthogonal newton solver, 2-layer model" << endl;
   /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_newton_TPFA.xml";

  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  // create a mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(&comm);
  factory.preference(pref);
  ParameterList mesh_list = parameter_list.get<ParameterList>("Mesh").get<ParameterList>("Unstructured");
  ParameterList factory_list = mesh_list.get<ParameterList>("Generate Mesh");
  Teuchos::RCP<Mesh> mesh(factory(factory_list, gm));

  Teuchos::ParameterList state_list = parameter_list.get<Teuchos::ParameterList>("State");
  State S(state_list);
  S.RegisterDomainMesh(mesh);
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(S));
  S.Setup();
  S.Initialize();
  FS->Initialize();

  // create Richards process kernel
  Richards_PK* RPK = new Richards_PK(parameter_list, FS);
  RPK->InitPK();


  
  RPK->InitSteadyState(0.0, 1e+4);

  // solve the problem
  // S.set_time(0.0);
  RPK->AdvanceToSteadyState(0.0, 1000.0);
  RPK->CommitState(FS);

  // derive dependent variable
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector  saturation(pressure);
  RPK->DeriveSaturationFromPressure(pressure, saturation);

  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(pressure, 0, "pressure");
  GMV::write_cell_data(saturation, 0, "saturation");
  GMV::close_data_file();

  // check the number of iteration
  // int numiter = RPK->bdf1_dae->total_non_iter;

  // BDF1Dae* bdf1 = RPK->time_intergrator();
  // cout << "nonlinear numiter " << bdf1->total_nonlinear_iter() << endl;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    // cout << (mesh->cell_centroid(c))[2] << " " << pressure[c] << endl;
    CHECK(pressure[c] > 4500.0 && pressure[c] < 101325.0);
  }
  
  delete RPK;
}

