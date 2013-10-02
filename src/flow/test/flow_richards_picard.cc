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
#include "SolverFnPicard.hh"
#include "SolverNewton.hh"
#include "Flow_State.hh"
#include "Richards_PK.hh"


/* **************************************************************** */
TEST(FLOW_RICHARDS_PICARD) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) cout << "Test: 2D crib model with Picard solver" << endl;

  /* read parameter list */
  string xmlFileName = "test/flow_richards_picard.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(&comm);
  factory.preference(pref);  
  ParameterList mesh_list = plist.get<ParameterList>("Mesh").get<ParameterList>("Unstructured");
  ParameterList factory_list = mesh_list.get<ParameterList>("Generate Mesh");
  Teuchos::RCP<Mesh> mesh(factory(factory_list, gm));

  // create flow state
  ParameterList state_list = plist.get<ParameterList>("State");
  State S(state_list);
  S.RegisterDomainMesh(mesh);
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(S));
  S.Setup();
  FS->Initialize();
  S.Initialize();

  // create Richards process kernel
  Teuchos::RCP<Richards_PK> RPK = Teuchos::rcp(new Richards_PK(plist, FS));
  RPK->InitPK();
  RPK->InitSteadyState(0.0, 1e-7);  // dT0 is not used

  // create the function class
  Teuchos::RCP<SolverFnPicard<Epetra_Vector> > fn = 
      Teuchos::rcp(new SolverFnPicard<Epetra_Vector>(mesh, RPK));

  // create the Solver
  const Teuchos::RCP<Epetra_Vector> u = RPK->get_solution();
  const Epetra_BlockMap& map = u->Map();
  Teuchos::ParameterList nlist = RPK->nonlin_solver_list_;

  Teuchos::RCP<AmanziSolvers::SolverNewton<Epetra_Vector, Epetra_BlockMap> > picard =
      Teuchos::rcp(new AmanziSolvers::SolverNewton<Epetra_Vector, Epetra_BlockMap>(nlist, fn, map));

  picard->Solve(u);

  RPK->CommitState(FS);

  // derive dependent variable
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector  saturation(pressure);
  RPK->DeriveSaturationFromPressure(pressure, saturation);

  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(pressure, 0, "pressure");
  GMV::write_cell_data(saturation, 0, "saturation");
  GMV::write_cell_data(*(*FS->permeability())(0), 0, "permeability_x");
  GMV::write_cell_data(*(*FS->permeability())(1), 0, "permeability_y");
  GMV::close_data_file();

  // check the pressure profile
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) CHECK(pressure[c] > 4500.0 && pressure[c] < 101325.0);
}
