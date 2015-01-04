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

#include "EnergyTwoPhase_PK.hh"
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "OperatorDiffusionFactory.hh"
#include "State.hh"

/* **************************************************************** */
TEST(ENERGY_2D_MATRIX) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  using namespace Amanzi::Energy;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: 2D thermal, homogeneous medium" << std::endl;

  // read parameter list 
  std::string xmlFileName = "test/energy_matrix_2D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::RCP<const Teuchos::ParameterList> plist = 
      Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 3, 3, gm);

  // create a simple state and populate it
  Amanzi::VerboseObject::hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->sublist("State");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  // initialize the Energy process kernel 
  EnergyTwoPhase_PK* EPK = new EnergyTwoPhase_PK(plist, S);
  EPK->Setup();
std::cout << "Passes EPK.Setup()" << std::endl;
  S->Setup();
std::cout << "Passed S.Setup()" << std::endl;
  S->InitializeFields();
std::cout << "Passed S.InitilizeFields()" << std::endl;
  EPK->InitializeFields();
std::cout << "Passed EPK.InitilizeField()" << std::endl;
  S->InitializeEvaluators();
std::cout << "Passed S.InitilizeEvaluators()" << std::endl;
  S->WriteDependencyGraph();
  S->CheckAllFieldsInitialized();
std::cout << "DONE" << std::endl;

  // modify the default state for the problem at hand 
  // create the initial temperature function 
  std::string passwd("thermal");
  Epetra_MultiVector& temperature = *S->GetFieldData("temperature", passwd)->ViewComponent("cell");
  temperature.PutScalar(273.0);

  // compute conductivity
  EPK->UpdateConductivityData(S.ptr());

  // create boundary data
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);
  std::vector<double> bc_mixed;
  
  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = 0.0;
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));
  
  // create diffusion operator 
  const Teuchos::ParameterList& elist = plist->sublist("Energy");
  Teuchos::ParameterList oplist = elist.sublist("operators")
                                       .sublist("diffusion operator")
                                       .sublist("preconditioner");
  OperatorDiffusionFactory opfactory;
  AmanziGeometry::Point g(2);
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc, oplist, g, 0);

  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_FACE + Operators::OPERATOR_SCHEMA_DOFS_CELL;
  double rho(1.0), mu(1.0);
  op->Setup(EPK->get_K(), Teuchos::null, Teuchos::null, rho, mu);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs();
  op->SymbolicAssembleMatrix(schema);
  op->AssembleMatrix(schema);

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"energy.gmv");
    GMV::start_data();
    GMV::write_cell_data(temperature, 0, "temperature");
    GMV::close_data_file();
  }

  delete EPK;
}
