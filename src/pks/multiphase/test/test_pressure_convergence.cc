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
#include "Mesh.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Pressure_PK.hh"
#include "OperatorDefs.hh"

// Some analytic data
// P(x,y) = sin(pi*x) * sin(2*pi*y)
double u_exact(const Amanzi::AmanziGeometry::Point& p) {
  double x = p[0];
  double y = p[1];
  return sin(M_PI*x) * sin(2*M_PI*y);
}

// Source S = 5*pi^2*sin(pi*x)*sin(2*pi*y)
double s_exact(const Amanzi::AmanziGeometry::Point& p) {
  double x = p[0];
  double y = p[1];
  return 5 * M_PI * M_PI * sin(M_PI*x) * sin(2*M_PI*y);
}

/* **************************************************************** */
TEST(MULTIPHASE_PRESSURE) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

std::vector<double> L2errors;
for (int nx = 8; nx < 150; nx *= 2) {
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: steady state 2-phase immiscible flow, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/test_pressure_pk1.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list, "PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, nx, nx, gm);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  /* create a simple state populate it */
  RCP<State> S = rcp(new State(*state_list));
  //RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  const Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);

  const Teuchos::RCP<const CompositeVectorSpace> tvs = Teuchos::rcp(new CompositeVectorSpace(*cvs));
  // create the global solution vector
  Teuchos::RCP<TreeVector> soln_ = Teuchos::rcp(new TreeVector(tvs));

  // Create pressure PK
  Multiphase::Pressure_PK* PPK = new Multiphase::Pressure_PK(*pk_tree_list, parameter_list, S, soln_);
  PPK->Setup();
  S->Setup();
  S->InitializeFields();

  // Init the PK and solve for pressure
  PPK->Initialize();
  S->CheckAllFieldsInitialized();
  PPK->AdvanceStep(0.0, 1.0, false);
  PPK->CommitState(S.ptr());

  std::string passwd("state");
  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);
  double p_error(0.0);
  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    double tmp = u_exact(xc);
    double volume = mesh->cell_volume(c);

    p_error += std::pow(tmp - p[0][c], 2.0) * volume;
  }
  L2errors.push_back(pow(p_error, 0.5));

  // get the data for viz
  Epetra_MultiVector& perm = *S->GetFieldData("permeability", passwd)->ViewComponent("cell", false);
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, "test_pressure_convergence.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::write_cell_data(perm, 0, "permeability");
    GMV::close_data_file();
  }
  delete PPK;
}
  // print out the errors
  for (std::vector<double>::iterator it = L2errors.begin(); it != L2errors.end(); it++) {
    std::cout << "L2Error = " << *it << "\n";
  }
}
