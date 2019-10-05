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
//#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
//#include "Teuchos_XMLParameterListHelpers.hpp"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "GMVMesh.hh"
#include "TreeVector.hh"

#include "State.hh"
#include "Pressure_PK.hh"
#include "OperatorDefs.hh"

/* **************************************************************** */

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

TEST(MULTIPHASE_PRESSURE_FUNCTIONAL) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: Pressure_PK::Functional, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/test_pressure_pk1.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list, "PK Tree");

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 40, 40, gm);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  /* create a simple state populate it */
  RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  // create the global solution vector
  Teuchos::RCP<TreeVector> soln_ = Teuchos::rcp(new TreeVector());

  // Create pressure PK
  Multiphase::Pressure_PK* PPK = new Multiphase::Pressure_PK(*pk_tree_list, parameter_list, S, soln_);
  PPK->Setup();
  S->Setup();
  S->InitializeFields();

  /* modify the default state for the problem at hand */
  std::string passwd("state");
  const Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell", false);
  for (int c = 0; c < K.MyLength(); c++) {
    K[0][c] = 1.0;
    K[1][c] = 1.0;
  }

  const Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);
  Teuchos::RCP<CompositeVector> uold = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> unew = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> functional = Teuchos::rcp(new CompositeVector(*cvs));
  uold->PutScalar(0.0);
  unew->PutScalar(0.0);
  functional->PutScalar(1.0);

  Teuchos::RCP<TreeVector> u_old = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> u_new = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f = Teuchos::rcp(new TreeVector());
  u_old->SetData(uold);
  u_new->SetData(unew);
  f->SetData(functional);

  Epetra_MultiVector& u_new_cell = *unew->ViewComponent("cell", true);
  for (int c = 0; c < u_new_cell.MyLength(); c++) 
  {
    const Point& xc = mesh->cell_centroid(c);
    u_new_cell[0][c] = u_exact(xc) ;
  }

  // Init the PK and solve for pressure
  PPK->Initialize();
  S->CheckAllFieldsInitialized();
  PPK->Functional(0.0, 1.0, u_old, u_new, f);
  
  const Teuchos::RCP<CompositeVector>& residual = f->Data();
  Epetra_MultiVector residual_cell = *residual->ViewComponent("cell");
  for (int i = 0; i < residual_cell.MyLength(); i++) {
    CHECK_CLOSE(0.0, residual_cell[0][i], 0.0001);
  }

  // get the data for viz
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_pressure_pk_functional.gmv");
    GMV::start_data();
    GMV::write_cell_data(residual_cell, 0, "residual");
    GMV::close_data_file();
  }

  delete PPK;
}
