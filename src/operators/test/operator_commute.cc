/*
  This is the operators component of the Amanzi code. 

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
#include "mfd3d_diffusion.hh"
#include "tensor.hh"

#include "LinearOperatorFactory.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorAdvection.hh"


/* *****************************************************************
* This test verified that operators can be computed in arbitrary
* order. In addition, the factory of operators in used.
* **************************************************************** */
TEST(ADVECTION_DIFFUSION_COMMUTE) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Commuting of advection-duffusion operators." << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_commute.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 40, 40, gm);

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);

  // create boundary data
  std::vector<int> bc_model(nfaces_wghost, OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);
  std::vector<double> bc_mixed;

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = xf[1] * xf[1];
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create operator map 
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create velocity field
  CompositeVector u(*cvs);
  Epetra_MultiVector& uf = *u.ViewComponent("face");
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  Point vel(4.0, 4.0);
  for (int f = 0; f < nfaces; f++) {
    uf[0][f] = vel * mesh->face_normal(f);
  }

  // create advection operator
  Teuchos::RCP<OperatorAdvection> op1 = Teuchos::rcp(new OperatorAdvection(cvs, 0));
  op1->Init();
  op1->Setup(u);
  op1->UpdateMatrices(u);

  // add the diffusion operator. It is the last due to BCs.
  int schema_base = Operators::OPERATOR_SCHEMA_BASE_CELL;
  int schema_dofs = Operators::OPERATOR_SCHEMA_DOFS_FACE + 
                    Operators::OPERATOR_SCHEMA_DOFS_CELL;
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                      .get<Teuchos::ParameterList>("diffusion operator");

  Teuchos::RCP<OperatorDiffusion> op2 = Teuchos::rcp(new OperatorDiffusion(*op1, olist, bc));
  op2->Setup(K, Teuchos::null, Teuchos::null, rho, mu);
  op2->UpdateMatrices(Teuchos::null, Teuchos::null);
  op2->ApplyBCs();

  // create a preconditioner 
  op2->SymbolicAssembleMatrix(schema_dofs);
  op2->AssembleMatrix(schema_dofs);

  // make reverse assembling: diffusion + advection
  Teuchos::RCP<OperatorDiffusion> op3 = Teuchos::rcp(new OperatorDiffusion(cvs, olist, bc));
  op3->Init();
  op3->Setup(K, Teuchos::null, Teuchos::null, rho, mu);
  op3->UpdateMatrices(Teuchos::null, Teuchos::null);

  Teuchos::RCP<OperatorAdvection> op4 = Teuchos::rcp(new OperatorAdvection(*op3));
  op4->Setup(u);
  op4->UpdateMatrices(u);
  op4->ApplyBCs();

  // create a preconditioner 
  op4->SymbolicAssembleMatrix(schema_dofs);
  op4->AssembleMatrix(schema_dofs); 

  // compare matrices
  int n, nrows = op2->A()->NumMyRows();

  for (int i = 0; i < nrows; ++i) {
    double *val2, *val4; 
    op2->A()->ExtractMyRowView(i, n, val2); 
    op4->A()->ExtractMyRowView(i, n, val4); 
    for (int k = 0; k < n; ++k) CHECK_CLOSE(val2[k], val4[k], 1e-10);
  }
}
