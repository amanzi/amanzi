/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
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

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "mfd3d_diffusion.hh"
#include "LinearOperatorFactory.hh"
#include "Tensor.hh"

// Operators
#include "Operator_FaceCell.hh"
#include "OperatorAdvection.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionMFD.hh"
#include "OperatorDiffusionFV.hh"
#include "OperatorDiffusionFactory.hh"

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
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 40, 40, gm);

  // modify diffusion coefficient.
  // -- since rho=mu=1.0, we do not need additional special steps.
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
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

  // create the global operator space
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create velocity field
  CompositeVector u(*cvs);
  Epetra_MultiVector& uf = *u.ViewComponent("face");
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  Point vel(4.0, 4.0);
  for (int f = 0; f < nfaces; f++) {
    uf[0][f] = vel * mesh->face_normal(f);
  }

  // create the global op
  Teuchos::ParameterList plist1;
  Teuchos::RCP<Operator> global_op = Teuchos::rcp(new Operator_FaceCell(cvs, plist1));
  
  // create advection operator
  Teuchos::ParameterList alist;
  Teuchos::RCP<OperatorAdvection> op1 = Teuchos::rcp(new OperatorAdvection(alist, global_op));
  op1->Setup(u);
  op1->UpdateMatrices(u);

  // add the diffusion operator
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                      .get<Teuchos::ParameterList>("diffusion operator mfd");
  Teuchos::RCP<OperatorDiffusion> op2 = Teuchos::rcp(new OperatorDiffusionMFD(olist, global_op));
  op2->SetBCs(bc, bc);
  op2->Setup(K, Teuchos::null, Teuchos::null);
  op2->UpdateMatrices(Teuchos::null, Teuchos::null);

  // create a preconditioner
  op1->ApplyBCs(bc, true);
  op2->ApplyBCs(true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // make reverse assembling: diffusion + advection
  Teuchos::ParameterList plist2;
  Teuchos::RCP<Operator> global_op2 = Teuchos::rcp(new Operator_FaceCell(cvs, plist2));
  
  Teuchos::ParameterList olist2 = plist.get<Teuchos::ParameterList>("PK operator")
                                       .get<Teuchos::ParameterList>("diffusion operator mfd");
  Teuchos::RCP<OperatorDiffusion> op3 = Teuchos::rcp(new OperatorDiffusionMFD(olist2, global_op2));
  op3->SetBCs(bc, bc);
  op3->Setup(K, Teuchos::null, Teuchos::null);
  op3->UpdateMatrices(Teuchos::null, Teuchos::null);

  Teuchos::ParameterList alist2;
  Teuchos::RCP<OperatorAdvection> op4 = Teuchos::rcp(new OperatorAdvection(alist2, global_op2));
  op4->Setup(u);
  op4->UpdateMatrices(u);

  // create a preconditioner
  op3->ApplyBCs(true, true);
  op4->ApplyBCs(bc, true);
  global_op2->SymbolicAssembleMatrix();
  global_op2->AssembleMatrix();

  // compare matrices
  int n, nrows = global_op->A()->NumMyRows();

  for (int i = 0; i < nrows; ++i) {
    double *val2, *val4; 
    global_op->A()->ExtractMyRowView(i, n, val2); 
    global_op2->A()->ExtractMyRowView(i, n, val4); 
    for (int k = 0; k < n; ++k) CHECK_CLOSE(val2[k], val4[k], 1e-10);
  }
}


TEST(ADVECTION_DIFFUSION_COMMUTE_FV) {
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
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 40, 40, gm);

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do need additional special steps.
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
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

  // create the global operator space
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create velocity field
  CompositeVector u(*cvs);
  Epetra_MultiVector& uf = *u.ViewComponent("face");
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  Point vel(4.0, 4.0);
  for (int f = 0; f < nfaces; f++) {
    uf[0][f] = vel * mesh->face_normal(f);
  }

  // create the global op
  Teuchos::ParameterList plist1;
  Teuchos::RCP<Operator> global_op = Teuchos::rcp(new Operator_Cell(cvs, plist1, OPERATOR_SCHEMA_DOFS_CELL));
  
  // create advection operator
  Teuchos::ParameterList alist;
  Teuchos::RCP<OperatorAdvection> op1 = Teuchos::rcp(new OperatorAdvection(alist, global_op));
  op1->Setup(u);
  op1->UpdateMatrices(u);

  // add the diffusion operator
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                      .get<Teuchos::ParameterList>("diffusion operator fv");
  Teuchos::RCP<OperatorDiffusion> op2 = Teuchos::rcp(new OperatorDiffusionFV(olist, global_op));
  op2->SetBCs(bc, bc);
  op2->Setup(K, Teuchos::null, Teuchos::null);
  op2->UpdateMatrices(Teuchos::null, Teuchos::null);

  // create a preconditioner
  op1->ApplyBCs(bc, true);
  op2->ApplyBCs(false, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // make reverse assembling: diffusion + advection
  Teuchos::ParameterList plist2;
  Teuchos::RCP<Operator> global_op2 = Teuchos::rcp(new Operator_Cell(cvs, plist2, OPERATOR_SCHEMA_DOFS_CELL));
  
  Teuchos::ParameterList olist2 = plist.get<Teuchos::ParameterList>("PK operator")
                                      .get<Teuchos::ParameterList>("diffusion operator fv");
  Teuchos::RCP<OperatorDiffusion> op3 = Teuchos::rcp(new OperatorDiffusionFV(olist2, global_op2));
  op3->SetBCs(bc, bc);
  op3->Setup(K, Teuchos::null, Teuchos::null);
  op3->UpdateMatrices(Teuchos::null, Teuchos::null);

  Teuchos::ParameterList alist2;
  Teuchos::RCP<OperatorAdvection> op4 = Teuchos::rcp(new OperatorAdvection(alist2, global_op2));
  op4->Setup(u);
  op4->UpdateMatrices(u);

  // create a preconditioner
  op3->ApplyBCs(false, true);
  op4->ApplyBCs(bc, true);
  global_op2->SymbolicAssembleMatrix();
  global_op2->AssembleMatrix();

  // compare matrices
  int n, nrows = global_op->A()->NumMyRows();

  for (int i = 0; i < nrows; ++i) {
    double *val2, *val4; 
    global_op->A()->ExtractMyRowView(i, n, val2); 
    global_op2->A()->ExtractMyRowView(i, n, val4); 
    for (int k = 0; k < n; ++k) CHECK_CLOSE(val2[k], val4[k], 1e-10);
  }
}
