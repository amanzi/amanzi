/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
#include "Tensor.hh"

// Amanzi::Operators
#include "Operator_FaceCell.hh"
#include "OperatorDefs.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_DiffusionMFD.hh"


/* *****************************************************************
 * This test verified that operators can be computed in arbitrary
 * order. In addition, the factory of operators is used.
 ***************************************************************** */
TEST(ADVECTION_DIFFUSION_COMMUTE)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0)
    std::cout << "\nTest: Commuting of advection-duffusion operators."
              << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_commute.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 40, 40);

  // modify diffusion coefficient.
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned =
    mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost =
    mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc =
    Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
        fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = xf[1] * xf[1];
    }
  }

  // create the global operator space
  Teuchos::RCP<CompositeVectorSpace> cvs =
    Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create velocity field
  Teuchos::RCP<CompositeVector> u = Teuchos::rcp(new CompositeVector(*cvs));
  Epetra_MultiVector& uf = *u->ViewComponent("face");
  int nfaces =
    mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  Point vel(4.0, 4.0);
  for (int f = 0; f < nfaces; f++) { uf[0][f] = vel * mesh->face_normal(f); }

  // create the global op
  Teuchos::ParameterList plist1;
  Teuchos::RCP<Operator> global_op =
    Teuchos::rcp(new Operator_FaceCell(cvs, plist1));

  // create advection operator
  Teuchos::ParameterList alist;
  Teuchos::RCP<PDE_AdvectionUpwind> op1 =
    Teuchos::rcp(new PDE_AdvectionUpwind(alist, global_op));
  op1->SetBCs(bc, bc);
  op1->Setup(*u);
  op1->UpdateMatrices(u.ptr());

  // add the diffusion operator
  Teuchos::ParameterList olist =
    plist.sublist("PK operator").sublist("diffusion operator mfd");
  auto op2 = Teuchos::rcp(new PDE_DiffusionMFD(olist, global_op));
  op2->Init(olist);
  op2->SetBCs(bc, bc);
  op2->Setup(K, Teuchos::null, Teuchos::null);
  op2->UpdateMatrices(Teuchos::null, Teuchos::null);

  // create a preconditioner
  op1->ApplyBCs(true, true, true);
  op2->ApplyBCs(true, true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // make reverse assembling: diffusion + advection
  Teuchos::ParameterList plist2;
  Teuchos::RCP<Operator> global_op2 =
    Teuchos::rcp(new Operator_FaceCell(cvs, plist2));

  Teuchos::ParameterList olist2 =
    plist.sublist("PK operator").sublist("diffusion operator mfd");
  auto op3 = Teuchos::rcp(new PDE_DiffusionMFD(olist2, global_op2));
  op3->Init(olist2);
  op3->SetBCs(bc, bc);
  op3->Setup(K, Teuchos::null, Teuchos::null);
  op3->UpdateMatrices(Teuchos::null, Teuchos::null);

  Teuchos::ParameterList alist2;
  Teuchos::RCP<PDE_AdvectionUpwind> op4 =
    Teuchos::rcp(new PDE_AdvectionUpwind(alist2, global_op2));
  op4->SetBCs(bc, bc);
  op4->Setup(*u);
  op4->UpdateMatrices(u.ptr());

  // create a preconditioner
  op3->ApplyBCs(true, true, true);
  op4->ApplyBCs(true, true, true);
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


TEST(ADVECTION_DIFFUSION_COMMUTE_FV)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0)
    std::cout << "\nTest: Commuting of advection-duffusion operators."
              << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_commute.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 40, 40);

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned =
    mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost =
    mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc =
    Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
        fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = xf[1] * xf[1];
    }
  }

  // create the global operator space
  Teuchos::RCP<CompositeVectorSpace> cvs =
    Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create velocity field
  Teuchos::RCP<CompositeVector> u = Teuchos::rcp(new CompositeVector(*cvs));
  Epetra_MultiVector& uf = *u->ViewComponent("face");
  int nfaces =
    mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  Point vel(4.0, 4.0);
  for (int f = 0; f < nfaces; f++) { uf[0][f] = vel * mesh->face_normal(f); }

  // create the global op
  Teuchos::ParameterList plist1;
  Teuchos::RCP<Operator> global_op =
    Teuchos::rcp(new Operator_Cell(cvs, plist1, OPERATOR_SCHEMA_DOFS_CELL));

  // create advection operator
  Teuchos::ParameterList alist;
  Teuchos::RCP<PDE_AdvectionUpwind> op1 =
    Teuchos::rcp(new PDE_AdvectionUpwind(alist, global_op));
  op1->SetBCs(bc, bc);
  op1->Setup(*u);
  op1->UpdateMatrices(u.ptr());

  // add the diffusion operator
  Teuchos::ParameterList olist =
    plist.sublist("PK operator").sublist("diffusion operator fv");
  Teuchos::RCP<PDE_Diffusion> op2 =
    Teuchos::rcp(new PDE_DiffusionFV(olist, global_op));
  op2->SetBCs(bc, bc);
  op2->Setup(K, Teuchos::null, Teuchos::null);
  op2->UpdateMatrices(Teuchos::null, Teuchos::null);

  // create a preconditioner
  op1->ApplyBCs(true, true, true);
  op2->ApplyBCs(false, true, false);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // make reverse assembling: diffusion + advection
  Teuchos::ParameterList plist2;
  Teuchos::RCP<Operator> global_op2 =
    Teuchos::rcp(new Operator_Cell(cvs, plist2, OPERATOR_SCHEMA_DOFS_CELL));

  Teuchos::ParameterList olist2 =
    plist.sublist("PK operator").sublist("diffusion operator fv");
  Teuchos::RCP<PDE_Diffusion> op3 =
    Teuchos::rcp(new PDE_DiffusionFV(olist2, global_op2));
  op3->SetBCs(bc, bc);
  op3->Setup(K, Teuchos::null, Teuchos::null);
  op3->UpdateMatrices(Teuchos::null, Teuchos::null);

  Teuchos::ParameterList alist2;
  Teuchos::RCP<PDE_AdvectionUpwind> op4 =
    Teuchos::rcp(new PDE_AdvectionUpwind(alist2, global_op2));
  op4->SetBCs(bc, bc);
  op4->Setup(*u);
  op4->UpdateMatrices(u.ptr());

  // create a preconditioner
  op3->ApplyBCs(false, true, false);
  op4->ApplyBCs(true, true, true);
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
