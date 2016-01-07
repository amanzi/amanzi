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
#include "LinearOperatorFactory.hh"
#include "mfd3d_diffusion.hh"
#include "Tensor.hh"

// Operiators
#include "Analytic00.hh"
#include "Analytic01.hh"
#include "BCs.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorDiffusionMFD.hh"
#include "TreeOperator.hh"


/* *****************************************************************
* This test runs one loop of the convergence test with two
* uncoupled diffusion problems on the diagonal of the tree operator.
****************************************************************** */
TEST(OPERATOR_UNCOUPLED) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "Test: 2D block uncoupled system of elliptic opeartors" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_tree.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create a mesh 
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 80, 80, gm);
  Teuchos::RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  // create diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  Analytic01 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);
  std::vector<double> bc_mixed;

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_value[f] = ana.pressure_exact(xf, 0.0);
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create source 
  CompositeVector source(*cvs);
  source.PutScalarMasterAndGhosted(0.0);
  
  Epetra_MultiVector& src = *source.ViewComponent("cell");
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    src[0][c] += ana.source_exact(xc, 0.0);
  }

  // populate the diffusion operator
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operators")
                                      .get<Teuchos::ParameterList>("mixed diffusion p-lambda");
  Teuchos::RCP<OperatorDiffusionMFD> op = Teuchos::rcp(new OperatorDiffusionMFD(olist, mesh));
  op->SetBCs(bc, bc);

  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // get and assemble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, false);
  op->ApplyBCs(true, true);

  // create the TreeOperator, which combines two copies of op1 in a block-diagonal operator.
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  Teuchos::RCP<TreeVectorSpace> cvs_as_tv = Teuchos::rcp(new TreeVectorSpace(cvs));
  tvs->PushBack(cvs_as_tv);
  tvs->PushBack(cvs_as_tv);
  Teuchos::RCP<TreeOperator> tree_op = Teuchos::rcp(new TreeOperator(tvs));

  tree_op->SetOperatorBlock(0, 0, global_op);
  tree_op->SetOperatorBlock(1, 1, global_op);

  // assemble the tree operator
  tree_op->SymbolicAssembleMatrix();
  tree_op->AssembleMatrix();
  
  // create preconditioner
  Teuchos::ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  tree_op->InitPreconditioner("Hypre AMG", slist);
  
  // solve the problem
  Teuchos::ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
  AmanziSolvers::LinearOperatorFactory<TreeOperator, TreeVector, TreeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<TreeOperator, TreeVector, TreeVectorSpace> >
      solver = factory.Create("AztecOO CG", lop_list, tree_op);

  // create a solution vector, rhs vector
  Teuchos::RCP<CompositeVector> rhs_comp = Teuchos::rcp(new CompositeVector(*global_op->rhs()));
  Teuchos::RCP<TreeVector> rhs_comp_tv = Teuchos::rcp(new TreeVector());
  rhs_comp_tv->SetData(rhs_comp);
  TreeVector rhs;
  rhs.PushBack(rhs_comp_tv);
  rhs.PushBack(rhs_comp_tv);

  TreeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver->ApplyInverse(rhs, solution);

  // calculate pressure errors
  Epetra_MultiVector& pA = *solution.SubVector(0)->Data()->ViewComponent("cell", false);
  double pnormA, pl2_errA, pinf_errA;
  ana.ComputeCellError(pA, 0.0, pnormA, pl2_errA, pinf_errA);

  Epetra_MultiVector& pB = *solution.SubVector(1)->Data()->ViewComponent("cell", false);
  double pnormB, pl2_errB, pinf_errB;
  ana.ComputeCellError(pB, 0.0, pnormB, pl2_errB, pinf_errB);
  
  // calculate flux errors
  CompositeVector flux(*cvs);
  Epetra_MultiVector& flx = *flux.ViewComponent("face", true);
  double unorm, ul2_err, uinf_err;

  op->UpdateFlux(*solution.SubVector(0)->Data(), flux);
  flux.ScatterMasterToGhosted();
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_errA /= pnormA;
    pl2_errB /= pnormB;
    ul2_err /= unorm;
    printf("L2(p)=%9.6f,%9.6f  Inf(p)=%9.6f,%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f itr=%3d\n", 
           pl2_errA, pl2_errB, pinf_errA, pinf_errB, ul2_err, uinf_err, solver->num_itrs()); 

    CHECK(pl2_errA < 0.15);
    CHECK(pl2_errB < 0.15);    
    CHECK(ul2_err < 0.15);
  }
}


/* *****************************************************************
* This test runs one loop of the convergence test with two
* coupled diffusion problems. Each block of the super matrix is
* a diffusion operator.
****************************************************************** */
void RunTest(std::string op_list_name) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "Test: 2D block coupled system of elliptic operators" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_tree.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create a mesh 
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Mesh> mesh;
  if (op_list_name == "mixed diffusion p-lambda") {
    mesh = meshfactory("test/median32x33.exo", gm);
  } else {
    mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 7, 5, gm);
  }

  // create boundary data
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);
  std::vector<double> bc_mixed;

  Analytic00 ana(mesh, 1.0, 2.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_value[f] = ana.pressure_exact(xf, 0.0);
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create discretization space
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create problem data: coefficients k1, k2, and zero source.
  Teuchos::RCP<CompositeVector> k1 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> k2 = Teuchos::rcp(new CompositeVector(*cvs));
  k1->PutScalar(1.0);
  k2->PutScalar(0.5);
  
  // populate the diagonal Laplace operators
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operators")
                                      .get<Teuchos::ParameterList>(op_list_name);
  Operators::OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op00 = opfactory.Create(olist, mesh, bc);
  op00->SetScalarCoefficient(k1, Teuchos::null);
  op00->UpdateMatrices(Teuchos::null, Teuchos::null);

  Teuchos::RCP<OperatorDiffusion> op11 = opfactory.Create(olist, mesh, bc);
  op11->SetScalarCoefficient(k1, Teuchos::null);
  op11->UpdateMatrices(Teuchos::null, Teuchos::null);

  // populate the off-diagonal Laplace operators
  Teuchos::RCP<OperatorDiffusion> op01 = opfactory.Create(olist, mesh, bc);
  op01->SetScalarCoefficient(k2, Teuchos::null);
  op01->UpdateMatrices(Teuchos::null, Teuchos::null);

  Teuchos::RCP<OperatorDiffusion> op10 = opfactory.Create(olist, mesh, bc);
  op10->SetScalarCoefficient(k2, Teuchos::null);
  op10->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update right-hand side (ZERO) and apply boundary conditions
  op00->ApplyBCs(true, true);
  op11->ApplyBCs(true, true);
  op01->ApplyBCs(false, false);
  op10->ApplyBCs(false, true);

  // create the TreeOperator, which combines four operators in a 2x2 block operator.
  Teuchos::RCP<const CompositeVectorSpace> cvs2 = Teuchos::rcpFromRef(op00->global_operator()->DomainMap());
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  Teuchos::RCP<TreeVectorSpace> cvs_as_tv = Teuchos::rcp(new TreeVectorSpace(cvs2));
  tvs->PushBack(cvs_as_tv);
  tvs->PushBack(cvs_as_tv);
  Teuchos::RCP<TreeOperator> tree_op = Teuchos::rcp(new TreeOperator(tvs));

  tree_op->SetOperatorBlock(0, 0, op00->global_operator());
  tree_op->SetOperatorBlock(1, 1, op11->global_operator());
  tree_op->SetOperatorBlock(0, 1, op01->global_operator());
  tree_op->SetOperatorBlock(1, 0, op10->global_operator());

  // assemble the tree operator
  tree_op->SymbolicAssembleMatrix();
  tree_op->AssembleMatrix();
  
  // create preconditioner
  Teuchos::ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  tree_op->InitPreconditioner("Hypre AMG", slist);
  
  // solve the problem
  // -- create an iterative solver
  Teuchos::ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
  AmanziSolvers::LinearOperatorFactory<TreeOperator, TreeVector, TreeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<TreeOperator, TreeVector, TreeVectorSpace> >
      solver = factory.Create("AztecOO CG", lop_list, tree_op);

  // -- create a rhs vector
  TreeVector rhs;
  Teuchos::RCP<CompositeVector> rhs_cv0, rhs_cv1;
  Teuchos::RCP<TreeVector> rhs_tv0 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> rhs_tv1 = Teuchos::rcp(new TreeVector());

  rhs_cv0 = op00->global_operator()->rhs();
  rhs_cv1 = op01->global_operator()->rhs();
  rhs_cv0->Update(1.0, *rhs_cv1, 1.0);
  rhs_tv0->SetData(rhs_cv0);
  rhs.PushBack(rhs_tv0);

  rhs_cv0 = op10->global_operator()->rhs();
  rhs_cv1 = op11->global_operator()->rhs();
  rhs_cv1->Update(1.0, *rhs_cv0, 1.0);
  rhs_tv1->SetData(rhs_cv1);
  rhs.PushBack(rhs_tv1);

  // -- run iterative solver
  TreeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver->ApplyInverse(rhs, solution);

  // calculate solution errors
  Epetra_MultiVector& pA = *solution.SubVector(0)->Data()->ViewComponent("cell", false);
  double pnormA, pl2_errA, pinf_errA;
  ana.ComputeCellError(pA, 0.0, pnormA, pl2_errA, pinf_errA);

  Epetra_MultiVector& pB = *solution.SubVector(1)->Data()->ViewComponent("cell", false);
  double pnormB, pl2_errB, pinf_errB;
  ana.ComputeCellError(pB, 0.0, pnormB, pl2_errB, pinf_errB);
  
  // calculate flux errors
  CompositeVector flux(*cvs);
  Epetra_MultiVector& flx = *flux.ViewComponent("face", true);
  double unorm, ul2_err, uinf_err;

  op00->UpdateFlux(*solution.SubVector(0)->Data(), flux);
  flux.ScatterMasterToGhosted();
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_errA /= pnormA;
    pl2_errB /= pnormB;
    ul2_err /= unorm;
    printf("L2(p)=%9.6f,%9.6f  Inf(p)=%9.6f,%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f itr=%3d\n", 
           pl2_errA, pl2_errB, pinf_errA, pinf_errB, ul2_err, uinf_err, solver->num_itrs()); 

    CHECK(pl2_errA < 1e-8);
    CHECK(pl2_errB < 1e-8);    
    CHECK(ul2_err < 1e-8);
  }
}


/* *****************************************************************
* p-lambda discretization
****************************************************************** */
TEST(COUPLED_SYSTEM_PLAMDA) {
  RunTest("mixed diffusion p-lambda");
}


/* *****************************************************************
* cell-centered discretization
****************************************************************** */
TEST(COUPLED_SYSTEM_CELL_CENTERED) {
  RunTest("mixed diffusion cell-centered");
}
