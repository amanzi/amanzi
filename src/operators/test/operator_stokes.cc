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
#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "LinearOperatorFactory.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "Accumulation.hh"
#include "AdvectionRiemann.hh"
#include "Elasticity.hh"
#include "TreeOperator.hh"

#include "AnalyticElasticity02.hh"
#include "Verification.hh"

/* *****************************************************************
* Stokes model: exactness test.
***************************************************************** */
TEST(OPERATOR_STOKES_EXACTNESS) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D Stokes: exactness test" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_stokes.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 16, 20, Teuchos::null);

  // modify diffusion coefficient
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  AnalyticElasticity02 ana(mesh);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  // -- on faces
  Teuchos::RCP<BCs> bcf = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, SCHEMA_DOFS_SCALAR));
  std::vector<int>& bcf_model = bcf->bc_model();
  std::vector<double>& bcf_value = bcf->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    const Point& normal = mesh->face_normal(f);
    double area = mesh->face_area(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bcf_model[f] = OPERATOR_BC_DIRICHLET;
      bcf_value[f] = (ana.velocity_exact(xf, 0.0) * normal) / area;
    }
  }

  // -- at nodes
  Point xv(2);
  Teuchos::RCP<BCs> bcv = Teuchos::rcp(new BCs(mesh, AmanziMesh::NODE, SCHEMA_DOFS_SCALAR));
  std::vector<int>& bcv_model = bcv->bc_model();
  std::vector<Point>& bcv_value = bcv->bc_value_point();

  for (int v = 0; v < nnodes_wghost; ++v) {
    mesh->node_get_coordinates(v, &xv);

    if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
        fabs(xv[1]) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bcv_model[v] = OPERATOR_BC_DIRICHLET;
      bcv_value[v] = ana.velocity_exact(xv, 0.0);
    }
  }

  // create elasticity operator 
  Teuchos::ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("elasticity operator");
  Teuchos::RCP<Elasticity> op00 = Teuchos::rcp(new Elasticity(op_list, mesh));
  op00->SetTensorCoefficient(K);
  op00->SetBCs(bcf, bcf);
  op00->AddBCs(bcv, bcv);

  // create divergence operator
  op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("divergence operator");
  Teuchos::RCP<AdvectionRiemann> op10 = Teuchos::rcp(new AdvectionRiemann(op_list, mesh));
  op10->SetBCs(bcf, bcf);
  op10->AddBCs(bcv, bcv);

  // create pressure block (for preconditioner)
  Teuchos::RCP<Accumulation> pc11 = Teuchos::rcp(new Accumulation(AmanziMesh::CELL, mesh));
  Teuchos::RCP<Operator> global11 = pc11->global_operator();

  // create tree operator
  Teuchos::RCP<Operator> global00 = op00->global_operator();
  Teuchos::RCP<Operator> global10 = op10->global_operator();

  const CompositeVectorSpace& cvs = global00->DomainMap();
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(cvs))));
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(op10->global_operator()->RangeMap()))));

  Teuchos::RCP<TreeOperator> op = Teuchos::rcp(new Operators::TreeOperator(tvs));
  op->SetOperatorBlock(0, 0, global00);
  op->SetOperatorBlock(1, 0, global10);
  op->SetOperatorBlock(0, 1, global10, true);

  // create and initialize state variables.
  TreeVector solution(*tvs);
  solution.PutScalar(0.0);

  // create source (for velocity block only)
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("node");

  for (int v = 0; v < nnodes; v++) {
    mesh->node_get_coordinates(v, &xv);
    Point tmp(ana.source_exact(xv, 0.0));
    for (int k = 0; k < 2; ++k) src[k][v] = tmp[k];
  }

  // populate the elasticity block
  op00->UpdateMatrices();
  global00->UpdateRHS(source, true);
  op00->ApplyBCs(true, true);
  global00->SymbolicAssembleMatrix();
  global00->AssembleMatrix();

  // populate the divergence block
  op10->UpdateMatrices(*solution.SubVector(0)->Data());
  op10->ApplyBCs(false, true);

  // populate pressure block (for preconditioner)
  CompositeVector vol(global11->DomainMap());
  vol.PutScalar(1.0);
  pc11->AddAccumulationTerm(vol, 1.0, "cell");
  global11->SymbolicAssembleMatrix();
  global11->AssembleMatrix();

  // create preconditoner, identity is the default one. 
  Teuchos::ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
  global00->InitPreconditioner("Hypre AMG", slist);

  global11->InitPreconditioner("Diagonal", slist); 

  Teuchos::RCP<TreeOperator> pc = Teuchos::rcp(new Operators::TreeOperator(tvs));
  pc->SetOperatorBlock(0, 0, op00->global_operator());
  pc->SetOperatorBlock(1, 1, pc11->global_operator());
  pc->InitBlockDiagonalPreconditioner();

  // Test SPD properties of the matrix and preconditioner.
  VerificationTV ver1(op), ver2(pc);
  ver1.CheckMatrixSPD();
  ver2.CheckPreconditionerSPD();

  // solve the problem
  Teuchos::ParameterList lop_list = plist.get<Teuchos::ParameterList>("solvers");
  AmanziSolvers::LinearOperatorFactory<TreeOperator, TreeVector, TreeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<TreeOperator, TreeVector, TreeVectorSpace> >
     solver = factory.Create("GMRES", lop_list, op, pc);

  TreeVector rhs(*tvs);
  *rhs.SubVector(0)->Data() = *global00->rhs();
  *rhs.SubVector(1)->Data() = *global10->rhs();

  int ierr = solver->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "elasticity solver (" << solver->name() 
              << "): ||r||=" << solver->residual() << " itr=" << solver->num_itrs()
              << " code=" << solver->returned_code() << std::endl;
  }

  // compute velocity error
  double unorm, ul2_err, uinf_err;
  ana.ComputeNodeError(*solution.SubVector(0)->Data(), 0.0, unorm, ul2_err, uinf_err);

  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(*solution.SubVector(1)->Data(), 0.0, pnorm, pl2_err, pinf_err);

  if (MyPID == 0) {
    ul2_err /= unorm;
    printf("L2(u)=%12.8g  Inf(u)=%12.8g  L2(p)=%12.8g  Inf(p)=%12.8g  itr=%3d\n",
        ul2_err, uinf_err, pl2_err, pinf_err, solver->num_itrs());

    CHECK(ul2_err < 0.01);
    CHECK(pl2_err < 0.05);
    CHECK(solver->num_itrs() < 60);
  }

  if (MyPID == 0) {
    const Epetra_MultiVector& u = *solution.SubVector(0)->Data()->ViewComponent("node");
    GMV::open_data_file(*mesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_node_data(u, 0, "velocity_x");
    GMV::write_node_data(u, 1, "velocity_y");
    GMV::close_data_file();
  }
}


