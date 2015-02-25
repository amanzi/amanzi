/*
  This is the operator component of the Amanzi code. 

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
#include "GMVMesh.hh"
#include "LinearOperatorFactory.hh"

#include "tensor.hh"
#include "mfd3d_diffusion.hh"

#include "BCs.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionMFD.hh"

#include "TreeOperator.hh"

#include "Analytic01.hh"

/* *****************************************************************
 * This test simply runs one iteration of the convergence test with two
 * uncoupled diffusion problems on the diagonal of the tree operator.
 * **************************************************************** */
TEST(OPERATOR_MIXED_DIFFUSION) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: 2D steady-state elliptic solver, mixed discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::hide_line_prefix = true;

  // create a mesh 
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 40, 40, gm);
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
  double rho(1.0), mu(1.0);

  // create boundary data
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  Point xv(2);
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

  // MAIN LOOP
  double factor = 1.0;
    
  // populate the diffusion operator
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operators")
      .get<Teuchos::ParameterList>("mixed diffusion");
  Teuchos::RCP<OperatorDiffusionMFD> op = Teuchos::rcp(new OperatorDiffusionMFD(olist, mesh));
  op->SetBCs(bc);

  op->set_factor(factor);  // for developers only
  op->Setup(K, Teuchos::null, Teuchos::null, rho, mu);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // get and assmeble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, false);
  op->ApplyBCs();

  // create the TreeOperator, which combines two copies of op1 in a block-diagonal operator.
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  Teuchos::RCP<TreeVectorSpace> cvs_as_tv = Teuchos::rcp(new TreeVectorSpace(cvs));
  tvs->PushBack(cvs_as_tv);
  tvs->PushBack(cvs_as_tv);
  Teuchos::RCP<TreeOperator> tree_op = Teuchos::rcp(new TreeOperator(tvs));

  tree_op->SetOperatorBlock(0,0,global_op);
  tree_op->SetOperatorBlock(1,1,global_op);

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
    printf("scale=%7.4g  L2(p)=%9.6f,%9.6f  Inf(p)=%9.6f,%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f itr=%3d\n", 
           factor, pl2_errA, pl2_errB, pinf_errA, pinf_errB, ul2_err, uinf_err, solver->num_itrs()); 

    CHECK(pl2_errA < 0.15);
    CHECK(pl2_errB < 0.15);    
    CHECK(ul2_err < 0.15);
  }
}

