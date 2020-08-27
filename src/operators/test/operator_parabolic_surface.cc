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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_DiffusionMFD.hh"

#include "Verification.hh"


/* *****************************************************************
* This test replaces tensor and boundary conditions by continuous
* functions. This is a prototype for future solvers.
* **************************************************************** */
void RunTest(std::string op_list_name) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Singular-perturbed Laplace Beltrami solver" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_laplace_beltrami.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.sublist("Regions Closed");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh = meshfactory.create("test/sphere.exo");

  // extract surface mesh
  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top surface"));
  RCP<Mesh> surfmesh = meshfactory.create(mesh, setnames, AmanziMesh::FACE);

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the diffsuion coefficient.
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = surfmesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
  }

  // create boundary data (no mixed bc)
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(surfmesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  bc->bc_model();  // allocate internal
  bc->bc_value();  // memory

  // create diffusion operator 
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(surfmesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create source and add it to the operator
  CompositeVector source(*cvs);
  source.PutScalarMasterAndGhosted(0.0);
  
  Epetra_MultiVector& src = *source.ViewComponent("cell");
  for (int c = 0; c < 20; c++) {
    if (MyPID == 0) src[0][c] = 1.0;
  }

  // add accumulation terms
  CompositeVector solution(*cvs);
  solution.PutScalar(0.0);  // solution at time T=0

  CompositeVector phi(*cvs);
  phi.PutScalar(0.2);

  double dT = 10.0;

  // add the diffusion operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist(op_list_name);
  PDE_DiffusionMFD op(olist, surfmesh);
  op.Init(olist);
  op.SetBCs(bc, bc);
  op.Setup(K, Teuchos::null, Teuchos::null);
  op.UpdateMatrices(Teuchos::null, Teuchos::null);

  // get the global operator
  Teuchos::RCP<Operator> global_op = op.global_operator();

  // add accumulation terms
  PDE_Accumulation op_acc(AmanziMesh::CELL, global_op);
  op_acc.AddAccumulationDelta(solution, phi, phi, dT, "cell");

  // apply BCs and assemble
  global_op->UpdateRHS(source, false);
  op.ApplyBCs(true, true, true);

  // create preconditoner
  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"), "PCG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD();
  ver.CheckPreconditionerSPD(1e-11);

  CompositeVector rhs = *global_op->rhs();
  solution.PutScalar(0.0);
  global_op->ApplyInverse(rhs, solution);

  // ver.CheckResidual(solution, 1.0e-12);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual() 
              << " itr=" << global_op->num_itrs()
              << " code=" << global_op->returned_code() << std::endl;
  }

  // repeat the above without destroying the operators.
  solution.PutScalar(0.0);
  global_op->rhs()->PutScalar(0.);

  op.UpdateMatrices(Teuchos::null, Teuchos::null);
  op_acc.AddAccumulationDelta(solution, phi, phi, dT, "cell");

  global_op->UpdateRHS(source, false);
  op.ApplyBCs(true, true, true);

  global_op->ComputeInverse();
  global_op->ApplyInverse(rhs, solution);

  int num_itrs = global_op->num_itrs();
  CHECK(num_itrs < 10);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual() 
              << " itr=" << num_itrs
              << " code=" << global_op->returned_code() << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("cell");
    GMV::open_data_file(*surfmesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::close_data_file();
  }
}


TEST(LAPLACE_BELTRAMI_CLOSED) {
  RunTest("diffusion operator");
}


TEST(LAPLACE_BELTRAMI_CLOSED_SFF) {
  RunTest("diffusion operator Sff");
}
