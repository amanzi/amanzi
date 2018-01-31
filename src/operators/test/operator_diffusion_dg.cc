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
#include "GMVMesh.hh"
#include "CompositeVector.hh"
#include "LinearOperatorPCG.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"

// Operators
#include "Analytic00.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionDG.hh"


/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_DG) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, discontinuous Galerkin" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK,STKMESH}));
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 10, 10, gm);
  // RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  Analytic00 ana(mesh, 1.0, 1.0, 1);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator dg");
  Teuchos::RCP<PDE_DiffusionDG> op = Teuchos::rcp(new PDE_DiffusionDG(op_list, mesh));
  // op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create source and add it to the operator

  // populate the diffusion operator
  op->SetTensorCoefficient(K);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update the source term
  Teuchos::RCP<Operator> global_op = op->global_operator();
  // global_op->UpdateRHS(source, true);

  // apply BCs (primary=true, eliminate=true) and assemble
  op->ApplyBCs(true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create preconditoner using the base operator class
  ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
  global_op->InitPreconditioner("Hypre AMG", slist);

  // Test SPD properties of the preconditioner.
  CompositeVector a(cvs), ha(cvs), b(cvs), hb(cvs);
  a.Random();
  b.Random();
  global_op->ApplyInverse(a, ha);
  global_op->ApplyInverse(b, hb);

  double ahb, bha, aha, bhb;
  a.Dot(hb, &ahb);
  b.Dot(ha, &bha);
  a.Dot(ha, &aha);
  b.Dot(hb, &bhb);

  if (MyPID == 0) {
    std::cout << "Preconditioner:\n"
              << "  Symmetry test: " << ahb << " = " << bha << std::endl;
    std::cout << "  Positivity test: " << aha << " " << bhb << std::endl;
  }
  CHECK_CLOSE(ahb, bha, 1e-12 * fabs(ahb));
  CHECK(aha > 0.0);
  CHECK(bhb > 0.0);

  // solve the problem
  ParameterList lop_list = plist.sublist("solvers")
                                .sublist("AztecOO CG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  solver.Init(lop_list);

  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver.ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << solver.residual() 
              << " itr=" << solver.num_itrs()
              << " code=" << solver.returned_code() << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("node");
    GMV::open_data_file(*mesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_node_data(p, 0, "solution");
    GMV::close_data_file();
  }

  CHECK(solver.num_itrs() < 10);

  // compute pressure error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  if (MyPID == 0) {
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_err, pinf_err, solver.num_itrs());

    CHECK(pl2_err < 3e-3);
    CHECK(solver.num_itrs() < 10);
  }
}

