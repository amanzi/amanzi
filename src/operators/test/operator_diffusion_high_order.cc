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
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "Tensor.hh"

// Operators
#include "Analytic00.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"


/* *****************************************************************
* Exactness test for high-order Lagrange elements.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_HIGH_ORDER) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, high-order Lagrange" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK, STKMESH}));
  // RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 4, 4, gm);
  RCP<const Mesh> mesh = meshfactory("test/median15x16.exo", gm);

  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  // create boundary data (no mixed bc)
  Analytic00 ana(mesh, 1.0, 1.0);

  Point xv(2);
  Teuchos::RCP<BCs> bc_v = Teuchos::rcp(new BCs(mesh, AmanziMesh::NODE, DOF_Type::SCALAR));
  std::vector<int>& bc_model_v = bc_v->bc_model();
  std::vector<double>& bc_value_v = bc_v->bc_value();

  for (int v = 0; v < nnodes_wghost; v++) {
    mesh->node_get_coordinates(v, &xv);
    if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
        fabs(xv[1]) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bc_model_v[v] = OPERATOR_BC_DIRICHLET;
      bc_value_v[v] = ana.pressure_exact(xv, 0.0);
    }
  }

  Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, DOF_Type::SCALAR));
  std::vector<int>& bc_model_f = bc_f->bc_model();
  std::vector<double>& bc_value_f = bc_f->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model_f[f] = OPERATOR_BC_DIRICHLET;
      bc_value_f[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator")
                               .sublist("diffusion operator high-order");
  Teuchos::RCP<PDE_Abstract> op = Teuchos::rcp(new PDE_Abstract(op_list, mesh));
  op->AddBCs(bc_v, bc_v);
  op->AddBCs(bc_f, bc_f);
  
  // populate the diffusion operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->Init();
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs(true, true);

  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create preconditoner using the base operator class
  ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
  global_op->InitPreconditioner("Hypre AMG", slist);

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
  }

  // compute pressure error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("node", false);

  double pnorm, pl2_err, pinf_err, hnorm, ph1_err;
  ana.ComputeNodeError(p, 0.0, pnorm, pl2_err, pinf_err, hnorm, ph1_err);

  if (MyPID == 0) {
    pl2_err /= pnorm;
    ph1_err /= hnorm;
    printf("L2(p)=%9.6f  H1(p)=%9.6f  itr=%3d\n", pl2_err, ph1_err, solver.num_itrs());

    CHECK(pl2_err < 1e-5 && ph1_err < 2e-5);
  }
}

