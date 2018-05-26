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
#include "LinearOperatorGMRES.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_AdvectionUpwind.hh"

#include "Analytic00.hh"


/* *****************************************************************
* Verify convergence of advection-diffusion solver for various BCs.
* **************************************************************** */
template<class Analytic>
void AdvectionDiffusion2D() 
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Advection-duffusion in 2D" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_advdiff.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an MSTK mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({Framework::MSTK}));
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 10, 10, gm);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // populate diffusion coefficient
  Analytic ana(mesh, 1.0, 1.0, 1);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, DOF_Type::SCALAR));

  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = xf[1] * xf[1];
    }
  }

  // create diffusion operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("diffusion operator");
  Teuchos::RCP<PDE_Diffusion> op_diff =
      Teuchos::rcp(new PDE_DiffusionMFD(olist, (Teuchos::RCP<const AmanziMesh::Mesh>) mesh));
  op_diff->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op_diff->global_operator()->DomainMap();

  // set up the diffusion operator
  op_diff->Setup(K, Teuchos::null, Teuchos::null);
  op_diff->UpdateMatrices(Teuchos::null, Teuchos::null);

  // get the global operator
  Teuchos::RCP<Operator> global_op = op_diff->global_operator();

  // create an advection operator  
  olist = plist.sublist("PK operator").sublist("advection operator");
  Teuchos::RCP<PDE_AdvectionUpwind> op_adv = Teuchos::rcp(new PDE_AdvectionUpwind(olist, global_op));

  // get a flux field
  Teuchos::RCP<CompositeVector> u = Teuchos::rcp(new CompositeVector(cvs));
  Epetra_MultiVector& uf = *u->ViewComponent("face");
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  Point vel(4.0, 4.0, 0.0);
  for (int f = 0; f < nfaces; f++) {
    uf[0][f] = vel * mesh->face_normal(f);
  }

  op_adv->Setup(*u);
  op_adv->UpdateMatrices(u.ptr());

  // Add an accumulation term.
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);  // solution at time T=0

  CompositeVector phi(cvs);
  phi.PutScalar(0.2);

  double dT = 0.02;
  Teuchos::RCP<PDE_Accumulation> op_acc = Teuchos::rcp(new PDE_Accumulation(AmanziMesh::CELL, global_op));
  op_acc->AddAccumulationDelta(solution, phi, phi, dT, "cell");

  // BCs and assemble
  op_diff->ApplyBCs(true, true);
  op_adv->ApplyBCs(bc, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // Create a preconditioner.
  ParameterList slist = plist.sublist("preconditioners");
  global_op->InitPreconditioner("Hypre AMG", slist);

  // Solve the problem.
  ParameterList lop_list = plist.sublist("solvers")
                                .sublist("AztecOO CG").sublist("gmres parameters");
  AmanziSolvers::LinearOperatorGMRES<Operator, CompositeVector, CompositeVectorSpace>
     solver(global_op, global_op);
  solver.Init(lop_list);

  CompositeVector& rhs = *global_op->rhs();
  int ierr = solver.ApplyInverse(rhs, solution);

  int num_itrs = solver.num_itrs();
  CHECK(num_itrs > 5 && num_itrs < 15);

  if (MyPID == 0) {
    std::cout << "pressure solver (gmres): ||r||=" << solver.residual() 
              << " itr=" << solver.num_itrs()
              << " code=" << solver.returned_code() << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("cell");
    GMV::open_data_file(*mesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::close_data_file();
  }
}

TEST(ADVECTION_DIFFUSION_2D) {
  AdvectionDiffusion2D<Analytic00>();
}
