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
#include "LinearOperatorFactory.hh"
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

// Operators
#include "AnalyticDG01.hh"
#include "AnalyticDG02.hh"

#include "OperatorAudit.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionDG.hh"


/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* **************************************************************** */
void OperatorDiffusionDG(std::string solver_name) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D elliptic problem, dG method, solver: " << solver_name << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK,STKMESH}));
  // RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 2, gm);
  RCP<const Mesh> mesh = meshfactory("test/median7x8_filtered.exo", gm);
  // RCP<const Mesh> mesh = meshfactory("test/triangular8_clockwise.exo", gm);

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // modify diffusion coefficient
  auto Kc = std::make_shared<std::vector<WhetStone::Tensor> >();
  auto Kf = std::make_shared<std::vector<double> >();

  AnalyticDG02 ana(mesh, 2);

  for (int c = 0; c < ncells_wghost; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Ktmp = ana.Tensor(xc, 0.0);
    Kc->push_back(Ktmp);
  }

  for (int f = 0; f < nfaces_wghost; f++) {
    double area = mesh->face_area(f);
    Kf->push_back(40.0 / area);
  }

  // create boundary data. We use full Taylor expansion of boundary data in
  // the vicinity of domain boundary.
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator")
                               .sublist("diffusion operator dg");
  int order = op_list.get<int>("method order");
  int nk = (order + 1) * (order + 2) / 2;

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double> >& bc_value = bc->bc_value_vector(nk);

  WhetStone::Polynomial coefs;
  WhetStone::DenseVector data;

  for (int f = 0; f < nfaces_wghost; ++f) {
    const Point& xf = mesh->face_centroid(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;

      ana.SolutionTaylor(xf, 0.0, coefs);
      coefs.GetPolynomialCoefficients(data);

      for (int i = 0; i < data.NumRows(); ++i) {
        bc_value[f][i] = data(i);
      }
    } else if (fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_NEUMANN;

      ana.SolutionTaylor(xf, 0.0, coefs);
      coefs.GetPolynomialCoefficients(data);

      for (int i = 0; i < data.NumRows(); ++i) {
        bc_value[f][i] = data(i);
      }
    }
  }

  // create diffusion operator 
  // -- primary term
  Teuchos::RCP<PDE_DiffusionDG> op = Teuchos::rcp(new PDE_DiffusionDG(op_list, mesh));
  auto global_op = op->global_operator();

  // -- boundary conditions
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create source and add it to the operator
  CompositeVector src(cvs);
  Epetra_MultiVector& src_c = *src.ViewComponent("cell");

  WhetStone::Polynomial pc(2, order);
  WhetStone::NumericalIntegration numi(mesh, false);

  for (int c = 0; c < ncells; ++c) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    ana.SourceTaylor(xc, 0.0, coefs);
    coefs.set_origin(xc);

    for (auto it = pc.begin(); it.end() <= pc.end(); ++it) {
      int n = it.PolynomialPosition();
      int k = it.MonomialSetOrder();

      double factor = numi.MonomialNaturalScales(c, k);
      WhetStone::Polynomial cmono(2, it.multi_index(), factor);
      cmono.set_origin(xc);      

      WhetStone::Polynomial tmp = coefs * cmono;      

      src_c[n][c] = numi.IntegratePolynomialCell(c, tmp);
    }
  }

  // populate the diffusion operator
  op->SetProblemCoefficients(Kc, Kf);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update the source term
  *global_op->rhs() = src;

  // apply BCs (primary=true, eliminate=true) and assemble
  op->ApplyBCs(true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // Test SPD properties of the matrix.
  Operators::CheckMatrixSymmetry(global_op->A());
  Operators::CheckMatrixCoercivity(global_op->A());

  // create preconditoner using the base operator class
  ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->InitializePreconditioner(slist);
  global_op->UpdatePreconditioner();

  // solve the problem
  ParameterList lop_list = plist.sublist("solvers");
  AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> solverfactory;
  auto solver = solverfactory.Create(solver_name, lop_list, global_op, global_op);

  CompositeVector& rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << solver->residual() 
              << " itr=" << solver->num_itrs()
              << " code=" << solver->returned_code() << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("cell");
    GMV::open_data_file(*mesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::write_cell_data(p, 1, "gradx");
    GMV::write_cell_data(p, 1, "grady");
    GMV::close_data_file();
  }

  CHECK(solver->num_itrs() < 200);

  // compute pressure error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean);

  if (MyPID == 0) {
    printf("Mean:  L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_mean, pinf_mean, solver->num_itrs());
    printf("Total: L2(p)=%9.6f  Inf(p)=%9.6f\n", pl2_err, pinf_err);

    CHECK(pl2_err < 3e-2);
  }
}

TEST(OPERATOR_DIFFUSION_DG) {
  OperatorDiffusionDG("AztecOO CG");
  OperatorDiffusionDG("Amesos1");
  OperatorDiffusionDG("Amesos2");
}
