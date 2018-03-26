/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  DG methods for linear advection equation c p + u . \nabla p = f.
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
#include "LinearOperatorGMRES.hh"
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"

// Operators
#include "AnalyticDG01.hh"
#include "AnalyticDG02.hh"

#include "OperatorAudit.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"


/* *****************************************************************
* This tests exactness/accuracy of the advection scheme.
* **************************************************************** */
TEST(OPERATOR_ADVECTION_DG) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D advection problem, discontinuous Galerkin" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_advection_dg.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK,STKMESH}));
  // RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 10, 10, gm);
  RCP<const Mesh> mesh = meshfactory("test/median7x8.exo", gm, true, true);

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // create global operator 
  // -- flux term
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator")
                               .sublist("flux operator");
  auto op_flux = Teuchos::rcp(new PDE_AdvectionRiemann(op_list, mesh));
  auto global_op = op_flux->global_operator();

  // -- volumetric advection term
  op_list = plist.get<Teuchos::ParameterList>("PK operator")
                 .sublist("advection operator");
  auto op_adv = Teuchos::rcp(new PDE_Abstract(op_list, global_op));

  // -- reaction term
  op_list = plist.get<Teuchos::ParameterList>("PK operator")
                 .sublist("reaction operator");
  auto op_reac = Teuchos::rcp(new PDE_Reaction(op_list, global_op));
  double Kreac = op_list.get<double>("coef");

  // populate problem data
  AnalyticDG02 ana(mesh, 2);

  // -- reaction coefficient
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);

  auto K = Teuchos::rcp(new CompositeVector(*cvs));
  auto Kc = K->ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost; c++) {
    const Point& xc = mesh->cell_centroid(c);
    (*Kc)[0][c] = Kreac;
  }

  // -- velocity coefficient
  int d = mesh->space_dimension();
  auto velc = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost));
  auto velf = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost));

  WhetStone::VectorPolynomial v(d, 2);
  v[0].Reshape(d, 2, true);
  v[0](1, 0) = 1.0;
  v[0](2, 0) = -1.0;

  v[1].Reshape(d, 2, true);
  v[1](1, 1) = 1.0;
  v[1](2, 2) = -1.0;
  
  for (int c = 0; c < ncells_wghost; ++c) {
    (*velc)[c] = v;
    (*velc)[c] *= -1.0;
  }

  for (int f = 0; f < nfaces_wghost; ++f) {
    (*velf)[f] = v * mesh->face_normal(f);
  }

  // -- source term
  WhetStone::Polynomial divv(d, 2);
  divv(0, 0) = 2.0;
  divv(1, 0) = -2.0;
  divv(1, 1) = -2.0;

  Point tmp(d);
  WhetStone::Polynomial sol, src;
  WhetStone::VectorPolynomial grad(d, 0);

  ana.TaylorCoefficients(tmp, 0.0, sol);
  grad.Gradient(sol); 
  src = Kreac * sol + v * grad + divv * sol;

  int order = op_list.get<int>("method order");

  WhetStone::Polynomial pc(2, order);
  WhetStone::NumericalIntegration numi(mesh);

  Epetra_MultiVector& rhs_c = *global_op->rhs()->ViewComponent("cell");
  for (int c = 0; c < ncells; ++c) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    WhetStone::Polynomial src_copy(src);
    src_copy.ChangeOrigin(xc);

    for (auto it = pc.begin(); it.end() <= pc.end(); ++it) {
      int n = it.PolynomialPosition();
      int k = it.MonomialOrder();

      double factor = numi.MonomialNaturalScale(k, volume);
      WhetStone::Polynomial cmono(2, it.multi_index(), factor);
      cmono.set_origin(xc);      

      WhetStone::Polynomial tmp = src_copy * cmono;      

      rhs_c[n][c] = numi.IntegratePolynomialCell(c, tmp);
    }
  }

  // populate the global operator
  op_flux->Setup(velc, velf);
  op_flux->UpdateMatrices(velf.ptr());

  op_adv->SetupPolyVector(velc);
  op_adv->UpdateMatrices();

  op_reac->Setup(Kc);
  op_reac->UpdateMatrices(Teuchos::null);

  // create preconditoner
  ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");

  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();
  global_op->InitPreconditioner("Hypre AMG", slist);

  // solve the problem
  ParameterList lop_list = plist.sublist("solvers")
                                .sublist("GMRES").sublist("gmres parameters");
  AmanziSolvers::LinearOperatorGMRES<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  solver.Init(lop_list);

  CompositeVector& rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver.ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "dG solver (gmres): ||r||=" << solver.residual() 
              << " itr=" << solver.num_itrs()
              << " code=" << solver.returned_code() 
              << " order=" << order << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("cell");
    GMV::open_data_file(*mesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::write_cell_data(p, 1, "gradx");
    GMV::write_cell_data(p, 1, "grady");
    GMV::close_data_file();
  }

  CHECK(solver.num_itrs() < 2000);

  // compute solution error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  if (MyPID == 0) {
    std::cout << "\nEXACT solution: " << sol << std::endl;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_err, pinf_err, solver.num_itrs());
    CHECK(pl2_err < 1e-10);
  }
}

