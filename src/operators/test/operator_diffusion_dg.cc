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
#include "NumericalIntegration.hh"
#include "Tensor.hh"

// Operators
#include "AnalyticDG01.hh"
#include "AnalyticDG02.hh"

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
  // RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 10, 10, gm);
  RCP<const Mesh> mesh = meshfactory("test/median7x8_filtered.exo", gm);

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

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

  // create boundary data
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator")
                               .sublist("diffusion operator dg");
  int order = op_list.get<int>("method order");
  int nk = (order + 1) * (order + 2) / 2;

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double> >& bc_value = bc->bc_value_vector(nk);

  for (int f = 0; f < nfaces_wghost; ++f) {
    const Point& xf = mesh->face_centroid(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;

      WhetStone::Polynomial coefs;
      ana.TaylorCoefficients(xf, 0.0, coefs);

      WhetStone::DenseVector data;
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

  WhetStone::Polynomial coefs, pc(2, order);
  WhetStone::NumericalIntegration numi(mesh);

  for (int c = 0; c < ncells; ++c) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    ana.SourceTaylor(xc, 0.0, coefs);
    coefs.set_origin(xc);

    for (auto it = pc.begin(); it.end() <= pc.end(); ++it) {
      int n = it.PolynomialPosition();
      int k = it.MonomialOrder();

      double factor = numi.MonomialNaturalScale(k, volume);
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

  CompositeVector& rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver.ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << solver.residual() 
              << " itr=" << solver.num_itrs()
              << " code=" << solver.returned_code() << std::endl;

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

  // compute pressure error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  if (MyPID == 0) {
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_err, pinf_err, solver.num_itrs());

    CHECK(pl2_err < 3e-2);
  }
}

