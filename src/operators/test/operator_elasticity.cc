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
#include "LinearOperatorGMRES.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "PDE_Elasticity.hh"

#include "MeshDeformation.hh"
#include "AnalyticElasticity01.hh"
#include "AnalyticElasticity03.hh"
#include "Verification.hh"

/* *****************************************************************
* Elasticity model: exactness test.
***************************************************************** */
void RunTest(const std::string& operator_name, bool bc_on_faces) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elasticity: exactness test for " << operator_name << std::endl;

  // read parameter list
  // -- it specifies details of the mesh, elasticity operator, and solver
  std::string xmlFileName = "test/operator_elasticity.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create the MSTK mesh framework 
  // -- geometric model is not created. Instead, we specify boundary conditions
  // -- using centroids of mesh faces.
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 4, 5);
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist(operator_name);

  // -- general information about mesh
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  // select an analytic solution for error calculations and setup of
  // boundary conditions
  AnalyticElasticity01 ana(mesh);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }

  // populate boundary conditions: type (called model) and value
  // -- normal component of velocity on boundary faces (a scalar)
  Teuchos::RCP<BCs> bcf = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
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

  // -- full velocity at boundary nodes (a vector)
  Point xv(2);
  Teuchos::RCP<BCs> bcv = Teuchos::rcp(new BCs(mesh, AmanziMesh::NODE, WhetStone::DOF_Type::POINT));
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

  // create and initialize a PDE 
  Teuchos::RCP<PDE_Elasticity> op = Teuchos::rcp(new PDE_Elasticity(op_list, mesh));
  op->SetBCs(bcv, bcv);
  if (bc_on_faces) op->AddBCs(bcf, bcf);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create and initialize solution
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);

  // create source 
  CompositeVector source(cvs);
  if (cvs.HasComponent("node")) {
  Epetra_MultiVector& src = *source.ViewComponent("node");

  for (int v = 0; v < nnodes; v++) {
    mesh->node_get_coordinates(v, &xv);
    Point tmp(ana.source_exact(xv, 0.0));
    for (int k = 0; k < 2; ++k) src[k][v] = tmp[k];
  }
  }

  // populate the elasticity operator
  op->SetTensorCoefficient(K);
  op->UpdateMatrices();

  // get and assemble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, true);  // FIXME
  op->ApplyBCs(true, true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create preconditoner using the base operator class
  Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->InitializePreconditioner(slist);
  global_op->UpdatePreconditioner();

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD(true, true);
  ver.CheckPreconditionerSPD(1e-12, true, true);

  // solve the problem
  Teuchos::ParameterList lop_list = plist.sublist("solvers")
                                         .sublist("PCG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      pcg(global_op, global_op);
  pcg.Init(lop_list);

  CompositeVector& rhs = *global_op->rhs();
  pcg.ApplyInverse(rhs, solution);

  ver.CheckResidual(solution, 1.0e-12);

  if (MyPID == 0) {
    std::cout << "elasticity solver (pcg): ||r||=" << pcg.residual() 
              << " itr=" << pcg.num_itrs()
              << " code=" << pcg.returned_code() << std::endl;
  }

  // compute velocity error
  if (cvs.HasComponent("node")) {
  double unorm, ul2_err, uinf_err;
  ana.VectorNodeError(solution, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    ul2_err /= unorm;
      printf("L2(u)=%12.8f  Inf(u)=%12.8f  itr=%3d\n",
        ul2_err, uinf_err, pcg.num_itrs());

    CHECK(ul2_err < 0.1);
    CHECK(pcg.num_itrs() < 15);
  }
  }
}


TEST(OPERATOR_ELASTICITY_EXACTNESS) {
  RunTest("elasticity operator 1", true);
  RunTest("elasticity operator 2", false);
}


/* *****************************************************************
* Elasticity model: convergence test.
***************************************************************** */
void RunTestLocalStress(const std::string& filename) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elasticity: convergence test: local stress scheme" << std::endl;

  // read parameter list
  // -- it specifies details of the mesh, elasticity operator, and solver
  Teuchos::ParameterXMLFileReader xmlreader("test/operator_elasticity.xml");
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create the MSTK mesh framework 
  Teuchos::ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  // Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 2, 2);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(filename);
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("elasticity operator local stress");

  // DeformMesh(mesh, 1, 1.0);

  // -- general information about mesh
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // select an analytic solution for error calculations and setup of
  // boundary conditions
  AnalyticElasticity03 ana(mesh);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells_wghost; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }

  // populate boundary conditions: type (called model) and value
  // -- full velocity on boundary faces
  Teuchos::RCP<BCs> bcf = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::VECTOR));
  std::vector<int>& bcf_model = bcf->bc_model();
  std::vector<std::vector<double> >& bcf_value = bcf->bc_value_vector(4);

  const auto& fmap = mesh->face_map(true);
  const auto& bmap = mesh->exterior_face_map(true);

  AmanziMesh::Entity_ID_List nodes;
  Point xv(2), uf(2);

  for (int bf = 0; bf < bmap.NumMyElements(); ++bf) {
    int f = fmap.LID(bmap.GID(bf));
    const Point& xf = mesh->face_centroid(f);

    mesh->face_get_nodes(f, &nodes);

    bcf_model[f] = OPERATOR_BC_DIRICHLET;

    for (int i = 0; i < nodes.size(); ++i) {
      mesh->node_get_coordinates(nodes[i], &xv);
      uf = ana.velocity_exact((xf + xv) / 2, 0.0);
      bcf_value[f][2 * i] = uf[0];
      bcf_value[f][2 * i + 1] = uf[1];
    }
  }

  // create and initialize a PDE 
  Teuchos::RCP<PDE_Elasticity> op = Teuchos::rcp(new PDE_Elasticity(op_list, mesh));
  op->AddBCs(bcf, bcf);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create and initialize solution
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);

  // create source 
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) {
    double volume = mesh->cell_volume(c);
    const Point& xc = mesh->cell_centroid(c);
    auto tmp = ana.source_exact(xc, 0.0);
    for (int k = 0; k < 2; ++k) src[k][c] = tmp[k] * volume;
  }

  // populate the elasticity operator
  op->SetTensorCoefficient(K);
  op->UpdateMatrices();

  // get and assmeble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, true);
  op->ApplyBCs(true, true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create preconditoner using the base operator class
  Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->InitializePreconditioner(slist);
  global_op->UpdatePreconditioner();

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  // ver.CheckSpectralBounds(3);

  // solve the problem
  Teuchos::ParameterList lop_list = plist.sublist("solvers")
                                         .sublist("GMRES").sublist("gmres parameters");
  AmanziSolvers::LinearOperatorGMRES<Operator, CompositeVector, CompositeVectorSpace>
      gmres(global_op, global_op);
  gmres.Init(lop_list);

  CompositeVector& rhs = *global_op->rhs();
  int ierr = gmres.ApplyInverse(rhs, solution);

  ver.CheckResidual(solution, 1.0e-12);

  if (MyPID == 0) {
    std::cout << "elasticity solver (gmres): ||r||=" << gmres.residual() 
              << " itr=" << gmres.num_itrs()
              << " code=" << gmres.returned_code() << std::endl;
  }

  // compute velocity error
  double unorm, ul2_err, uinf_err;
  ana.VectorCellError(solution, 0.0, unorm, ul2_err, uinf_err);
  int ndofs = 2 * mesh->cell_map(false).NumGlobalElements();

  if (MyPID == 0) {
    ul2_err /= unorm;
    printf("L2(u)=%12.10f  Inf(u)=%12.10f  itr=%3d  DOFs=%6d\n", ul2_err, uinf_err, gmres.num_itrs(), ndofs);

    CHECK(ul2_err < 0.1);
    // CHECK(gmres.num_itrs() < 15);
  }
}

TEST(OPERATOR_ELASTICITY_LOCAL_STRESS) {
  // RunTestLocalStress("test/triangular8.exo");
  RunTestLocalStress("test/triangular8.exo");
  RunTestLocalStress("test/triangular16.exo");
  RunTestLocalStress("test/triangular32.exo");
  RunTestLocalStress("test/triangular64.exo");
  RunTestLocalStress("test/triangular128.exo");

  RunTestLocalStress("test/mesh_poly20x20.exo");
  RunTestLocalStress("test/mesh_poly40x40.exo");
  RunTestLocalStress("test/mesh_poly80x80.exo");
  RunTestLocalStress("test/mesh_poly160x160.exo");
  RunTestLocalStress("test/mesh_poly320x320.exo");

  RunTestLocalStress("test/median15x16.exo");
  RunTestLocalStress("test/median32x33.exo");
  RunTestLocalStress("test/median63x64.exo");
  RunTestLocalStress("test/median127x128.exo");
  RunTestLocalStress("test/median255x256.exo");

  RunTestLocalStress("test/circle_quad10.exo");
  RunTestLocalStress("test/circle_quad20.exo");
  RunTestLocalStress("test/circle_quad40.exo");
  RunTestLocalStress("test/circle_quad80.exo");
  RunTestLocalStress("test/circle_quad160.exo");
}
