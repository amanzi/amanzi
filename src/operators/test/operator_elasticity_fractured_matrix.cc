/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

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
#include "Tensor.hh"
#include "WhetStoneMeshUtils.hh"

// Amanzi::Operators
#include "PDE_ElasticityFracturedMatrix.hh"

#include "AnalyticElasticity01.hh"
#include "Verification.hh"


/* *****************************************************************
* Elasticity model: exactness test.
***************************************************************** */
template<class Analytic>
double
RunTest(int icase, double mu, double lambda, double tol = 1e-10)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTEST: 3D elasticity: exactness test" << std::endl;

  // read parameter list
  // -- it specifies details of the mesh, elasticity operator, and solver
  std::string xmlFileName = "test/operator_elasticity_fractured_matrix.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("elasticity operator");

  if (icase == 3) {
    plist.sublist("regions")
      .sublist("fracture")
      .sublist("region: plane")
      .set<Teuchos::Array<double>>("normal", { 1.0, 1.0, 1.0 });
  }

  Teuchos::ParameterList region_list = plist.sublist("regions");
  auto gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  // create the MSTK mesh framework
  auto mesh_list = Teuchos::rcp(new Teuchos::ParameterList());
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);

  MeshFactory meshfactory(comm, gm, mesh_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh;
  if (icase == 1) mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 6, 6, 6);
  else mesh = meshfactory.create("test/tetrahedra.exo");

  // -- general information about mesh
  int ncells_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  int nnodes_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  int nfaces_fracture =
    mesh->getSetEntities("fracture", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL)
      .size();
  std::cout << "Mesh: cells=" << ncells_owned << " fracture faces=" << nfaces_fracture << std::endl;

  // select an analytic solution for error calculations and setup of
  // boundary conditions
  Analytic ana(mesh, mu, lambda);

  auto K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }

  // create a PDE: operator and boundary conditions
  // -- XML list speficies discretization method and location of degrees of freedom
  // -- (called schema). This seems redundant but only when use a low-order method.
  auto op = Teuchos::rcp(new PDE_ElasticityFracturedMatrix(op_list, mesh));
  op->Init(op_list);

  // populate boundary conditions: type (called model) and value
  // -- full velocity at boundary nodes (a vector)
  auto bcv = Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT));
  std::vector<int>& bcv_model = bcv->bc_model();
  std::vector<Point>& bcv_value = bcv->bc_value_point();

  for (int v = 0; v < nnodes_wghost; ++v) {
    auto xv = mesh->getNodeCoordinate(v);

    if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 || fabs(xv[1]) < 1e-6 ||
        fabs(xv[1] - 1.0) < 1e-6 || fabs(xv[2]) < 1e-6 || fabs(xv[2] - 1.0) < 1e-6) {
      bcv_model[v] = OPERATOR_BC_DIRICHLET;
      bcv_value[v] = ana.velocity_exact(xv, 0.0);
    }
  }
  op->AddBCs(bcv, bcv);

  // -- normal component of velocity on boundary faces (a scalar)
  auto bcf =
    Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bcf_model = bcf->bc_model();
  std::vector<double>& bcf_value = bcf->bc_value();

  int dir;
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->getFaceCentroid(f);
    const Point& normal = getFaceNormalExterior(*mesh, f, &dir);
    double area = mesh->getFaceArea(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
        fabs(xf[1] - 1.0) < 1e-6 || fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) {
      bcf_model[f] = OPERATOR_BC_DIRICHLET;
      bcf_value[f] = (ana.velocity_exact(xf, 0.0) * normal) / area;
    }
  }
  for (int i = 0; i < 3; ++i) op->AddBCs(bcf, bcf); // FIXME logic of BCs should be simplified

  // create and initialize solution
  auto global_op = op->global_operator();

  const CompositeVectorSpace& cvs = global_op->DomainMap();
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);

  // create source
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("node", true);
  src.PutScalar(0.0);

  for (int c = 0; c < ncells_owned; ++c) {
    double vol = mesh->getCellVolume(c);

    auto nodes = mesh->getCellNodes(c);
    int nnodes = nodes.size();

    for (int n = 0; n < nnodes; ++n) {
      int v = nodes[n];
      auto xv = mesh->getNodeCoordinate(v);
      Point tmp(ana.source_exact(xv, 0.0));
      for (int k = 0; k < 3; ++k) src[k][v] += tmp[k] * (vol / nnodes);
    }
  }
  source.GatherGhostedToMaster("node");

  // populate the elasticity operator
  op->SetTensorCoefficient(K);
  op->UpdateMatrices();

  VerificationCV ver(global_op);
  ver.CheckMatrixSPD(true, true);

  // get and assemble the global operator (volume was included)
  global_op->UpdateRHS(source, true);
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->set_inverse_parameters(
    "Hypre AMG", plist.sublist("preconditioners"), "PCG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  ver.CheckPreconditionerSPD(1e-12, true, true);

  CompositeVector& rhs = *global_op->rhs();
  global_op->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "elasticity solver (PCG): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code()
              << std::endl;
  }

  // verify shear strain
  if (icase < 3) {
    for (int c = 0; c < ncells_owned; ++c) {
      auto Tc = op->ComputeCellStrain(solution, c);
      CHECK_CLOSE(Tc.Trace(), -1.0, 1e-12);
    }
  }

  // compute displacement error
  double unorm, ul2_err, uinf_err;
  ana.VectorNodeError(solution, 0.0, unorm, ul2_err, uinf_err);
  ul2_err /= unorm;

  if (MyPID == 0) {
    printf("L2(u)=%12.8g  Inf(u)=%12.8g  itr=%3d\n", ul2_err, uinf_err, global_op->num_itrs());

    CHECK(ul2_err < tol);
    CHECK(global_op->num_itrs() < 15);
  }

  return ul2_err;
}


TEST(OPERATOR_ELASTICITY_FRACTURED_MATRIX_EXACTNESS)
{
  RunTest<AnalyticElasticity01>(1, 1.0, 0.0);
}

TEST(OPERATOR_ELASTICITY_FRACTURED_MATRIX_TETRAHEDRA)
{
  RunTest<AnalyticElasticity01>(2, 1.0, 0.0);
}

TEST(OPERATOR_ELASTICITY_FRACTURED_MATRIX_SYMMETRY)
{
  RunTest<AnalyticElasticity01>(3, 1.0, 0.0, 0.1);
}
