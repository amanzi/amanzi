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
#include "PDE_Elasticity.hh"

#include "AnalyticElasticity01.hh"
#include "AnalyticElasticity03.hh"
#include "Verification.hh"


/* *****************************************************************
* Elasticity model: exactness test.
***************************************************************** */
template<class Analytic>
double
RunTest(int icase,
        const std::string& solver,
        double mu,
        double lambda,
        bool flag,
        double tol = 1e-10,
        int scale = 1)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0)
    std::cout << "\n======================================\n"
              << "TEST: 2D elasticity: exactness test #" << icase
              << "\n======================================" << std::endl;

  // read parameter list
  // -- it specifies details of the mesh, elasticity operator, and solver
  std::string xmlFileName = "test/operator_elasticity.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("elasticity operator");

  if (icase == 2 || icase == 3 || icase == 4) {
    op_list.sublist("schema").set<std::string>("method", "elasticity");
  }

  // create the MSTK mesh framework
  // -- geometric model is not created. Instead, we specify boundary conditions
  // -- using centroids of mesh faces.
  double Lx(1.0), Ly(2.0);
  auto mesh_list = Teuchos::rcp(new Teuchos::ParameterList());
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);
  MeshFactory meshfactory(comm, Teuchos::null, mesh_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, Lx, Ly, 4 * scale, 5 * scale);

  // -- general information about mesh
  int ncells_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  int nnodes_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  // select an analytic solution for error calculations and setup of
  // boundary conditions
  Analytic ana(mesh, mu, lambda, flag);

  auto K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }

  // create a PDE: operator and boundary conditions
  // -- XML list speficies discretization method and location of degrees of freedom
  // -- (called schema). This seems redundant but only when use a low-order method.
  Teuchos::RCP<PDE_Elasticity> op = Teuchos::rcp(new PDE_Elasticity(op_list, mesh));
  op->Init(op_list);

  // populate boundary conditions: type (called model) and value
  // -- normal component of velocity on boundary faces (a scalar)
  int ndir(0), nshear(0), nkinematic(0);
  if (icase == 1) {
    Teuchos::RCP<BCs> bcf =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    std::vector<int>& bcf_model = bcf->bc_model();
    std::vector<double>& bcf_value = bcf->bc_value();

    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& xf = mesh->getFaceCentroid(f);
      const Point& normal = mesh->getFaceNormal(f);
      double area = mesh->getFaceArea(f);

      if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - Lx) < 1e-6 || fabs(xf[1]) < 1e-6 ||
          fabs(xf[1] - Ly) < 1e-6) {
        bcf_model[f] = OPERATOR_BC_DIRICHLET;
        bcf_value[f] = (ana.velocity_exact(xf, 0.0) * normal) / area;
        ndir++;
      }
    }
    op->AddBCs(bcf, bcf);
  }

  // -- full velocity at boundary nodes (a vector)
  if (icase == 1 || icase == 2) {
    Teuchos::RCP<BCs> bcv =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT));
    std::vector<int>& bcv_model = bcv->bc_model();
    std::vector<Point>& bcv_value = bcv->bc_value_point();

    for (int v = 0; v < nnodes_wghost; ++v) {
      auto xv = mesh->getNodeCoordinate(v);

      if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - Lx) < 1e-6 || fabs(xv[1]) < 1e-6 ||
          fabs(xv[1] - Ly) < 1e-6) {
        bcv_model[v] = OPERATOR_BC_DIRICHLET;
        bcv_value[v] = ana.velocity_exact(xv, 0.0);
        ndir++;
      }
    }
    op->AddBCs(bcv, bcv);
  }

  // -- full velocity at boundary nodes (a vector)
  if (icase == 3) {
    double eps(1e-4);
    auto bcf =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    std::vector<int>& bcf_model = bcf->bc_model();
    std::vector<double>& bcf_value = bcf->bc_value();

    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& xf = mesh->getFaceCentroid(f);
      const Point& normal = mesh->getFaceNormal(f);
      const Point& tau = mesh->getEdgeVector(f);
      double area = mesh->getFaceArea(f);

      if ((fabs(xf[1]) < 1e-6 && xf[0] > eps && xf[0] < Lx - eps) ||
          (fabs(xf[1] - Ly) < 1e-6 && xf[0] > eps && xf[0] < Lx - eps)) {
        bcf_model[f] = OPERATOR_BC_SHEAR_STRESS;
        bcf_value[f] = ((ana.stress_exact(xf, 0.0) * tau) * normal) / area / area;
        nshear++;
      }
    }
    op->AddBCs(bcf, bcf);

    auto bcv =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::SCALAR));
    std::vector<int>& bcv_model = bcv->bc_model();
    std::vector<double>& bcv_value = bcv->bc_value();

    for (int v = 0; v < nnodes_wghost; ++v) {
      auto xv = mesh->getNodeCoordinate(v);

      if ((fabs(xv[1]) < 1e-6 && xv[0] > eps && xv[0] < Lx - eps) ||
          (fabs(xv[1] - Ly) < 1e-6 && xv[0] > eps && xv[0] < Lx - eps)) {
        auto normal = WhetStone::getNodeUnitNormal(*mesh, v);
        bcv_model[v] = OPERATOR_BC_KINEMATIC;
        bcv_value[v] = ana.velocity_exact(xv, 0.0) * normal;
        nkinematic++;
      }
    }
    op->AddBCs(bcv, bcv);

    auto bcp =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT));
    std::vector<int>& bcp_model = bcp->bc_model();
    std::vector<Point>& bcp_value = bcp->bc_value_point();

    for (int v = 0; v < nnodes_wghost; ++v) {
      auto xv = mesh->getNodeCoordinate(v);
      if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - Lx) < 1e-6) {
        bcp_model[v] = OPERATOR_BC_DIRICHLET;
        bcp_value[v] = ana.velocity_exact(xv, 0.0);
        ndir++;
      }
    }
    op->AddBCs(bcp, bcp);
  }

  if (icase == 4) {
    auto bcv =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT));
    std::vector<int>& bcv_model = bcv->bc_model();
    std::vector<AmanziGeometry::Point>& bcv_value = bcv->bc_value_point();

    for (int v = 0; v < nnodes_wghost; ++v) {
      auto xv = mesh->getNodeCoordinate(v);

      if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - Lx) < 1e-6 || fabs(xv[1]) < 1e-6) {
        bcv_model[v] = OPERATOR_BC_DIRICHLET;
        bcv_value[v] = ana.velocity_exact(xv, 0.0);
        ndir++;
      }
    }
    op->AddBCs(bcv, bcv);

    auto bcf =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::POINT));
    std::vector<int>& bcf_model = bcf->bc_model();
    std::vector<AmanziGeometry::Point>& bcf_value = bcf->bc_value_point();

    for (int f = 0; f < nfaces_wghost; ++f) {
      auto xf = mesh->getFaceCentroid(f);

      if (fabs(xf[1] - Ly) < 1e-6) {
        double area = mesh->getFaceArea(f);
        auto normal = mesh->getFaceNormal(f);
        bcf_model[f] = OPERATOR_BC_NORMAL_STRESS;
        bcf_value[f] = (ana.stress_exact(xf, 0.0) * normal) / area;
        nshear++;
      }
    }
    op->AddBCs(bcf, bcf);
  }

  // create and initialize solution
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();
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
      for (int k = 0; k < 2; ++k) src[k][v] += tmp[k] * (vol / nnodes);
    }
  }
  source.GatherGhostedToMaster("node");

  // populate the elasticity operator
  op->SetTensorCoefficient(K);
  op->UpdateMatrices();

  // get and assemble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, true); // FIXME volume is missing but RHS is zero
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->set_inverse_parameters(
    "Hypre AMG", plist.sublist("preconditioners"), solver, plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  if (icase == 1 || icase == 2) {
    ver.CheckMatrixSPD(true, true);
    ver.CheckPreconditionerSPD(1e-12, true, true);
  }

  CompositeVector& rhs = *global_op->rhs();
  global_op->ApplyInverse(rhs, solution);

  if (icase == 1 || icase == 2) {
    ver.CheckResidual(solution, 1.0e-13);
  }

  if (MyPID == 0) {
    std::cout << ana.Tensor(mesh->getCellCentroid(0), 0.0) << std::endl;
    std::cout << "elasticity solver (" << solver << "): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code()
              << std::endl
              << "BCs: noslip: " << ndir << ", kinematic: " << nkinematic
              << ", shear stress: " << nshear << std::endl;
  }

  // compute velocity error
  double unorm, ul2_err, uinf_err;
  ana.VectorNodeError(solution, 0.0, unorm, ul2_err, uinf_err);
  ul2_err /= unorm;

  if (MyPID == 0) {
    printf("L2(u)=%12.8g  Inf(u)=%12.8g  itr=%3d\n", ul2_err, uinf_err, global_op->num_itrs());

    CHECK(ul2_err < tol);
    CHECK(global_op->num_itrs() < 15);
  }

  // compute volumetric strain error
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh);
  cvs2.SetGhosted(true);
  cvs2.SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  CompositeVector e(cvs2);
  op->ComputeVolumetricStrain(solution, e);
  const auto& e_c = *e.ViewComponent("cell");

  double el2_err(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    const auto& xc = mesh->getCellCentroid(c);
    double e_exact = ana.volumetric_strain_exact(xc, 0.0);
    el2_err += std::pow(e_exact - e_c[0][c], 2);
  }
  ana.GlobalOp("sum", &el2_err, 1);
  el2_err = std::sqrt(el2_err / ncells_owned);
  CHECK(el2_err < tol);

  return ul2_err;
}


TEST(OPERATOR_ELASTICITY_EXACTNESS_BERNARDI_RAUGEL)
{
  RunTest<AnalyticElasticity01>(1, "PCG", 1.0, 0.0, false);
}

TEST(OPERATOR_ELASTICITY_EXACTNESS_ELASTICITY_SCALAR_COEFFICIENT)
{
  RunTest<AnalyticElasticity01>(2, "PCG", 2.0, 0.0, false);
}

TEST(OPERATOR_ELASTICITY_EXACTNESS_ELASTICITY_SCALAR_TENSOR)
{
  RunTest<AnalyticElasticity01>(2, "PCG", 2.0, 0.0, true);
}

TEST(OPERATOR_ELASTICITY_EXACTNESS_ELASTICITY_FULL_TENSOR)
{
  RunTest<AnalyticElasticity01>(2, "PCG", 2.0, 1.0, true);
}

TEST(OPERATOR_ELASTICITY_EXACTNESS_KINEMATIC_SCALAR_COEFFICIENT)
{
  RunTest<AnalyticElasticity01>(3, "GMRES", 2.0, 0.0, false);
}

TEST(OPERATOR_ELASTICITY_EXACTNESS_KINEMATIC_SCALAR_TENSOR)
{
  RunTest<AnalyticElasticity01>(3, "GMRES", 2.0, 0.0, true);
}

TEST(OPERATOR_ELASTICITY_EXACTNESS_KINEMATIC_FULL_TENSOR)
{
  RunTest<AnalyticElasticity01>(3, "GMRES", 2.0, 1.0, true);
}

TEST(OPERATOR_ELASTICITY_EXACTNESS_NORMAL_STRESS_FULL_TENSOR)
{
  RunTest<AnalyticElasticity01>(4, "GMRES", 2.0, 1.0, true);
}

TEST(OPERATOR_CONVERGENCE_ELASTICITY_FULL_TENSOR_BERNARDI_RAUGEL)
{
  double E(1e+10), nu(0.2);
  double mu = E / (2 * (1 + nu));
  // double lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // 3D
  double lambda = E * nu / ((1 + nu) * (1 - nu)); // 2D
  double err1 = RunTest<AnalyticElasticity03>(1, "PCG", mu, lambda, true, 5e-2, 5);
  double err2 = RunTest<AnalyticElasticity03>(1, "PCG", mu, lambda, true, 5e-2, 10);

  CHECK(err1 / err2 > 3.8);
}

TEST(OPERATOR_CONVERGENCE_ELASTICITY_FULL_TENSOR_NODAL)
{
  double E(1e+10), nu(0.2);
  double mu = E / (2 * (1 + nu));
  double lambda = E * nu / ((1 + nu) * (1 - nu)); // 2D
  double err1 = RunTest<AnalyticElasticity03>(2, "PCG", mu, lambda, true, 5e-2, 5);
  double err2 = RunTest<AnalyticElasticity03>(2, "PCG", mu, lambda, true, 5e-2, 10);

  CHECK(err1 / err2 > 3.8);
}
