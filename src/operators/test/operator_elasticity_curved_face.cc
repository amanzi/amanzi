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
#include "PDE_ElasticityCurvedFace.hh"

#include "AnalyticElasticity01.hh"

/* *****************************************************************
* Elasticity model: exactness test.
***************************************************************** */
void
RunTest(double mu, double lambda, const std::string& filename = "")
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTEST: Elasticity with curved faces \"" << filename << "\"\n";

  // read parameter list
  // -- it specifies details of the mesh, elasticity operator, and solver
  std::string xmlFileName = "test/operator_elasticity_curved_face.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create the MSTK mesh framework
  // -- geometric model is not created. Instead, we specify boundary conditions
  // -- using centroids of mesh faces.
  auto mesh_list = Teuchos::rcp(new Teuchos::ParameterList());
  mesh_list->set<bool>("request faces", true);
  MeshFactory meshfactory(comm, Teuchos::null, mesh_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  int d;
  Teuchos::RCP<Mesh> mesh;
  if (filename == "") {
    mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 5, 5);
    d = 2;
  } else {
    mesh = meshfactory.create(filename);
    d = 3;
  }

  // -- general information about mesh
  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int ncells_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  int nfaces_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  // select an analytic solution for error calculations and setup of
  // boundary conditions
  AnalyticElasticity01 ana(mesh, mu, lambda);

  auto K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K->push_back(Kc);
  }

  // create a PDE: operator and boundary conditions
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("elasticity operator");
  auto pde = Teuchos::rcp(new PDE_ElasticityCurvedFace(op_list, mesh));

  // populate boundary conditions: type (called model) and value
  int ndir(0), nshear(0), nkinematic(0);
  AmanziGeometry::Point val0(2), val1(2), val2(2);

  auto bcf = Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::VECTOR));
  std::vector<int>& bcf_model = bcf->bc_model();
  std::vector<std::vector<double>>& bcf_value = bcf->bc_value_vector(d);

  const auto fmap = mesh->getMap(AmanziMesh::Entity_kind::FACE, true);
  const auto bfmap = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true);

  for (int n = 0; n < bfmap.NumMyElements(); ++n) {
    int f = fmap.LID(bfmap.GID(n));
    const auto& xf = mesh->getFaceCentroid(f);
    auto yf = (*pde->get_bf())[f];

    bcf_model[f] = OPERATOR_BC_DIRICHLET;
    val0 = ana.velocity_exact(yf, 0.0);
    for (int k = 0; k < d; ++k) bcf_value[f][k] = val0[k];
    ndir++;
  }
  pde->AddBCs(bcf, bcf);

  // create and initialize solution
  const CompositeVectorSpace& cvs = pde->global_operator()->DomainMap();
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);

  // create source
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells_wghost; ++c) {
    const auto& xc = mesh->getCellCentroid(c);
    Point tmp(ana.source_exact(xc, 0.0));
    for (int k = 0; k < d; ++k) src[k][c] = tmp[k];
  }

  // populate the elasticity operator
  pde->Setup(K, false);
  pde->UpdateMatrices();

  // get and assemble the global operator
  Teuchos::RCP<Operator> global_op = pde->global_operator();
  global_op->UpdateRHS(source, true); // FIXME
  pde->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->set_inverse_parameters(
    "Hypre AMG", plist.sublist("preconditioners"), "PCG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  CompositeVector& rhs = *global_op->rhs();
  global_op->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "elasticity solver (PCG): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code() 
              << " size=" << rhs.GetLocalElementCount() << std::endl
              << "BCs: noslip: " << ndir << ", kinematic: " << nkinematic
              << ", shear stress: " << nshear << std::endl;
  }

  // compute velocity error
  double unorm, ul2_err, uinf_err;
  ana.VectorCellError(solution, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    ul2_err /= unorm;
    printf("L2(u)=%12.8g  Inf(u)=%12.8g  itr=%3d\n", ul2_err, uinf_err, global_op->num_itrs());
    // CHECK(ul2_err < 1e-10);
  }
}


TEST(OPERATOR_ELASTICITY_CURVED_FACE_2D)
{
  RunTest(1.0, 0.1);
}

TEST(OPERATOR_ELASTICITY_CURVED_FACE_3D)
{
  RunTest(1.0, 0.2, "test/random3D_05.exo");
  RunTest(1.0, 0.0, "test/sphere.exo");
}

