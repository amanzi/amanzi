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
#include "PDE_Abstract.hh"
#include "SurfaceCoordinateSystem.hh"

#include "AnalyticElasticity01.hh"
#include "Verification.hh"

/* *****************************************************************
* Elasticity model: exactness test.
***************************************************************** */
void
RunTest(int d, double mu, double lambda, const std::string& method)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0)
    std::cout << "\nTEST: " << d << "D elasticity exactness test for \"" << method << "\"\n";

  // read parameter list
  // -- it specifies details of the mesh, elasticity operator, and solver
  std::string xmlFileName = "test/operator_elasticity_weak_symmetry.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  plist.sublist("PK operator").sublist("elasticity operator").sublist("schema")
       .set<std::string>("method", method);

  // create the MSTK mesh framework
  // -- geometric model is not created. Instead, we specify boundary conditions
  // -- using centroids of mesh faces.
  auto mesh_list = Teuchos::rcp(new Teuchos::ParameterList());
  mesh_list->set<bool>("request faces", true);
  if (d == 3) mesh_list->set<bool>("request edges", true);
  MeshFactory meshfactory(comm, Teuchos::null, mesh_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  Teuchos::RCP<const Mesh> mesh;
  if (d == 2) {
    mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 5, 6);
  } else {
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 6, 4);
  }

  // -- general information about mesh
  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
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
  // -- XML list speficies discretization method and location of degrees of freedom
  // -- (called schema). This seems redundant but only when use a low-order method.
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("elasticity operator");
  Teuchos::RCP<PDE_Abstract> pde = Teuchos::rcp(new PDE_Abstract(op_list, mesh));

  // populate boundary conditions: type (called model) and value
  // -- normal component of velocity on boundary faces (a scalar)
  int ndir(0), nshear(0), nkinematic(0);
  AmanziGeometry::Point val0(d), val1(d), val2(d);

  Teuchos::RCP<BCs> bcf =
    Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::VECTOR));
  std::vector<int>& bcf_model = bcf->bc_model();
  int ndofs = (method == "elasticity weak symmetry BdV") ? d * d : d;
  std::vector<std::vector<double>>& bcf_value = bcf->bc_value_vector(ndofs);

  for (int f = 0; f < nfaces_wghost; ++f) {
    const auto& xf = mesh->getFaceCentroid(f);
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || 
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        (d == 3 && (fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6))) {
      bcf_model[f] = OPERATOR_BC_DIRICHLET;
      val0 = ana.velocity_exact(xf, 0.0);
      for (int k = 0; k < d; ++k) bcf_value[f][k] = val0[k];

      if (ndofs == d * d) {
        auto coordsys = std::make_shared<AmanziGeometry::SurfaceCoordinateSystem>(xf, normal);
        const auto& tau = (*coordsys->tau())[0];

        val1 = ana.velocity_exact(xf - tau / 2, 0.0);
        val2 = ana.velocity_exact(xf + tau / 2, 0.0);
        for (int k = 0; k < d; ++k) {
          bcf_value[f][k + 2] = (val2[k] - val1[k]) / 6 / 4; // FIXME
        }
      }

      ndir++;
    }
  }
  pde->AddBCs(bcf, bcf);

  // create and initialize solution
  const CompositeVectorSpace& cvs = pde->global_operator()->DomainMap();
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);

  // create source
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) {
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

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD(true, true);
  ver.CheckPreconditionerSPD(1e-12, true, true);

  CompositeVector& rhs = *global_op->rhs();
  global_op->ApplyInverse(rhs, solution);

  ver.CheckResidual(solution, 1.0e-13);

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

    CHECK(ul2_err < 1e-10);
    CHECK(global_op->num_itrs() < 20);
  }
}


TEST(OPERATOR_ELASTICITY_WEAK_SYMMETRY_EXACTNESS)
{
  RunTest(2, 1.0, 0.2, "elasticity weak symmetry BdV");
  RunTest(2, 1.0, 0.1, "elasticity weak symmetry");

  RunTest(3, 1.0, 0.1, "elasticity weak symmetry");
  // RunTest(3, 1.0, 0.1, "elasticity weak symmetry BdV");
}
