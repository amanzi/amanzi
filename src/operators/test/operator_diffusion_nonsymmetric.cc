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
#include "EpetraExt_RowMatrixOut.h"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "Tensor.hh"

// Operators
#include "Analytic05.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"
#include "UpwindSecondOrder.hh"


/* *****************************************************************
* Non-symmetric diffusion tensor.
***************************************************************** */
TEST(OPERATOR_DIFFUSION_NONSYMMETRIC)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, non-symmetric tensor" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list =
    plist.sublist("PK operator").sublist("diffusion operator nonsymmetric");

  // create a mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 20, 20);

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the diffusion tensor
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  Analytic05 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc =
    Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    } else if (fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      double area = mesh->getFaceArea(f);
      const Point& normal = mesh->getFaceNormal(f);
      bc_model[f] = OPERATOR_BC_NEUMANN;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;
    }
  }

  // create diffusion operator
  auto op = Teuchos::rcp(new PDE_DiffusionMFD(op_list, mesh));
  op->Init(op_list);
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(cvs));
  solution->PutScalar(0.0);

  // create source
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    src[0][c] = ana.source_exact(xc, 0.0);
  }

  // populate the diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // get and assmeble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, false);
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->set_inverse_parameters(
    "Hypre AMG", plist.sublist("preconditioners"), "Belos GMRES", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  CompositeVector& rhs = *global_op->rhs();
  global_op->ApplyInverse(rhs, *solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (belos gmres): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code()
              << std::endl;
  }

  // compute pressure error
  Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // calculate flux error
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(cvs));
  Epetra_MultiVector& flx = *flux->ViewComponent("face", true);

  op->UpdateFlux(solution.ptr(), flux.ptr());
  double unorm, ul2_err, uinf_err;
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm;
    ul2_err /= unorm;
    printf("L2(p)=%12.8g  Inf(p)=%12.8g  L2(u)=%12.8g  Inf(u)=%12.8g  itr=%3d\n",
           pl2_err,
           pinf_err,
           ul2_err,
           uinf_err,
           global_op->num_itrs());

    CHECK(pl2_err < 0.03 && ul2_err < 0.1);
    CHECK(global_op->num_itrs() < 15);
  }
}
