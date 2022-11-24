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
#include "CompositeVector.hh"
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "Analytic02.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionFracturedMatrix.hh"


/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* **************************************************************** */
void
TestDiffusionFracturedMatrix(double gravity)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: 3D fractured matrix problem: gravity=" << gravity << "\n";

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_fractured_matrix.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  ParameterList region_list = plist.sublist("regions");
  auto gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  // create a mesh framework
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // modify diffusion coefficient
  WhetStone::Tensor Knull;
  auto Kc = Teuchos::rcp(new std::vector<WhetStone::Tensor>());

  Analytic02 ana(mesh, Point(1.0, 2.0, 0.0), gravity, Knull);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Ktmp = ana.TensorDiffusivity(xc, 0.0);
    Kc->push_back(Ktmp);
  }

  // create boundary data.
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator");

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; ++f) {
    const Point& xf = mesh->face_centroid(f);

    // external boundary
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
        fabs(xf[1] - 1.0) < 1e-6 || fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }

    // internal boundary
    if (fabs(xf[2] - 0.5) < 1e-6) {
      bc_model[f] = OPERATOR_BC_NEUMANN;
      bc_value[f] = 0.0;
    }
  }

  // create diffusion operator
  double rho = 1.0;
  AmanziGeometry::Point gvec(3);
  gvec[2] = -gravity;

  if (gravity > 0.0) op_list.set<bool>("gravity", true);

  auto op = Teuchos::rcp(new PDE_DiffusionFracturedMatrix(op_list, mesh, rho, gvec));
  op->Init(op_list);
  auto global_op = op->global_operator();

  // -- boundary conditions
  op->SetBCs(bc, bc);

  // populate the diffusion operator
  op->Setup(Kc, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // apply BCs (primary=true, eliminate=true) and assemble
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  // ParameterList slist = plist.sublist("preconditioners").sublist("identity");
  global_op->set_inverse_parameters(
    "Hypre AMG", plist.sublist("preconditioners"), "AztecOO CG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();


  CompositeVector& rhs = *global_op->rhs();
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(rhs));
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(rhs));
  solution->PutScalar(0.0);

  global_op->ApplyInverse(rhs, *solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code()
              << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution->ViewComponent("cell");
    GMV::open_data_file(*mesh, (std::string) "operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::close_data_file();
  }

  CHECK(global_op->num_itrs() < 200);

  // compute pressure error
  Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // calculate flux error. To reuse the standard tools, we need to
  // collapse flux on fracture interface
  Epetra_MultiVector& flx_long = *flux->ViewComponent("face", true);
  Epetra_MultiVector flx_short(mesh->face_map(false), 1);
  double unorm, ul2_err, uinf_err;

  op->UpdateFlux(solution.ptr(), flux.ptr());

  const auto& fmap = *flux->Map().Map("face", true);
  for (int f = 0; f < nfaces; ++f) {
    int g = fmap.FirstPointInElement(f);
    flx_short[0][f] = flx_long[0][g];
  }

  ana.ComputeFaceError(flx_short, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  itr=%3d\n",
           pl2_err,
           pinf_err,
           ul2_err,
           uinf_err,
           global_op->num_itrs());

    CHECK(pl2_err < 1e-10);
    CHECK(ul2_err < 1e-10);
  }
}

TEST(OPERATOR_DIFFUSION_FRACTURED_MATRIX)
{
  TestDiffusionFracturedMatrix(0.0);
  TestDiffusionFracturedMatrix(9.8);
}
