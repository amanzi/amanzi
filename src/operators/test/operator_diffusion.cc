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
#include "EpetraExt_RowMatrixOut.h"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "Tensor.hh"

// Operators
#include "Analytic02.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionMFDwithGravity.hh"
#include "UpwindSecondOrder.hh"
#include "Verification.hh"


/* *****************************************************************
* Exactness test for mixed diffusion solver.
***************************************************************** */
void RunTestDiffusionMixed(int dim, double gravity, std::string pc_name = "Hypre AMG") {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: " << dim << "D elliptic solver, exactness test" 
                            << " for mixed discretization, g=" << gravity << std::endl;

  // read parameter list
  // -- it specifies details of the mesh, diffusion operator, and solver
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // provide at lest one framework to the mesh factory. The first available
  // -- framework will be used
  MeshFactory meshfactory(comm);
  RCP<const Mesh> mesh;
  if (dim == 2) {
    meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
    // mesh = meshfactory.create("test/median15x16.exo");
    mesh = meshfactory.create("test/circle_quad10.exo");
    // mesh = meshfactory.create("test/circle_poly10.exo");
  } else {
    meshfactory.set_preference(Preference({AmanziMesh::Framework::SIMPLE}));
    if (comm->NumProc() > 1) meshfactory.set_preference(Preference({Framework::MSTK}));
    mesh = meshfactory.create(0.0,0.0,0.0, 1.0,1.0,1.0, 4, 5, 6);
  }

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  Analytic02 ana(mesh, gravity);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();
  // std::vector<double>& bc_mixed = bc->bc_mixed();

  const auto& fmap = mesh->face_map(true);
  const auto& bmap = mesh->exterior_face_map(true);

  for (int bf = 0; bf < bmap.NumMyElements(); ++bf) {
    int f = fmap.LID(bmap.GID(bf));
    const Point& xf = mesh->face_centroid(f);
    double area = mesh->face_area(f);
    bool flag;
    Point normal = ana.face_normal_exterior(f, &flag);

    if (xf[0] < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;
/*
    } else if (xf[1] < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_MIXED;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;

      double tmp = ana.pressure_exact(xf, 0.0);
      bc_mixed[f] = 1.0;
      bc_value[f] -= bc_mixed[f] * tmp;
*/
    } else {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create diffusion operator 
  double rho(1.0);
  AmanziGeometry::Point g(dim);
  g[dim - 1] = -gravity;

  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator mixed");
  auto op = Teuchos::rcp(new PDE_DiffusionMFDwithGravity(op_list, mesh, rho, g));
  op->Init(op_list);
  op->SetBCs(bc, bc);

  // set up the diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // get and assemble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->set_inverse_parameters(pc_name, plist.sublist("preconditioners"),
          "AztecOO CG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // Test SPD properties of the preconditioner.
  VerificationCV ver(global_op);
  ver.CheckPreconditionerSPD();


  CompositeVector rhs = *global_op->rhs();
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(rhs));
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(rhs));
  solution->PutScalar(0.0);

  global_op->ApplyInverse(rhs, *solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual() 
              << " itr=" << global_op->num_itrs()
              << " code=" << global_op->returned_code() << std::endl;
  }

  // compute pressure error
  Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // calculate flux error
  Epetra_MultiVector& flx = *flux->ViewComponent("face", true);
  double unorm, ul2_err, uinf_err;

  op->UpdateFlux(solution.ptr(), flux.ptr());
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    ul2_err /= unorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  itr=%3d\n",
        pl2_err, pinf_err, ul2_err, uinf_err, global_op->num_itrs());

    CHECK(pl2_err < 1e-12 && ul2_err < 1e-12);
    if (pc_name != "identity") CHECK(global_op->num_itrs() < 10);
  }
}


TEST(OPERATOR_DIFFUSION_MIXED) {
  RunTestDiffusionMixed(2, 0.0, "identity");
  RunTestDiffusionMixed(3, 0.0);
}

TEST(OPERATOR_DIFFUSION_MIXED_wGRAVITY) {
  RunTestDiffusionMixed(2, 0.1);
}


/* *****************************************************************
* Exactness test for cell-based diffusion solver.
***************************************************************** */
TEST(OPERATOR_DIFFUSION_CELL_EXACTNESS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, exactness" 
                            << " test for cell-based discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a geometric model and mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 15, 8);

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Analytic02 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    WhetStone::Tensor Kc(2, 2);

    Kc(0, 0) = 3.0;
    Kc(1, 1) = 1.0;
    Kc(0, 1) = 0.0;
    Kc(1, 0) = 0.0;

    K->push_back(Kc);
  }

  // create boundary data.
  Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc_f->bc_model();
  std::vector<double>& bc_value = bc_f->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);

    if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = 3.0;
    /*
    } else if (fabs(xf[1]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_MIXED;
      bc_value[f] = 2.0;

      double tmp = ana.pressure_exact(xf, 0.0);
      bc_mixed[f] = 1.0;
      bc_value[f] -= bc_mixed[f] * tmp;
    */
    } else if (fabs(xf[0] - 1.0) < 1e-6 || 
               fabs(xf[1] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create diffusion operator 
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator cell");
  Teuchos::RCP<PDE_Diffusion> op = Teuchos::rcp(new PDE_DiffusionFV(op_list, mesh));
  op->SetBCs(bc_f, bc_f);

  // set up the diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // get and assmeble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  op->ApplyBCs(true, true, true);
  
  // create preconditoner using the base operator class
  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"),
          "AztecOO CG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  global_op->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual() 
              << " itr=" << global_op->num_itrs()
              << " code=" << global_op->returned_code() << std::endl;
  }

  // compute pressure error
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_err, pinf_err, global_op->num_itrs());

    CHECK(pl2_err < 1e-5);
    CHECK(global_op->num_itrs() < 10);
  }
}

