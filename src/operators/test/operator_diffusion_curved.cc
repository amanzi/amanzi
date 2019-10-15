/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
#include "Tensor.hh"

// Operators
#include "Analytic02.hh"

#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "Verification.hh"


/* *****************************************************************
 * Exactness test for diffusion solver on meshes with curved faces.
 ***************************************************************** */
void
RunTestDiffusionCurved()
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int getRank = comm->getRank();
  if (getRank == 0)
    std::cout << "\nTest: elliptic solver, mesh with curved faces\n";

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_curved.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a randomized mesh
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK, Framework::STK }));
  // RCP<const Mesh> mesh = meshfactory.create(0.0,0.0,0.0, 1.0,1.0,1.0, 2,2,2);
  RCP<const Mesh> mesh = meshfactory.create("test/random3D_05.exo");

  // populate diffusion coefficient using the problem with analytic solution.
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned =
    mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned =
    mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost =
    mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Analytic02 ana(mesh);

  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // populate boundary data: The discretization method uses 3 DOFs (moment)
  // on each mesh face which require to specify 3 boundary data of type double
  // for each mesh face.
  Teuchos::RCP<BCs> bc =
    Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double>>& bc_value = bc->bc_value_vector(3);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    const Point& normal = mesh->face_normal(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
        fabs(xf[1] - 1.0) < 1e-6 || fabs(xf[2]) < 1e-6 ||
        fabs(xf[2] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f][0] = ana.pressure_exact(xf, 0.0);
      bc_value[f][1] = 0.0;
      bc_value[f][2] = 0.0;
    }
  }

  // create diffusion operator
  // -- use the abstract operator based on factory of MFD discretization
  // methods.
  ParameterList op_list =
    plist.sublist("PK operator").sublist("diffusion operator curved");
  Teuchos::RCP<PDE_Abstract> op = Teuchos::rcp(new PDE_Abstract(op_list, mesh));
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // -- set up diffusivity coefficient and populate local matrices.
  op->Setup(K);
  op->UpdateMatrices();

  // -- apply boundary conditions
  op->ApplyBCs(true, true, true);

  // -- assemble the global matrix
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create a preconditoner using the global matrix
  ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->InitializePreconditioner(slist);
  global_op->UpdatePreconditioner();

  // Test SPD properties of the preconditioner.
  VerificationCV ver(global_op);
  ver.CheckPreconditionerSPD();

  // solve the problem
  // -- create iand initialize a linear operator: PCG
  ParameterList lop_list =
    plist.sublist("solvers").sublist("PCG").sublist("pcg parameters");
  AmanziSolvers::
    LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  solver.Init(lop_list);

  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs), flux(rhs);
  solution.putScalar(0.0);

  // -- run PCG
  int ierr = solver.applyInverse(rhs, solution);

  if (getRank == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << solver.residual()
              << " itr=" << solver.num_itrs()
              << " code=" << solver.returned_code() << std::endl;
  }

  ver.CheckResidual(solution, 3.0e-7);

  // Post-processing
  // -- compute pressure error
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // -- coumpute flux error (work in progress)
  Epetra_MultiVector& flx = *flux.ViewComponent("face", true);
  double unorm(1.0), ul2_err(0.0), uinf_err(0.0);

  // op->UpdateFlux(solution, flux);
  // ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  double totvol(0.0);
  AmanziGeometry::Point center(3);
  for (int c = 0; c < ncells_owned; c++) {
    double vol = mesh->cell_volume(c, false);
    totvol += vol;
    center += mesh->cell_centroid(c) * vol;
  }
  double tmp_in[4], tmp_out[4] = { totvol, center[0], center[1], center[2] };
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 4, tmp_out, tmp_in);
  totvol = tmp_in[0];
  for (int i = 0; i < 3; ++i) center[i] = tmp_in[i + 1] / totvol;
  CHECK_CLOSE(1.0, totvol, 1e-12);

  if (getRank == 0) {
    std::cout << "Domain center:" << center << std::endl;
    std::cout << "Volume error: " << 1.0 - totvol << std::endl;

    pl2_err /= pnorm;
    ul2_err /= unorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  itr=%3d\n",
           pl2_err,
           pinf_err,
           ul2_err,
           uinf_err,
           solver.num_itrs());

    CHECK(pl2_err < 0.001 && ul2_err < 0.01);
    CHECK(solver.num_itrs() < 1000);
  }
}


TEST(OPERATOR_DIFFUSION_CURVED)
{
  RunTestDiffusionCurved();
}
