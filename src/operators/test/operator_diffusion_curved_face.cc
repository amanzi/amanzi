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

// Operators
#include "Analytic00.hh"
#include "Analytic02.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionCurvedFace.hh"


/* *****************************************************************
* Exactness test for diffusion solver on meshes with curved faces.
***************************************************************** */
void
RunTestDiffusionCurved(const std::string& filename, int icase)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: elliptic solver, mesh with curved faces, new algorithm\n";

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_curved_face.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a randomized mesh
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<const Mesh> mesh = meshfactory.create(0.0,0.0,0.0, 1.0,1.0,1.0, 2,2,2);
  RCP<const Mesh> mesh = meshfactory.create(filename);

  // populate diffusion coefficient using the problem with analytic solution.
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Analytic00 ana(mesh, 1);
  // Analytic02 ana(mesh);

  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create diffusion operator
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion curved face");
  auto op = Teuchos::rcp(new PDE_DiffusionCurvedFace(op_list, mesh));

  // -- boundary data
  int nneu(0);
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  const auto fmap = mesh->face_map(true);
  const auto bfmap = mesh->exterior_face_map(true);

  // for (int f = 0; f < nfaces_wghost; ++f) {
  for (int n = 0; n < bfmap.NumMyElements(); ++n) {
    int f = fmap.LID(bfmap.GID(n));
    const Point& xf = mesh->face_centroid(f);

    auto yf = (*op->get_bf())[f];
    bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value[f] = ana.pressure_exact(yf, 0.0);

    // overwrite boundary data on case-by-case basis
    if (icase == 1) {
      if (fabs(xf[1]) < 1e-6) {
        double area = mesh->face_area(f);
        const Point& normal = mesh->face_normal(f);
        bc_model[f] = OPERATOR_BC_NEUMANN;
        bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;
        nneu++;
      }
    }
  }
  op->SetBCs(bc, bc);

  // -- set up diffusivity coefficient and populate local matrices.
  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices();

  // -- apply boundary conditions
  op->ApplyBCs(true, true, true);

  // -- assemble the global matrix
  Teuchos::RCP<Operator> global_op = op->global_operator();

  // create a preconditoner using the global matrix
  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // create a preconditoner using the global matrix
  global_op->set_inverse_parameters(
    "Hypre AMG", plist.sublist("preconditioners"), "PCG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs), flux(rhs);
  solution.PutScalar(0.0);

  // -- run PCG
  global_op->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code()
              << ",  Neumann BCs: " << nneu << std::endl;
  }

  // Post-processing
  // -- compute pressure error
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // -- coumpute flux error (work in progress)
  // Epetra_MultiVector& flx = *flux.ViewComponent("face", true);
  double unorm(1.0), ul2_err(0.0), uinf_err(0.0);
  // op->UpdateFlux(solution, flux);
  // ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm;
    ul2_err /= unorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  itr=%3d\n",
           pl2_err,
           pinf_err,
           ul2_err,
           uinf_err,
           global_op->num_itrs());

    CHECK(pl2_err < 1e-10 && ul2_err < 1e-10);
    CHECK(global_op->num_itrs() < 1000);
  }
}


TEST(OPERATOR_DIFFUSION_CURVED)
{
  RunTestDiffusionCurved("test/random3D_05.exo", 1);
  RunTestDiffusionCurved("test/hexes.exo", 1);
  RunTestDiffusionCurved("test/sphere.exo", 0);
}
