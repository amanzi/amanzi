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
#include "MFD3D_Diffusion_Edge.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"

// Operators
#include "Analytic02.hh"
#include "Analytic01b.hh"

#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "Verification.hh"

/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* The degrees of freedom are on mesh edges.
* **************************************************************** */
template <class Analytic>
void
TestDiffusionEdges(int dim, double tol, std::string filename)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0)
    std::cout << "\nTest: elliptic solver, edge discretization, d=" << dim << ", file: " << filename
              << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  Teuchos::RCP<GeometricModel> gm;
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK, Framework::STK }));
  RCP<const Mesh> mesh;
  if (dim == 2) mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 6, 7, true, true);
  // mesh = meshfactory.create("test/median32x33.exo", true, true);
  else
    mesh = meshfactory.create(filename, true, true);
  // mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 1, 1, true, true);

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  int nedges_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::EDGE, AmanziMesh::Parallel_type::ALL);

  Analytic ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc =
    Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::EDGE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  int d = dim - 1;
  for (int e = 0; e < nedges_wghost; ++e) {
    const Point& xe = mesh->getEdgeCentroid(e);
    if (fabs(xe[0]) < 1e-6 || fabs(xe[0] - 1.0) < 1e-6 || fabs(xe[1]) < 1e-6 ||
        fabs(xe[1] - 1.0) < 1e-6 || fabs(xe[d]) < 1e-6 || fabs(xe[d] - 1.0) < 1e-6) {
      bc_model[e] = OPERATOR_BC_DIRICHLET;
      bc_value[e] = ana.pressure_exact(xe, 0.0);
    }
  }

  // create diffusion operator
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator edge");
  auto op = Teuchos::rcp(new PDE_Abstract(op_list, mesh));
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create source and add it to the operator
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("edge", true);
  src.PutScalar(0.0);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    double volume = mesh->getCellVolume(c);

    AmanziMesh::Entity_ID_List edges;
    edges = mesh->getCellEdges(c);
    int nedges = edges.size();

    for (int k = 0; k < nedges; k++) {
      int e = edges[k];
      src[0][e] += ana.source_exact(xc, 0.0) * volume / nedges;
    }
  }
  source.GatherGhostedToMaster();

  // populate the diffusion operator
  op->Setup(K, false);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update the source term
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, true);

  // apply BCs (primary=true, eliminate=true) and assemble
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->set_inverse_parameters(
    "Hypre AMG", plist.sublist("preconditioners"), "AztecOO CG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // Test SPD properties of the preconditioner.
  VerificationCV ver(global_op);
  ver.CheckPreconditionerSPD(1.0e-11);

  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  global_op->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code()
              << std::endl;
  }
  CHECK(global_op->num_itrs() < 10);

  // compute error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("edge", false);

  double pnorm, pl2_err, pinf_err, hnorm, ph1_err;
  ana.ComputeEdgeError(p, 0.0, pnorm, pl2_err, pinf_err, hnorm, ph1_err);

  if (MyPID == 0) {
    pl2_err /= pnorm;
    ph1_err /= hnorm;
    printf("L2(p)=%9.6f  H1(p)=%9.6f  itr=%3d\n", pl2_err, ph1_err, global_op->num_itrs());

    CHECK(pl2_err < tol && ph1_err < tol);
  }
}


TEST(OPERATOR_DIFFUSION_EDGES)
{
  TestDiffusionEdges<Analytic02>(2, 1e-12, "");
  TestDiffusionEdges<Analytic02>(3, 1e-12, "test/tetrahedra.exo");
  TestDiffusionEdges<Analytic02>(3, 1e-12, "test/hexes.exo");
}
