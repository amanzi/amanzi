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
#include "GMVMesh.hh"
#include "Tensor.hh"

// Operators
#include "Analytic01.hh"
#include "Analytic02.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_DiffusionMFD.hh"
#include "Verification.hh"


/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_NODAL)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, nodal discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 30, 30);
  RCP<const Mesh> mesh = meshfactory.create("test/median15x16.exo");

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);

  Analytic01 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Point xv(2);
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  const auto& fmap = mesh->getMap(AmanziMesh::Entity_kind::FACE,true);
  const auto& bmap = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);

  for (int bf = 0; bf < bmap.NumMyElements(); ++bf) {
    int f = fmap.LID(bmap.GID(bf));
    AmanziMesh::Entity_ID_List nodes;
    nodes = mesh->getFaceNodes(f);
    for (int n = 0; n < nodes.size(); ++n) {
      int v = nodes[n];
      xv = mesh->getNodeCoordinate(v);
      bc_model[v] = OPERATOR_BC_DIRICHLET;
      bc_value[v] = ana.pressure_exact(xv, 0.0);
    }
  }

  // create diffusion operator
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator nodal");
  auto op = Teuchos::rcp(new PDE_DiffusionMFD(op_list, mesh));
  op->Init(op_list);
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create source and add it to the operator
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("node", true);
  src.PutScalar(0.0);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    double volume = mesh->getCellVolume(c);

    AmanziMesh::Entity_ID_List nodes;
    nodes = mesh->getCellNodes(c);
    int nnodes = nodes.size();

    for (int k = 0; k < nnodes; k++) {
      int v = nodes[k];
      src[0][v] += ana.source_exact(xc, 0.0) * volume / nnodes;
    }
  }
  source.GatherGhostedToMaster();

  // populate the diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);
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
  ver.CheckPreconditionerSPD();
  ver.CheckSpectralBounds();

  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  global_op->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code()
              << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("node");
    GMV::open_data_file(*mesh, (std::string) "operators.gmv");
    GMV::start_data();
    GMV::write_node_data(p, 0, "solution");
    GMV::close_data_file();
  }

  CHECK(global_op->num_itrs() < 10);

  // compute pressure error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("node", false);

  double pnorm, pl2_err, pinf_err, hnorm, ph1_err;
  ana.ComputeNodeError(p, 0.0, pnorm, pl2_err, pinf_err, hnorm, ph1_err);

  if (MyPID == 0) {
    pl2_err /= pnorm;
    ph1_err /= hnorm;
    printf("L2(p)=%9.6f  H1(p)=%9.6f  itr=%3d\n", pl2_err, ph1_err, global_op->num_itrs());

    CHECK(pl2_err < 2e-2 && ph1_err < 7e-2);
    CHECK(global_op->num_itrs() < 10);
  }
}


/* *****************************************************************
* Exactness test for mixed diffusion solver.
* NOTE. Mixed boundary condition requires to use mass matrix. We
*       lump it which leads to a small error.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_NODAL_EXACTNESS)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0)
    std::cout << "\nTest: 2D elliptic solver, exactness"
              << " test for nodal discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 4, 4);
  RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the diffusion coefficient.
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::ALL);

  Analytic02 ana(mesh);

  for (int c = 0; c < ncells_wghost; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }
  double rho(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data (no mixed bc)
  Point xv(2);
  Teuchos::RCP<BCs> bc_v =
    Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model_v = bc_v->bc_model();
  std::vector<double>& bc_value_v = bc_v->bc_value();

  for (int v = 0; v < nnodes_wghost; v++) {
    xv = mesh->getNodeCoordinate(v);
    if (fabs(xv[0] - 1.0) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bc_model_v[v] = OPERATOR_BC_DIRICHLET;
      bc_value_v[v] = ana.pressure_exact(xv, 0.0);
    }
  }

  Teuchos::RCP<BCs> bc_f =
    Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model_f = bc_f->bc_model();
  std::vector<double>& bc_value_f = bc_f->bc_value();
  std::vector<double>& bc_mixed_f = bc_f->bc_mixed();

  int nn = 0;
  int nm = 0;
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6) {
      nn++;
      bc_model_f[f] = OPERATOR_BC_NEUMANN;
      bc_value_f[f] = -(ana.velocity_exact(xf, 0.0))[0]; // We assume exterior normal.
    } else if (fabs(xf[1]) < 1e-6) {
      nm++;
      bc_model_f[f] = OPERATOR_BC_MIXED;
      bc_value_f[f] = -(ana.velocity_exact(xf, 0.0))[1]; // We assume exterior normal.

      double tmp = ana.pressure_exact(xf, 0.0);
      bc_mixed_f[f] = 1.0;
      bc_value_f[f] -= bc_mixed_f[f] * tmp;
    }
  }

  // create diffusion operator
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator nodal");
  PDE_DiffusionFactory opfactory;
  Teuchos::RCP<PDE_Diffusion> op = opfactory.Create(op_list, mesh, bc_v, rho, g);
  op->AddBCs(bc_f, bc_f);

  // populate the diffusion operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->Init();
  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->set_inverse_parameters(
    "Hypre AMG", plist.sublist("preconditioners"), "AztecOO CG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();


  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  global_op->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code()
              << std::endl;
  }

  // compute pressure error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("node", false);

  double pnorm, pl2_err, pinf_err, hnorm, ph1_err;
  ana.ComputeNodeError(p, 0.0, pnorm, pl2_err, pinf_err, hnorm, ph1_err);

  if (MyPID == 0) {
    pl2_err /= pnorm;
    ph1_err /= hnorm;
    printf("L2(p)=%9.6f  H1(p)=%9.6f  itr=%3d\n", pl2_err, ph1_err, global_op->num_itrs());

    CHECK(pl2_err < 1e-5 && ph1_err < 2e-5);
    CHECK(global_op->num_itrs() < 10);
  }
}
