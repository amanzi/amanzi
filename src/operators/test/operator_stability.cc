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
#include "BCs.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"


/* *****************************************************************
* This test replaves tensor and boundary conditions by continuous
* functions. It analyzes accuracy of the MFD discretization with
* respect to scaling the stability term.
**************************************************************** */
TEST(OPERATOR_MIXED_DIFFUSION)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0)
    std::cout << "Test: 2D steady-state elliptic solver, mixed discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_stability.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 40, 40);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // create diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the diffusion coefficient.
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  Analytic01 ana(mesh);

  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  Point xv(2);
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
        fabs(xf[1] - 1.0) < 1e-6) {
      bc_value[f] = ana.pressure_exact(xf, 0.0);
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    }
  }

  // create space for diffusion operator
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  cvs->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(*cvs));

  // create source
  CompositeVector source(*cvs);
  source.PutScalarMasterAndGhosted(0.0);

  Epetra_MultiVector& src = *source.ViewComponent("cell");
  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    src[0][c] += ana.source_exact(xc, 0.0);
  }

  // MAIN LOOP
  for (int n = 0; n < 240; n += 50) {
    double factor = pow(10.0, (double)(n - 50) / 100.0) / 2;

    // create the local diffusion operator
    Teuchos::ParameterList olist = plist.sublist("PK operators").sublist("mixed diffusion");
    PDE_DiffusionMFD op2(olist, mesh);
    op2.Init(olist);
    op2.SetBCs(bc, bc);

    int schema_dofs = op2.schema_dofs();
    int schema_prec_dofs = op2.schema_prec_dofs();
    CHECK(schema_dofs ==
          (Operators::OPERATOR_SCHEMA_BASE_CELL + Operators::OPERATOR_SCHEMA_DOFS_FACE +
           Operators::OPERATOR_SCHEMA_DOFS_CELL));
    CHECK(schema_prec_dofs ==
          (Operators::OPERATOR_SCHEMA_DOFS_FACE + Operators::OPERATOR_SCHEMA_DOFS_CELL));

    op2.set_factor(factor); // for developers only
    op2.Setup(K, Teuchos::null, Teuchos::null);
    op2.UpdateMatrices();

    // get and assemeble the global operator
    Teuchos::RCP<Operator> global_op = op2.global_operator();
    global_op->UpdateRHS(source, false);
    op2.ApplyBCs(true, true, true);

    global_op->set_inverse_parameters(
      "Hypre AMG", plist.sublist("preconditioners"), "PCG", plist.sublist("solvers"));
    global_op->InitializeInverse();
    global_op->ComputeInverse();

    CompositeVector& rhs = *global_op->rhs();
    global_op->ApplyInverse(rhs, *solution);

    // calculate pressure errors
    Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
    double pnorm, pl2_err, pinf_err;
    ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

    // calculate flux errors
    Epetra_MultiVector& flx = *flux->ViewComponent("face", true);
    double unorm, ul2_err, uinf_err;

    op2.UpdateFlux(solution.ptr(), flux.ptr());
    flux->ScatterMasterToGhosted();
    ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

    if (MyPID == 0) {
      pl2_err /= pnorm;
      ul2_err /= unorm;
      printf("scale=%7.4g  L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f itr=%3d\n",
             factor,
             pl2_err,
             pinf_err,
             ul2_err,
             uinf_err,
             global_op->num_itrs());

      CHECK(pl2_err < 0.15 && ul2_err < 0.15);
    }
  }

  Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
  GMV::open_data_file(*mesh, (std::string) "operators.gmv");
  GMV::start_data();
  GMV::write_cell_data(p, 0, "pressure");
  GMV::close_data_file();
}


/* *****************************************************************
* This test replaces tensor and boundary conditions by continuous
* functions. It analyzed accuracy of the MFd discretization with
* respect to scaling of the stability term.
**************************************************************** */
TEST(OPERATOR_NODAL_DIFFUSION)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0)
    std::cout << "\nTest: 2D steady-state elliptic solver, nodal discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_stability.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 5, 5);
  // RCP<const Mesh> mesh = meshfactory.create("test/median255x256.exo");

  // create diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the diffusion coefficient.
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nnodes_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  Analytic01 ana(mesh);

  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data (no mixed bc)
  Point xv(2);
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int v = 0; v < nnodes_wghost; v++) {
    xv = mesh->getNodeCoordinate(v);
    if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 || fabs(xv[1]) < 1e-6 ||
        fabs(xv[1] - 1.0) < 1e-6) {
      bc_value[v] = ana.pressure_exact(xv, 0.0);
      bc_model[v] = Operators::OPERATOR_BC_DIRICHLET;
    }
  }

  // create diffusion operator
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("node", AmanziMesh::Entity_kind::NODE, 1);

  CompositeVector solution(*cvs);
  solution.PutScalar(0.0);

  // create source
  CompositeVector source(*cvs);
  source.PutScalarMasterAndGhosted(0.0);

  Epetra_MultiVector& src = *source.ViewComponent("node", true);
  for (int v = 0; v < nnodes_wghost; v++) {
    xv = mesh->getNodeCoordinate(v);
    src[0][v] = ana.source_exact(xv, 0.0);
  }

  // MAIN LOOP
  for (int n = 0; n < 400; n += 110) {
    // double factor = pow(10.0, (double)(n - 50) / 100.0) / 2;
    double factor = pow(10.0, (double)(n - 150) / 100.0) / 2;

    // create the local diffusion operator
    Teuchos::ParameterList olist = plist.sublist("PK operators").sublist("nodal diffusion");
    PDE_DiffusionMFD op2(olist, mesh);
    op2.Init(olist);
    op2.SetBCs(bc, bc);

    int schema_dofs = op2.schema_dofs();
    CHECK(schema_dofs ==
          (Operators::OPERATOR_SCHEMA_BASE_CELL | Operators::OPERATOR_SCHEMA_DOFS_NODE));

    op2.set_factor(factor); // for developers only
    op2.Setup(K, Teuchos::null, Teuchos::null);
    op2.UpdateMatrices(Teuchos::null, Teuchos::null);

    // get and assemeble the global operator
    Teuchos::RCP<Operator> global_op = op2.global_operator();
    global_op->UpdateRHS(source, false);
    op2.ApplyBCs(true, true, true);

    global_op->set_inverse_parameters(
      "Hypre AMG", plist.sublist("preconditioners"), "PCG", plist.sublist("solvers"));
    global_op->InitializeInverse();
    global_op->ComputeInverse();

    CompositeVector& rhs = *global_op->rhs();
    solution.PutScalar(0.0);
    int ierr = global_op->ApplyInverse(rhs, solution);
    CHECK(ierr == 0);

    // calculate errors
#ifdef HAVE_MPI
    solution.ScatterMasterToGhosted();
#endif

    Epetra_MultiVector& sol = *solution.ViewComponent("node", true);
    double pnorm, pl2_err, pinf_err, hnorm, ph1_err;
    ana.ComputeNodeError(sol, 0.0, pnorm, pl2_err, pinf_err, hnorm, ph1_err);

    if (MyPID == 0) {
      pl2_err /= pnorm;
      ph1_err /= hnorm;
      double tmp = op2.nfailed_primary() * 100.0 / ncells_owned;
      printf("scale=%7.4g  L2(p)=%9.6f  Inf(p)=%9.6f  H1(p)=%9.6g  itr=%3d  nfailed=%4.1f\n",
             factor,
             pl2_err,
             pinf_err,
             ph1_err,
             global_op->num_itrs(),
             tmp);

      CHECK(pl2_err < 0.1 && ph1_err < 0.15);
    }
  }
}
