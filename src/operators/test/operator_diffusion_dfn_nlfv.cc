/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Tests for diffusion solver on manifolds, nonlinear FV schemes.
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "Analytic00b.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionNLFV.hh"


/* *****************************************************************
* TBW.
* **************************************************************** */
void
RunTest(int icase, double gravity, int nx = 10, double tol = 1e-12)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0)
    std::cout << "\nTest: FV scheme for diffusion in fractures, icase=" << icase << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_dfn_fv.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an SIMPLE mesh framework
  ParameterList region_list = plist->sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  auto flist = Teuchos::rcp(new Teuchos::ParameterList());
  flist->set<bool>("request faces", true);
  flist->set<bool>("request edges", true);
  MeshFactory meshfactory(comm, gm, flist);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx);
  std::string setname = "fracture 1";
  auto surfmesh_fw =
    Teuchos::rcp(new MeshExtractedManifold(mesh, setname, AmanziMesh::FACE, comm, gm, plist));
  Teuchos::RCP<Mesh> surfmesh = Teuchos::rcp(
    new Mesh(surfmesh_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));

  // modify diffusion coefficient
  int ncells_owned = surfmesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_owned = surfmesh->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_wghost = surfmesh->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);

  AmanziGeometry::Point v(3);
  Teuchos::RCP<AnalyticBase> ana;
  ana = Teuchos::rcp(new Analytic00b(surfmesh, 1.0, 2.0, 0.0, 1, v, gravity));

  // create boundary data (Dirichlet everywhere)
  auto bc = Teuchos::rcp(new BCs(surfmesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = surfmesh->getFaceCentroid(f);
    if (fabs(xf[2] - 0.5) < 1e-6 && (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
                                     fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6)) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana->pressure_exact(xf, 0.0);
    } else if (fabs(xf[1] - 0.5) < 1e-6 && (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
                                            fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6)) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana->pressure_exact(xf, 0.0);
    }
  }

  // allocate memory for solution
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(surfmesh)->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1)->SetOwned(false);

  auto solution = Teuchos::rcp(new CompositeVector(*cvs));
  solution->PutScalar(0.0);

  // create diffusion operator
  AmanziGeometry::Point gvec(0.0, 0.0, -gravity);
  Teuchos::ParameterList olist = plist->sublist("PK operator").sublist("diffusion operator");
  olist.set<bool>("gravity", (gravity > 0.0));

  auto op = Teuchos::rcp(new Operators::PDE_DiffusionNLFV(olist, surfmesh));
  op->SetBCs(bc, bc);
  op->SetScalarCoefficient(Teuchos::null, Teuchos::null);

  auto cvs2 = Operators::CreateManifoldCVS(surfmesh);
  auto flux = Teuchos::rcp(new CompositeVector(*cvs2));

  // create optional source term
  CompositeVector src(*cvs);
  auto& src_c = *src.ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = surfmesh->getCellCentroid(c);
    src_c[0][c] = ana->source_exact(xc, 0.0);
  }

  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->set_inverse_parameters(
    "Hypre AMG", plist->sublist("preconditioners"), "PCG", plist->sublist("solvers"));

  for (int loop = 0; loop < 2; ++loop) {
    // populate diffusion operator
    global_op->Init();

    // add source term
    global_op->UpdateRHS(src, false);

    op->UpdateMatrices(Teuchos::null, solution.ptr());
    op->ApplyBCs(true, true, true);

    // create preconditoner
    global_op->InitializeInverse();
    global_op->ComputeInverse();

    CompositeVector rhs = *global_op->rhs();
    global_op->ApplyInverse(rhs, *solution);
  }

  // post-processing
  op->UpdateFlux(solution.ptr(), flux.ptr());

  // statistics
  int ndofs = global_op->A()->NumGlobalRows();
  double fnorm;
  global_op->rhs()->Norm2(&fnorm);
  if (MyPID == 0) {
    std::cout << "pressure solver" << ": ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << "  ||f||=" << fnorm << "  #dofs=" << ndofs
              << std::endl;
  }

  // calculate error in potential
  double pnorm, l2_err, inf_err;
  Epetra_MultiVector& p = *solution->ViewComponent("cell");

  ana->ComputeCellError(p, 0.0, pnorm, l2_err, inf_err);
  CHECK(l2_err < tol * pnorm);

  // calculate flux error
  double unorm, ul2_err, uinf_err;
  Epetra_MultiVector& flux_f = *flux->ViewComponent("face", true);

  ana->ComputeFaceError(flux_f, 0.0, unorm, ul2_err, uinf_err);
  // CHECK(ul2_err < tol * unorm);

  if (MyPID == 0) {
    l2_err /= pnorm;
    printf("rel norms: L2(p)=%9.6f Inf(p)=%9.6f L2(u)=%9.6g Inf(u)=%9.6f norms(p/u): %9.6f %9.6f\n",
           l2_err / pnorm,
           inf_err / pnorm,
           ul2_err / unorm,
           uinf_err / unorm,
           pnorm,
           unorm);
  }
}


TEST(DIFFUSION_FRACTURES_NLFV_NO_K)
{
  RunTest(0, 0.0);
}
