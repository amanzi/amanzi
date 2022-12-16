/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Tests for the diffusion solver on a fracture network.
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
#include "PDE_DiffusionFVonManifolds.hh"


/* *****************************************************************
* TBW.
* **************************************************************** */
void
RunTest(int icase, double gravity)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0)
    std::cout << "\nTest: FV scheme for diffusion in fractures g=" << gravity << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_dfn_fv.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an SIMPLE mesh framework
  ParameterList region_list = plist->sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, true, true);
  std::string setname("fractures");
  RCP<const Mesh> surfmesh = Teuchos::rcp(
    new MeshExtractedManifold(mesh, setname, AmanziMesh::FACE, comm, gm, plist, true, false));

  // modify diffusion coefficient
  int ncells_owned = surfmesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = surfmesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = surfmesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  double rho(1.0);
  AmanziGeometry::Point v(3);
  Analytic00b ana(surfmesh, 1.0, 2.0, 3.0, 1, v, gravity);

  // create boundary data (Dirichlet everywhere)
  Teuchos::RCP<BCs> bc =
    Teuchos::rcp(new BCs(surfmesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = surfmesh->face_centroid(f);
    if (fabs(xf[2] - 0.5) < 1e-6 && (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
                                     fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6)) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    } else if (fabs(xf[1] - 0.5) < 1e-6 && (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
                                            fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6)) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create solution
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(surfmesh)->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1)->SetOwned(false);

  auto solution = Teuchos::rcp(new CompositeVector(*cvs));
  solution->PutScalar(0.0);

  // create diffusion operator
  AmanziGeometry::Point gvec(0.0, 0.0, -gravity);
  Teuchos::ParameterList olist = plist->sublist("PK operator").sublist("diffusion operator");
  olist.set<bool>("gravity", (gravity > 0.0));

  auto op = Teuchos::rcp(new Operators::PDE_DiffusionFVonManifolds(olist, surfmesh));
  op->SetBCs(bc, bc);

  auto cvs2 = Operators::CreateManifoldCVS(surfmesh);
  if (icase == 0) {
    op->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  } else if (icase == 1) {
    auto k = Teuchos::rcp(new CompositeVector(*cvs2));
    k->PutScalar(1.0);
    op->SetScalarCoefficient(k, Teuchos::null);
    op->SetDensity(rho);
    op->SetGravity(gvec);
  }

  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->Init();

  // populate diffusion operator
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs(true, true, true);

  // create preconditoner
  global_op->set_inverse_parameters(
    "Hypre AMG", plist->sublist("preconditioners"), "PCG", plist->sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  CompositeVector rhs = *global_op->rhs();
  global_op->ApplyInverse(rhs, *solution);

  // post-processing
  auto flux = Teuchos::rcp(new CompositeVector(*cvs2));
  op->UpdateFlux(solution.ptr(), flux.ptr());

  // statistics
  int ndofs = global_op->A()->NumGlobalRows();
  double a;
  rhs.Norm2(&a);
  if (MyPID == 0) {
    std::cout << "pressure solver"
              << ": ||r||=" << global_op->residual() << " itr=" << global_op->num_itrs()
              << "  ||f||=" << a << "  #dofs=" << ndofs << " code=" << global_op->returned_code()
              << std::endl;
  }

  // calculate error in potential
  double pnorm, l2_err, inf_err;
  Epetra_MultiVector& p = *solution->ViewComponent("cell");

  ana.ComputeCellError(p, 0.0, pnorm, l2_err, inf_err);
  CHECK(l2_err < 1e-12 * pnorm);

  // calculate flux error. To reuse the standard tools, we need to
  // collapse flux on fracture interface
  double unorm, ul2_err, uinf_err;
  Epetra_MultiVector& flx_long = *flux->ViewComponent("face", true);
  Epetra_MultiVector flx_short(surfmesh->face_map(false), 1);

  const auto& fmap = *flux->Map().Map("face", true);
  for (int f = 0; f < nfaces_owned; ++f) {
    int g = fmap.FirstPointInElement(f);
    flx_short[0][f] = flx_long[0][g];
  }

  ana.ComputeFaceError(flx_short, 0.0, unorm, ul2_err, uinf_err);
  CHECK(ul2_err < 1e-12 * unorm);

  if (MyPID == 0) {
    l2_err /= pnorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  norms=%9.6f %9.6f\n",
           l2_err,
           inf_err,
           ul2_err,
           uinf_err,
           pnorm,
           unorm);
  }
}


TEST(DIFFUSION_FRACTURES_FV_NO_K)
{
  RunTest(0, 0.0);
}

TEST(DIFFUSION_FRACTURES_FV_K)
{
  RunTest(1, 1.0);
}
