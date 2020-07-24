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
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "Tensor.hh"

// Operators
#include "Analytic03.hh"
#include "HeatConduction.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"
#include "UpwindFlux.hh"
#include "UpwindSecondOrder.hh"

#include "MyPDE_DiffusionMFD.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;

/* *****************************************************************
* Tests DivK diffusion solver with full tensor and source term.
* The model for kf is volime-weighted arithmetic average.
***************************************************************** */
template<class PDE, class UpwindClass>
void RunTestDiffusionDivK2D(std::string diffusion_list, std::string upwind_list) {
  auto comm = Amanzi::getDefaultComm();

  // parallel bug: twin component is used incorrectly in UpdateMatrices(). 
  // Scatter of little_k overrrides its ghost values. The subsequent 
  // algorithm uses the second item in the list returned by face_get_cells
  // as the twin component. We need to use global ids of cells for proper
  // ordering.  
  if (upwind_list == "upwind second-order" && comm->NumProc() > 1) return;

  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, divK discretization: \"" 
                            << diffusion_list << "\" + \"" << upwind_list << "\"\n";

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist(diffusion_list);

  // create an SIMPLE mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  // Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  std::string file = op_list.get<std::string>("file name", "test/random20.exo");
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(file);

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the diffusion tensor
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Analytic03 ana(mesh);

  const WhetStone::Tensor Kc(2, 1);
  Kc(0, 0) = 1.0;
  for (int c = 0; c < ncells; c++) K->push_back(Kc);

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create diffusion operator 
  auto op = Teuchos::rcp(new PDE(op_list, mesh));
  op->Init(op_list);
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(cvs));
  solution->PutScalar(0.0);

  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(cvs));
  Epetra_MultiVector& flx = *flux->ViewComponent("face", true);

  Point velocity(-1.0, 0.0);
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& normal = mesh->face_normal(f);
    flx[0][f] = velocity * normal;
  }
  
  // Create nonlinear coefficient.
  Teuchos::RCP<HeatConduction> knc = Teuchos::rcp(new HeatConduction(mesh));

  // Create upwind model
  Teuchos::ParameterList& ulist = plist.sublist("PK operator").sublist(upwind_list);
  UpwindClass upwind(mesh, knc);
  upwind.Init(ulist);

  knc->UpdateValues(*flux, bc_model, bc_value);  // 1st argument is not used
  upwind.Compute(*flux, *solution, bc_model, *knc->values());

  if (upwind_list == "upwind second-order") knc->UpdateValuesPostUpwind();

  // create source 
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    src[0][c] = ana.source_exact(xc, 0.0);
  }

  // populate the diffusion operator
  op->Setup(K, knc->values(), knc->derivatives());
  op->UpdateMatrices(flux.ptr());

  // get and assemble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, false);
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->InitializeInverse("Hypre AMG", plist.sublist("preconditioners"), "AztecOO CG", plist.sublist("solvers"));
  global_op->UpdateInverse();
  global_op->ComputeInverse();

  CompositeVector& rhs = *global_op->rhs();
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
  op->UpdateFlux(solution.ptr(), flux.ptr());
  double unorm, ul2_err, uinf_err;
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    ul2_err /= unorm;
    printf("L2(p)=%12.8g  Inf(p)=%12.8g  L2(u)=%12.8g  Inf(u)=%12.8g  itr=%3d\n",
        pl2_err, pinf_err, ul2_err, uinf_err, global_op->num_itrs());

    CHECK(pl2_err < 0.03 && ul2_err < 0.1);
    CHECK(global_op->num_itrs() < 10);
  }
}

TEST(OPERATOR_DIFFUSION_DIVK_AVERAGE_2D) {
  RunTestDiffusionDivK2D<PDE_DiffusionMFD, UpwindFlux<HeatConduction> >("diffusion operator divk", "upwind");
}

TEST(OPERATOR_DIFFUSION_DIVK_SECOND_ORDER) {
  RunTestDiffusionDivK2D<MyPDE_DiffusionMFD, UpwindSecondOrder<HeatConduction> >("diffusion operator second-order", "upwind second-order");
}


/* *****************************************************************
* Tests DivK diffusion solver with full tensor and source term.
* The model for kf is volime-weighted arithmetic average.
* 3D version
***************************************************************** */
TEST(OPERATOR_DIFFUSION_DIVK_AVERAGE_3D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 3D elliptic solver, divK discretization, average" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);
  // Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/mesh.exo");

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the nonlinear coefficient.
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Analytic03 ana(mesh);

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create diffusion operator 
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator divk");
  auto op = Teuchos::rcp(new PDE_DiffusionMFD(op_list, mesh));
  op->Init(op_list);
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(cvs));
  solution->PutScalar(0.0);

  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(cvs));
  Epetra_MultiVector& flx = *flux->ViewComponent("face", true);

  Point velocity(0.0, 0.0, 0.0);
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& normal = mesh->face_normal(f);
    flx[0][f] = velocity * normal;
  }
  
  // Create nonlinear coefficient.
  Teuchos::RCP<HeatConduction> knc = Teuchos::rcp(new HeatConduction(mesh));

  // Create upwind model
  Teuchos::ParameterList& ulist = plist.sublist("PK operator").sublist("upwind");
  UpwindFlux<HeatConduction> upwind(mesh, knc);
  upwind.Init(ulist);

  knc->UpdateValues(*flux, bc_model, bc_value);  // 1st argument is not used
  // upwind.Compute(*flux, *solution, bc_model, *knc->values());
  upwind.Compute(*flux, *solution, bc_model, *knc->values());

  // create source 
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    src[0][c] = ana.source_exact(xc, 0.0);
  }

  // populate the diffusion operator
  op->SetScalarCoefficient(knc->values(), knc->derivatives());
  op->UpdateMatrices(flux.ptr());

  // get and assmeble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->UpdateRHS(source, false);
  op->ApplyBCs(true, true, true);

  // create preconditoner using the base operator class
  global_op->InitializeInverse("Hypre AMG", plist.sublist("preconditioners"), "AztecOO CG", plist.sublist("solvers"));
  global_op->UpdateInverse();
  global_op->ComputeInverse();

  CompositeVector rhs = *global_op->rhs();
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
  op->UpdateFlux(solution.ptr(), flux.ptr());
  double unorm, ul2_err, uinf_err;
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    ul2_err /= unorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  itr=%3d\n",
        pl2_err, pinf_err, ul2_err, uinf_err, global_op->num_itrs());

    CHECK(pl2_err < 0.03 && ul2_err < 0.1);
    CHECK(global_op->num_itrs() < 10);
  }
}


