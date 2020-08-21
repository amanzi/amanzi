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

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_DiffusionNLFVwithGravity.hh"
#include "PDE_DiffusionNLFVwithBndFacesGravity.hh"

#include "Analytic01.hh"
#include "Analytic02.hh"

/* *****************************************************************
* Nonlinear finite volume scheme.
***************************************************************** */
template<class Analytic>
void RunTestDiffusionNLFV_DMP(double gravity, bool testing) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, NLFV with DMP, g=" << gravity << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator nlfv");

  // create an MSTK mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  // Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 4, 4);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/random10.exo");

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Analytic ana(mesh, gravity);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    } else if (fabs(xf[1] - 1.0) < 1e-6) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      double area = mesh->face_area(f);
      bc_model[f] = OPERATOR_BC_NEUMANN;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;
    }
  }

  // create diffusion operator 
  double rho(1.0);
  AmanziGeometry::Point g(0.0, -gravity);
  Teuchos::RCP<PDE_Diffusion> op = Teuchos::rcp(new PDE_DiffusionNLFVwithGravity(op_list, mesh, rho, g));
  Teuchos::RCP<Operator> global_op = op->global_operator();
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = global_op->DomainMap();

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(cvs));
  solution->PutScalar(0.0);

  // create source 
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    src[0][c] = ana.source_exact(xc, 0.0);
  }

  // populate the diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);
  for (int loop = 0; loop < 12; ++loop) {
    global_op->Init();
    op->UpdateMatrices(Teuchos::null, solution.ptr());

    // get and assmeble the global operator
    global_op->UpdateRHS(source, false);
    op->ApplyBCs(true, true, true);

    // create preconditoner using the base operator class
    global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"),
            "Belos GMRES", plist.sublist("solvers"));
    global_op->InitializeInverse();
    global_op->ComputeInverse();

    CompositeVector& rhs = *global_op->rhs();
    global_op->ApplyInverse(rhs, *solution);
    
    // compute pressure error
    Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
    double pnorm, pl2_err, pinf_err;
    ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

    // compute flux error
    CompositeVectorSpace cvs_tmp;
    cvs_tmp.SetMesh(mesh)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("face", AmanziMesh::FACE, 1);
    Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(cvs_tmp));
    Epetra_MultiVector& flx = *flux->ViewComponent("face", true);
   
    op->UpdateFlux(solution.ptr(), flux.ptr());
    double unorm, ul2_err, uinf_err;

    ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

    if (MyPID == 0) {
      pl2_err /= pnorm; 
      ul2_err /= unorm;
      printf("L2(p)=%10.4e  Inf(p)=%10.4e  L2(u)=%10.4e  Inf(u)=%10.4e  (itr=%d ||r||=%10.4e code=%d)\n",
          pl2_err, pinf_err, ul2_err, uinf_err,
          global_op->num_itrs(), global_op->residual(), global_op->returned_code());

      if (testing) {
        double factor = std::pow(1.6, loop);
        CHECK(pl2_err < 0.2 / factor && ul2_err < 0.4 / factor);
      }
      CHECK(global_op->num_itrs() < 15);
    }
  }
}


template<class Analytic>
void RunTestDiffusionNLFVwithBndFaces_DMP(double gravity, bool testing) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, NLFVwithBndFaces with DMP, g=" << gravity << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator")
                                        .sublist("diffusion operator nlfv");

  // create an MSTK mesh framework
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  //Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 4, 4);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/random10.exo");

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Analytic ana(mesh, gravity);
  Analytic ana_diff(mesh, 0.0);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    } else if (fabs(xf[1] - 1.0) < 1e-6) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      double area = mesh->face_area(f);
      bc_model[f] = OPERATOR_BC_NEUMANN;
      bc_value[f] = ana_diff.velocity_exact(xf, 0.0) * normal / area;
    }
  }

  // create diffusion operator 
  double rho(1.0);
  AmanziGeometry::Point g(0.0, -gravity);
  Teuchos::RCP<PDE_Diffusion> op = Teuchos::rcp(new PDE_DiffusionNLFVwithBndFacesGravity(op_list, mesh, rho, g));
  Teuchos::RCP<Operator> global_op = op->global_operator();
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = global_op->DomainMap();

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(cvs));
  solution->PutScalar(0.0);
 
  // create source 
  CompositeVector source(cvs), t1(cvs), t2(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    src[0][c] = ana.source_exact(xc, 0.0);
  }

  // populate the diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);

  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"),
          "Belos GMRES", plist.sublist("solvers"));
  global_op->InitializeInverse();
  
  for (int loop = 0; loop < 12; ++loop) {
    global_op->Init();
    op->UpdateMatrices(Teuchos::null, solution.ptr());
    // get and assmeble the global operator
    global_op->UpdateRHS(source, false);
    op->ApplyBCs(true, true, true);

    // create preconditoner using the base operator class
    global_op->ComputeInverse();

    CompositeVector& rhs = *global_op->rhs();
    global_op->ApplyInverse(rhs, *solution);
 
    // compute pressure error
    Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
    double pnorm, pl2_err, pinf_err;
    ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

    // compute flux error
    CompositeVectorSpace cvs_tmp;
    cvs_tmp.SetMesh(mesh)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("face", AmanziMesh::FACE, 1);
    Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(cvs_tmp));
    Epetra_MultiVector& flx = *flux->ViewComponent("face", true);

    op->UpdateFlux(solution.ptr(), flux.ptr());
    double unorm, ul2_err, uinf_err;

    ana_diff.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

    if (MyPID == 0) {
      pl2_err /= pnorm; 
      ul2_err /= unorm;
      printf("L2(p)=%10.4e  Inf(p)=%10.4e  L2(u)=%10.4e  Inf(u)=%10.4e  (itr=%d ||r||=%10.4e code=%d)\n",
          pl2_err, pinf_err, ul2_err, uinf_err,
          global_op->num_itrs(), global_op->residual(), global_op->returned_code());

      if (testing) {
        double factor = 0.5*std::pow(1.6, loop);
        CHECK(pl2_err < 0.2 / factor && ul2_err < 0.4 / factor);
      }
      CHECK(global_op->num_itrs() < 15);
    }
  }
}

TEST(OPERATOR_DIFFUSION_NLFV_DMP_02) {
  RunTestDiffusionNLFV_DMP<Analytic02>(0.0, true);
}

TEST(OPERATOR_DIFFUSION_NLFV_wGravity) {
  RunTestDiffusionNLFV_DMP<Analytic02>(2.7, true);
}

TEST(OPERATOR_DIFFUSION_NLFV_DMP_01) {
  RunTestDiffusionNLFV_DMP<Analytic01>(2.7, false);
}

TEST(OPERATOR_DIFFUSION_NLFVwithBndFaces_DMP_02) {
  RunTestDiffusionNLFVwithBndFaces_DMP<Analytic02>(0.0, true);
}

TEST(OPERATOR_DIFFUSION_NLFVwithBndFaces_wGravity) {
  RunTestDiffusionNLFVwithBndFaces_DMP<Analytic02>(2.7, true);
}

TEST(OPERATOR_DIFFUSION_NLFVwithBndFaces_DMP_01) {
  RunTestDiffusionNLFVwithBndFaces_DMP<Analytic01>(2.7, false);
}
