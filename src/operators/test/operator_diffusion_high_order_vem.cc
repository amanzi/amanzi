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
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "SurfaceCoordinateSystem.hh"

// Operators
#include "Analytic07b.hh"
#include "MeshDeformation.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "VEM_Diffusion_HighOrder.hh"


/* *****************************************************************
* Exactness test for high-order Raviart-Thomas element.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_HIGH_ORDER_RAVIART_THOMAS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  using namespace Amanzi::WhetStone;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, high-order Raviart-Thomas" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  Teuchos::RCP<GeometricModel> gm;
  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5, true, true);
  DeformMesh(mesh, 5, 1.0);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // model
  Analytic07b ana(mesh);

  // create boundary data
  ParameterList op_list = plist.sublist("PK operator")
                               .sublist("diffusion operator Raviart-Thomas");
  int order = op_list.sublist("schema").get<int>("method order");
  int nd = PolynomialSpaceDimension(2, order);

  Point xv(3), x0(3), x1(3);
  AmanziMesh::Entity_ID_List nodes;

  Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, DOF_Type::VECTOR));
  std::vector<int>& bc_model_f = bc_f->bc_model();
  std::vector<std::vector<double> >& bc_value_f = bc_f->bc_value_vector(nd);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) {

      bc_model_f[f] = OPERATOR_BC_DIRICHLET;
      bc_value_f[f][0] = ana.pressure_exact(xf, 0.0);
      bc_value_f[f][1] = 0.0;
      bc_value_f[f][2] = 0.0;
    }
  }

  RCP<const Mesh> mesh_tmp(mesh);
  Teuchos::RCP<PDE_Abstract> op = Teuchos::rcp(new PDE_Abstract(op_list, mesh_tmp));
  op->AddBCs(bc_f, bc_f);

  // create source term
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->Init();

  auto& rhs = *global_op->rhs();
  auto& rhs_c = *rhs.ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);
    rhs_c[0][c] = ana.source_exact(xc, 0.0) * volume;
  }

  // populate the diffusion operator
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs(true, true, true);

  // create preconditioner using the base operator class
  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"), "AztecOO CG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // solve the problem
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);
  global_op->ApplyInverse(rhs, solution);

  // compute pressure error
  auto& lambda = *solution.ViewComponent("face");
  auto& p = *solution.ViewComponent("cell");
  double pnorm, pl2_err, pinf_err;

  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // compute flux and flux error. Since, the abstract operator does not support
  // flux calculation at the moment, we delegate this to the discretization class
  Epetra_MultiVector flux(lambda);
  double unorm, ul2_err, uinf_err;

  ParameterList vem_list = op_list.sublist("schema");
  VEM_Diffusion_HighOrder vem(vem_list, mesh);

  double **pf, **uf;
  lambda.ExtractView(&pf);
  flux.ExtractView(&uf);

  for (int c = 0; c < ncells_owned; ++c) {
    const auto& xc = mesh->cell_centroid(c);
    auto K = ana.TensorDiffusivity(xc, 0.0);
    vem.UpdateFlux(c, K, pf, p[0], uf);
  }

  ana.ComputeFaceError(flux, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    ul2_err /= unorm;  // this is l2 norm (little l)
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6f  Inf(u)=%9.6f  size=%d  itr=%d\n", 
        pl2_err, pinf_err, ul2_err, uinf_err, global_op->A()->NumGlobalRows(), global_op->num_itrs());

    CHECK(pl2_err < 1e-2);
  }
}

