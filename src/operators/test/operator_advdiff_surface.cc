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
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_AdvectionUpwind.hh"

#include "Verification.hh"


/* *****************************************************************
* This test replaces tensor and boundary conditions by continuous
* functions. This is a prototype for future solvers.
* **************************************************************** */
TEST(ADVECTION_DIFFUSION_SURFACE) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Advection-duffusion on a surface" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_advdiff_surface.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an MSTK mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 40, 40, 5);
  RCP<const Mesh_MSTK> mesh_mstk = rcp_static_cast<const Mesh_MSTK>(mesh);

  // extract surface mesh
  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top surface"));

  RCP<Mesh> surfmesh = meshfactory.create(mesh_mstk, setnames, AmanziMesh::FACE);

  /* modify diffusion coefficient */
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = surfmesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = surfmesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(surfmesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = surfmesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = xf[1] * xf[1];
    }
  }

  // create diffusion operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("diffusion operator");
  auto op_diff = Teuchos::rcp(new PDE_DiffusionMFD(olist, (Teuchos::RCP<const AmanziMesh::Mesh>) surfmesh));
  op_diff->Init(olist);
  op_diff->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op_diff->global_operator()->DomainMap();

  // set up the diffusion operator
  op_diff->Setup(K, Teuchos::null, Teuchos::null);
  op_diff->UpdateMatrices(Teuchos::null, Teuchos::null);

  // get the global operator
  Teuchos::RCP<Operator> global_op = op_diff->global_operator();

  // create an advection operator  
  Teuchos::ParameterList alist;
  Teuchos::RCP<PDE_AdvectionUpwind> op_adv = Teuchos::rcp(new PDE_AdvectionUpwind(alist, global_op));
  op_adv->SetBCs(bc, bc);

  // get a flux field
  Teuchos::RCP<CompositeVector> u = Teuchos::rcp(new CompositeVector(cvs));
  Epetra_MultiVector& uf = *u->ViewComponent("face");
  int nfaces = surfmesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  Point vel(4.0, 4.0, 0.0);
  for (int f = 0; f < nfaces; f++) {
    uf[0][f] = vel * surfmesh->face_normal(f);
  }

  op_adv->Setup(*u);
  op_adv->UpdateMatrices(u.ptr());

  // Add an accumulation term.
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);  // solution at time T=0

  CompositeVector phi(cvs);
  phi.PutScalar(0.2);

  double dT = 0.02;
  Teuchos::RCP<PDE_Accumulation> op_acc = Teuchos::rcp(new PDE_Accumulation(AmanziMesh::CELL, global_op));
  op_acc->AddAccumulationDelta(solution, phi, phi, dT, "cell");

  // BCs and assemble
  op_diff->ApplyBCs(true, true, true);
  op_adv->ApplyBCs(true, true, true);

  // Create a preconditioner.
  ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"), "AztecOO CG", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD(false, true);
  ver.CheckPreconditionerSPD(1e-12, false, true);

  CompositeVector& rhs = *global_op->rhs();
  global_op->ApplyInverse(rhs, solution);

  int num_itrs = global_op->num_itrs();
  CHECK(num_itrs > 5 && num_itrs < 15);

  ver.CheckResidual(solution, 1.0e-12);

  if (MyPID == 0) {
    std::cout << "pressure solver (gmres): ||r||=" << global_op->residual() 
              << " itr=" << global_op->num_itrs()
              << " code=" << global_op->returned_code() << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("cell");
    GMV::open_data_file(*surfmesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::close_data_file();
  }
}
