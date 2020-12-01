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

// Operators
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "PDE_Accumulation.hh"
#include "PDE_DiffusionMFD.hh"

#include "Verification.hh"

namespace Amanzi{

class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) { 
    CompositeVectorSpace cvs;
    auto cmap = Amanzi::getMaps(*mesh_, AmanziMesh::CELL);
    auto fmap = Amanzi::getMaps(*mesh_, AmanziMesh::FACE);
    cvs.SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, cmap.first, cmap.second, 1)
      ->AddComponent("face", AmanziMesh::FACE, fmap.first, fmap.second, 1);

    values_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
    derivatives_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
  }
  ~HeatConduction() {};

  // main members
  void UpdateValues(const CompositeVector& u) { 
    const Epetra_MultiVector& uc = *u.ViewComponent("cell", true); 
    const Epetra_MultiVector& values_c = *values_->ViewComponent("cell", true); 

    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    for (int c = 0; c < ncells; c++) {
      values_c[0][c] = 0.3 + uc[0][c];
    }

    const Epetra_MultiVector& uf = *u.ViewComponent("face", true); 
    const Epetra_MultiVector& values_f = *values_->ViewComponent("face", true); 
    int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    for (int f = 0; f < nfaces; f++) {
      values_f[0][f] = 0.3 + uf[0][f];
    }

    derivatives_->PutScalar(1.0);
  }

  Teuchos::RCP<CompositeVector> values() { return values_; }
  Teuchos::RCP<CompositeVector> derivatives() { return derivatives_; }
   

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<CompositeVector> values_, derivatives_;
};

}  // namespace Amanzi


namespace {

/* *****************************************************************
* This test replaves tensor and boundary conditions by continuous
* functions. This is a prototype forheat conduction solvers.
* **************************************************************** */
void RunTest(std::string op_list_name) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Singular-perturbed nonlinear Laplace Beltrami solver" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_laplace_beltrami.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an MSTK mesh framework
  ParameterList region_list = plist.sublist("Regions Closed");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh = meshfactory.create("test/sphere.exo");

  // extract surface mesh
  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top surface"));

  RCP<Mesh> surfmesh = meshfactory.create(mesh, setnames, AmanziMesh::FACE);

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = surfmesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
  }

  // create boundary data (no mixed bc)
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(surfmesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  bc->bc_model();  // allocate
  bc->bc_value();  // memory

  // create solution map.
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(surfmesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(*cvs));
  solution->PutScalar(0.0);  // solution at time T=0

  CompositeVector phi(*cvs);
  phi.PutScalar(0.2);

  // create source and add it to the operator
  CompositeVector source(*cvs);
  source.PutScalarMasterAndGhosted(0.0);
  
  Epetra_MultiVector& src = *source.ViewComponent("cell");
  for (int c = 0; c < 20; c++) {
    if (MyPID == 0) src[0][c] = 1.0;
  }

  // Create nonlinear coefficient.
  Teuchos::RCP<HeatConduction> knc = Teuchos::rcp(new HeatConduction(surfmesh));

  // MAIN LOOP
  double dT = 1.0;
  for (int loop = 0; loop < 3; loop++) {
    // create diffusion operator
    solution->ScatterMasterToGhosted();
    knc->UpdateValues(*solution);

    Teuchos::ParameterList olist = plist.sublist("PK operator").sublist(op_list_name);
    PDE_DiffusionMFD op(olist, surfmesh);
    op.Init(olist);
    op.SetBCs(bc, bc);

    // get the global operator
    Teuchos::RCP<Operator> global_op = op.global_operator();
    global_op->Init();

    // populate diffusion operator
    op.Setup(K, knc->values(), knc->derivatives());
    op.UpdateMatrices(flux.ptr(), Teuchos::null);

    // add accumulation terms
    PDE_Accumulation op_acc(AmanziMesh::CELL, global_op);
    op_acc.AddAccumulationDelta(*solution, phi, phi, dT, "cell");

    // apply BCs and assemble
    global_op->UpdateRHS(source, false);
    op.ApplyBCs(true, true, true);
    
    // Test SPD properties of the matrix and preconditioner.
    VerificationCV ver(global_op);
    if (loop == 2) {
      global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"));
      global_op->InitializeInverse();
      global_op->ComputeInverse();
      ver.CheckMatrixSPD(true, true);
      ver.CheckPreconditionerSPD(1e-12, true, true);
    }

    // create solver (GMRES with Hypre preconditioner)
    CompositeVector rhs = *global_op->rhs();

    global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"), "Amanzi GMRES", plist.sublist("solvers"));
    global_op->InitializeInverse();
    global_op->ComputeInverse();
    global_op->ApplyInverse(rhs, *solution);

    if (op_list_name == "diffusion operator")
         ver.CheckResidual(*solution, 1.0e-12);

    int num_itrs = global_op->num_itrs();
    CHECK(num_itrs > 5 && num_itrs < 15);

    if (MyPID == 0) {
      double a;
      rhs.Norm2(&a);
      std::cout << "pressure solver (gmres): ||r||=" << global_op->residual() << " itr=" << num_itrs
                << "  ||f||=" << a 
                << " code=" << global_op->returned_code() << std::endl;
    }

    // derive diffusion flux.
    op.UpdateFlux(solution.ptr(), flux.ptr());

    // turn off the source
    source.PutScalar(0.0);
  }
 
  if (MyPID == 0) {
    // visualization
    const Epetra_MultiVector& p = *solution->ViewComponent("cell");
    GMV::open_data_file(*surfmesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::close_data_file();
  }
}

}  // end anonymous namespace


TEST(NONLINEAR_HEAT_CONDUCTION) {
  RunTest("diffusion operator");
}


TEST(NONLINEAR_HEAT_CONDUCTION_SFF) {
  RunTest("diffusion operator Sff");
}
