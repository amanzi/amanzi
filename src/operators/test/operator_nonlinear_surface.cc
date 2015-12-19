/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
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
#include "LinearOperatorFactory.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "mfd3d_diffusion.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "Operator.hh"
#include "OperatorAccumulation.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionMFD.hh"
#include "Verification.hh"

namespace Amanzi{

class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) { 
    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("face", AmanziMesh::FACE, 1);

    values_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
    derivatives_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs, true));
  }
  ~HeatConduction() {};

  // main members
  void UpdateValues(const CompositeVector& u) { 
    const Epetra_MultiVector& uc = *u.ViewComponent("cell", true); 
    const Epetra_MultiVector& values_c = *values_->ViewComponent("cell", true); 

    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    for (int c = 0; c < ncells; c++) {
      values_c[0][c] = 0.3 + uc[0][c];
    }

    const Epetra_MultiVector& uf = *u.ViewComponent("face", true); 
    const Epetra_MultiVector& values_f = *values_->ViewComponent("face", true); 
    int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
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

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Singular-perturbed nonlinear Laplace Beltrami solver" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_laplace_beltrami.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an MSTK mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions Closed");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory("test/sphere.exo", gm);
  RCP<const Mesh_MSTK> mesh_mstk = rcp_static_cast<const Mesh_MSTK>(mesh);

  // extract surface mesh
  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top surface"));

  RCP<Mesh> surfmesh = Teuchos::rcp(new Mesh_MSTK(*mesh_mstk, setnames, AmanziMesh::FACE));

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the nonlinear coefficient.
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = surfmesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = surfmesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
  }
  double rho(1.0), mu(1.0);

  // create boundary data (no mixed bc)
  std::vector<int> bc_model(nfaces_wghost, OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);
  std::vector<double> bc_mixed;
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create solution map.
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(surfmesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(*cvs));

  CompositeVector solution(*cvs);
  solution.PutScalar(0.0);  // solution at time T=0

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
    solution.ScatterMasterToGhosted();
    knc->UpdateValues(solution);

    Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                        .get<Teuchos::ParameterList>(op_list_name);
    OperatorDiffusionMFD op(olist, surfmesh);
    op.SetBCs(bc, bc);

    // get the global operator
    Teuchos::RCP<Operator> global_op = op.global_operator();
    global_op->Init();

    // populate diffusion operator
    op.Setup(K, knc->values(), knc->derivatives());
    op.UpdateMatrices(flux.ptr(), Teuchos::null);

    // add accumulation terms
    OperatorAccumulation op_acc(AmanziMesh::CELL, global_op);
    op_acc.AddAccumulationTerm(solution, phi, dT, "cell");

    // apply BCs and assemble
    global_op->UpdateRHS(source, false);
    op.ApplyBCs(true, true);
    global_op->SymbolicAssembleMatrix();
    global_op->AssembleMatrix();
    
    // create preconditoner
    ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
    global_op->InitPreconditioner("Hypre AMG", slist);

    // Test SPD properties of the matrix and preconditioner.
    if (loop == 2) {
      Verification ver(global_op);
      ver.CheckMatrixSPD(true, true);
      ver.CheckPreconditionerSPD(true, true);
    }

    // solve the problem
    ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
    AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
       solver = factory.Create("Amanzi GMRES", lop_list, global_op);

    CompositeVector rhs = *global_op->rhs();
    int ierr = solver->ApplyInverse(rhs, solution);

    int num_itrs = solver->num_itrs();
    CHECK(num_itrs > 5 && num_itrs < 15);

    if (MyPID == 0) {
      double a;
      rhs.Norm2(&a);
      std::cout << "pressure solver (" << solver->name() 
                << "): ||r||=" << solver->residual() << " itr=" << num_itrs
                << "  ||f||=" << a 
                << " code=" << solver->returned_code() << std::endl;
    }

    // derive diffusion flux.
    op.UpdateFlux(solution, *flux);

    // turn off the source
    source.PutScalar(0.0);
  }
 
  if (MyPID == 0) {
    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("cell");
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
