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

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "GMVMesh.hh"

#include "tensor.hh"
#include "mfd3d_diffusion.hh"

#include "LinearOperatorFactory.hh"
#include "OperatorDefs.hh"
#include "Operator.hh"
#include "OperatorDiffusionSurface.hh"
#include "OperatorSource.hh"

namespace Amanzi{

class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) { 
    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh_);
    cvs.SetGhosted(true);
    cvs.SetComponent("cell", AmanziMesh::CELL, 1);
    cvs.SetOwned(false);
    cvs.AddComponent("face", AmanziMesh::FACE, 1);

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


/* *****************************************************************
* This test replaves tensor and boundary conditions by continuous
* functions. This is a prototype forheat conduction solvers.
* **************************************************************** */
TEST(NONLINEAR_HEAT_CONDUCTION) {
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

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells_owned = surfmesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = surfmesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K.push_back(Kc);
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
    Teuchos::RCP<OperatorSource> op1 = Teuchos::rcp(new OperatorSource(cvs, 0));
    op1->Init();
    op1->UpdateMatrices(source);

    // add accumulation terms
    op1->AddAccumulationTerm(solution, phi, dT, "cell");

    // add diffusion operator
    solution.ScatterMasterToGhosted();
    knc->UpdateValues(solution);

    Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                        .get<Teuchos::ParameterList>("diffusion operator Sff");
    Teuchos::RCP<OperatorDiffusionSurface> op2 = Teuchos::rcp(new OperatorDiffusionSurface(*op1, olist, bc));

    int schema_dofs = op2->schema_dofs();
    CHECK(schema_dofs == Operators::OPERATOR_SCHEMA_DOFS_FACE + Operators::OPERATOR_SCHEMA_DOFS_CELL);
    int schema_prec_dofs = op2->schema_prec_dofs();
    CHECK(schema_prec_dofs == Operators::OPERATOR_SCHEMA_DOFS_FACE);

    op2->InitOperator(K, knc->values(), knc->derivatives(), rho, mu);
    op2->UpdateMatrices(flux, Teuchos::null);
    op2->ApplyBCs();
    op2->SymbolicAssembleMatrix(schema_prec_dofs);
    op2->AssembleMatrix(schema_prec_dofs);

    // create preconditoner
    ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
    op2->InitPreconditioner("Hypre AMG", slist);

    // solve the problem
    ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
    AmanziSolvers::LinearOperatorFactory<OperatorDiffusionSurface, CompositeVector, CompositeVectorSpace> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<OperatorDiffusionSurface, CompositeVector, CompositeVectorSpace> >
       solver = factory.Create("Amanzi GMRES", lop_list, op2);

    CompositeVector rhs = *op2->rhs();
    int ierr = solver->ApplyInverse(rhs, solution);

    if (MyPID == 0) {
      std::cout << "pressure solver (" << solver->name() 
                << "): ||r||=" << solver->residual() << " itr=" << solver->num_itrs()
                << " code=" << solver->returned_code() << std::endl;
    }

    // derive diffusion flux.
    op2->UpdateFlux(solution, *flux);
    // const Epetra_MultiVector& flux_data = *flux->ViewComponent("face");
    // std::cout << flux_data << std::endl;

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


