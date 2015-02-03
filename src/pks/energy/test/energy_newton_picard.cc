#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Epetra_MpiComm.h"

#include "Energy_PK.hh"
#include "MeshFactory.hh"
#include "SolverFnBase01.hh"
#include "SolverNewton.hh"

const double TemperatureSource = 100.0; 
const double TemperatureFloor = 0.02; 

SUITE(SOLVERS) {
using namespace Amanzi;

class Problem {
 public:
  Problem() {};
  ~Problem() {};

  void Setup(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::ParameterList& plist) {
    mesh_ = mesh;
    // create vector
    cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh);
    cvs->SetGhosted(true);
    cvs->SetComponent("cell", AmanziMesh::CELL, 1);
    cvs->SetOwned(false);
    cvs->AddComponent("face", AmanziMesh::FACE, 1);

    solution = Teuchos::rcp(new CompositeVector(*cvs)); 
    flux = Teuchos::rcp(new CompositeVector(*cvs)); 

    // create BCs
    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

    bc_model.assign(nfaces_wghost, Operators::OPERATOR_BC_NONE);
    bc_value.assign(nfaces_wghost, 0.0);

    for (int f = 0; f < nfaces_wghost; f++) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);

      if (fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = 0.0;
      } else if(fabs(xf[0]) < 1e-6) {
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = TemperatureSource;
      } else if(fabs(xf[0] - 3.0) < 1e-6) {
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = TemperatureFloor;
      }
    }
    Teuchos::RCP<Operators::BCs> bc =
        Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

    // constant accumulation term
    double dT = 0.02;
    CompositeVector phi(*solution);
    phi.PutScalar(0.2);

    // create statis diffusion data
    for (int c = 0; c < ncells_owned; c++) {
      WhetStone::Tensor Kc(2, 1);
      Kc(0, 0) = 1.0;
      K.push_back(Kc);
    }
    double rho(1.0), mu(1.0);

    // create temperature-dependent data
    k = Teuchos::rcp(new CompositeVector(*cvs));
    dkdT = Teuchos::rcp(new CompositeVector(*cvs));
    UpdateValues(*solution);

    // create diffusion operator
    Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("diffusion operator");
    op = Teuchos::rcp(new Operators::OperatorDiffusion(cvs, olist, bc));
    op->Init();

    int schema_dofs = op->schema_dofs();
    int schema_prec_dofs = op->schema_prec_dofs();

    op->Setup(K, k, dkdT, rho, mu);
    op->UpdateMatrices(flux, solution);
    op->AddAccumulationTerm(*solution, phi, dT, "cell");
    op->ApplyBCs();
    op->SymbolicAssembleMatrix(schema_prec_dofs);
    op->AssembleMatrix(schema_prec_dofs);

    // create preconditoner using the base operator class
    Teuchos::ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
    op->InitPreconditioner("Hypre AMG", slist);
  }

  void InitialGuess(Amanzi::CompositeVector& u) {
    u.PutScalar(1.0);
  }

  void UpdateValues(const CompositeVector& u) { 
    const Epetra_MultiVector& uc = *u.ViewComponent("cell", true); 
    const Epetra_MultiVector& kc = *k->ViewComponent("cell", true); 

    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    for (int c = 0; c < ncells; c++) {
      kc[0][c] = std::pow(uc[0][c], 3.0);
    }

    k->PutScalar(1.0);
    dkdT->PutScalar(1.0);
  }

 public:
  Teuchos::RCP<CompositeVectorSpace> cvs;
  Teuchos::RCP<Operators::OperatorDiffusion> op;

 public:
  Teuchos::RCP<CompositeVector> solution;
  Teuchos::RCP<CompositeVector> flux;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::vector<int> bc_model;
  std::vector<double> bc_value, bc_mixed;

  std::vector<WhetStone::Tensor> K;
  Teuchos::RCP<CompositeVector> k, dkdT;
};


/* ******************************************************************
* Test 1.
****************************************************************** */
TEST_FIXTURE(Problem, NEWTON_PICARD) {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  using namespace Amanzi::Energy;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "Test: nonlinear diffusion: Newton-Picard." << std::endl;

  // read parameter list 
  std::string xmlFileName = "test/energy_newton_picard.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);
 
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 3.0, 1.0, 30, 10, gm);

  // create problem
  Setup(mesh, plist);

  // create the Solver
  Teuchos::ParameterList& prec_list = plist.get<Teuchos::ParameterList>("Preconditioners");
  Teuchos::RCP<SolverFnBase01> fn = Teuchos::rcp(new SolverFnBase01(op, prec_list));
  Teuchos::ParameterList slist = plist.sublist("Newton parameters");

  Teuchos::RCP<AmanziSolvers::SolverNewton<CompositeVector, CompositeVectorSpace> > 
      newton_picard = Teuchos::rcp(new AmanziSolvers::SolverNewton<CompositeVector, CompositeVectorSpace>(slist));
  newton_picard->Init(fn, *cvs);

  // initial guess
  InitialGuess(*solution);

  // solve
  newton_picard->Solve(solution);
};

}  // SUITE

