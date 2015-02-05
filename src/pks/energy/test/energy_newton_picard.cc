#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Epetra_MpiComm.h"

#include "Energy_PK.hh"
#include "MeshFactory.hh"
#include "SolverNewton.hh"
#include "SolverFnBase.hh"

const double TemperatureSource = 100.0; 
const double TemperatureFloor = 0.02; 


namespace Amanzi {

class Problem : public Amanzi::AmanziSolvers::SolverFnBase<Amanzi::CompositeVector> {
 public:
  Problem() {};
  ~Problem() {};

  // Initialization requires global parameter list.
  void Init(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::ParameterList& plist) {
    mesh_ = mesh;
    plist_ = plist;

    // create vector
    cvs_ = Teuchos::rcp(new CompositeVectorSpace());
    cvs_->SetMesh(mesh);
    cvs_->SetGhosted(true);
    cvs_->SetComponent("cell", AmanziMesh::CELL, 1);
    cvs_->SetOwned(false);
    cvs_->AddComponent("face", AmanziMesh::FACE, 1);

    solution_ = Teuchos::rcp(new CompositeVector(*cvs_)); 
    flux_ = Teuchos::rcp(new CompositeVector(*cvs_)); 

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
    dT = 0.02;
    phi_ = Teuchos::rcp(new CompositeVector(*cvs_));
    phi_->PutScalar(0.2);

    // create static diffusion data
    for (int c = 0; c < ncells_owned; c++) {
      WhetStone::Tensor Kc(2, 1);
      Kc(0, 0) = 1.0;
      K.push_back(Kc);
    }
    double rho(1.0), mu(1.0);

    // create temperature-dependent data
    k = Teuchos::rcp(new CompositeVector(*cvs_));
    dkdT = Teuchos::rcp(new CompositeVector(*cvs_));
    UpdateValues(*solution_);

    // create diffusion operator
    Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("diffusion operator");
    op_ = Teuchos::rcp(new Operators::OperatorDiffusion(cvs_, olist, bc));
    op_->Init();

    int schema_dofs = op_->schema_dofs();
    int schema_prec_dofs = op_->schema_prec_dofs();

    op_->Setup(K, k, dkdT, rho, mu);
    op_->UpdateMatrices(flux_, solution_);
    op_->AddAccumulationTerm(*solution_, *phi_, dT, "cell");
    op_->ApplyBCs();
    op_->SymbolicAssembleMatrix(schema_prec_dofs);
    op_->AssembleMatrix(schema_prec_dofs);

    // create preconditoner
    Teuchos::ParameterList slist = plist.sublist("Preconditioners");
    op_->InitPreconditioner("Hypre AMG", slist);
  }

  void InitialGuess() { solution_->PutScalar(1.0); }

  void UpdateValues(const CompositeVector& u) { 
    const Epetra_MultiVector& uc = *u.ViewComponent("cell", true); 
    const Epetra_MultiVector& kc = *k->ViewComponent("cell", true); 

    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    for (int c = 0; c < ncells; c++) {
      kc[0][c] = std::pow(uc[0][c], 3.0);
    }

    k->PutScalar(1.0);
    dkdT->PutScalar(0.0);
  }

  // virtual member functions
  void Residual(const Teuchos::RCP<Amanzi::CompositeVector>& u,
                const Teuchos::RCP<Amanzi::CompositeVector>& f) {
    op_->Init();
    op_->UpdateMatrices(flux_, u);
    op_->AddAccumulationTerm(*u, *phi_, dT, "cell");
    op_->ApplyBCs();
    op_->ComputeNegativeResidual(*u, *f);
  }

  void ApplyPreconditioner(const Teuchos::RCP<const Amanzi::CompositeVector>& u,
                           const Teuchos::RCP<Amanzi::CompositeVector>& hu) {
    op_->ApplyInverse(*u, *hu);
  }

  double ErrorNorm(const Teuchos::RCP<const Amanzi::CompositeVector>& u,
                   const Teuchos::RCP<const Amanzi::CompositeVector>& du) {
    double norm_du, norm_u;
    du->NormInf(&norm_du);
    return norm_du;
  }

  void UpdatePreconditioner(const Teuchos::RCP<const Amanzi::CompositeVector>& up) {
    op_->UpdateFlux(*up, *flux_);

    // Calculate new matrix.
    op_->Init();
    op_->UpdateMatrices(flux_, up);
    op_->AddAccumulationTerm(*up, *phi_, dT, "cell");
    op_->ApplyBCs();

    // Assemble matrix and calculate preconditioner.
    int schema_prec_dofs = op_->schema_prec_dofs();
    op_->AssembleMatrix(schema_prec_dofs);

    Teuchos::ParameterList prec_list = plist_.sublist("Preconditioners");
    op_->InitPreconditioner("Hypre AMG", prec_list);
  }

  void ChangedSolution() {};

  // access
  CompositeVectorSpace& cvs() { return *cvs_; }
  Teuchos::RCP<CompositeVector> solution() { return solution_; }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;

  Teuchos::RCP<CompositeVectorSpace> cvs_;
  Teuchos::RCP<Operators::OperatorDiffusion> op_;
  std::vector<int> bc_model;
  std::vector<double> bc_value, bc_mixed;

  std::vector<WhetStone::Tensor> K;
  Teuchos::RCP<CompositeVector> k, dkdT;

  double dT;
  Teuchos::RCP<CompositeVector> phi_;

  Teuchos::RCP<CompositeVector> solution_;
  Teuchos::RCP<CompositeVector> flux_;
};

}  // namespace Amanzi


/* ******************************************************************
* Test 1.
****************************************************************** */
TEST(NEWTON_PICARD) {
  using namespace Amanzi;
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
  Teuchos::RCP<Problem> problem = Teuchos::rcp(new Problem());
  problem->Init(mesh, plist);

  // create the Solver
  Teuchos::ParameterList& prec_list = plist.get<Teuchos::ParameterList>("Preconditioners");
  Teuchos::ParameterList slist = plist.sublist("Newton parameters");

  Teuchos::RCP<AmanziSolvers::SolverNewton<CompositeVector, CompositeVectorSpace> > 
      newton_picard = Teuchos::rcp(new AmanziSolvers::SolverNewton<CompositeVector, CompositeVectorSpace>(slist));
  newton_picard->Init(problem, problem->cvs());

  // initial guess
  problem->InitialGuess();

  // solve
  newton_picard->Solve(problem->solution());
};


