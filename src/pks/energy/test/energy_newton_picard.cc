#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Epetra_MpiComm.h"

#include "Energy_PK.hh"
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "SolverNewton.hh"
#include "SolverFnBase.hh"
#include "UpwindStandard.hh"

const double TemperatureSource = 1.0; 
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

    // create generic vector
    cvs_ = Teuchos::rcp(new CompositeVectorSpace());
    cvs_->SetMesh(mesh);
    cvs_->SetGhosted(true);
    cvs_->SetComponent("cell", AmanziMesh::CELL, 1);
    cvs_->SetOwned(false);
    cvs_->AddComponent("face", AmanziMesh::FACE, 1);

    solution_ = Teuchos::rcp(new CompositeVector(*cvs_)); 
    solution0_ = Teuchos::rcp(new CompositeVector(*cvs_)); 
    flux_ = Teuchos::rcp(new CompositeVector(*cvs_)); 
    InitialGuess();

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
    dT = 5e-3;
    phi_ = Teuchos::rcp(new CompositeVector(*cvs_));
    phi_->PutScalar(1.0);

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

    // Create upwind model
    Teuchos::ParameterList& ulist = plist.sublist("PK operator").sublist("upwind");
    Teuchos::RCP<Problem> problem = Teuchos::rcp(new Problem());
    upwind_ = Teuchos::rcp(new Operators::UpwindStandard<Problem>(mesh_, problem));
    upwind_->Init(ulist);

    // Update values
    UpdateValues(*solution_);

    // create diffusion operator
    Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("diffusion operator");
    op_ = Teuchos::rcp(new Operators::OperatorDiffusion(cvs_, olist, bc));
    op_->Init();

    int schema_dofs = op_->schema_dofs();
    int schema_prec_dofs = op_->schema_prec_dofs();

    op_->Setup(K, k, dkdT, rho, mu);
    op_->UpdateMatrices(flux_, solution_);
    op_->AddAccumulationTerm(*solution0_, *phi_, dT, "cell");
    op_->ApplyBCs();
    op_->SymbolicAssembleMatrix(schema_prec_dofs);
    op_->AssembleMatrix(schema_prec_dofs);

    // create preconditoner
    Teuchos::ParameterList slist = plist.sublist("Preconditioners");
    op_->InitPreconditioner("Hypre AMG", slist);
  }

  double Conduction(int c, double T) const {
    ASSERT(T > 0.0);
    return T * T * T;
  }

  double ConductionDerivative(int c, double T) const {
    ASSERT(T > 0.0);
    return 3 * T * T;
  }

  void InitialGuess() { 
    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    Epetra_MultiVector& sol_c = *solution_->ViewComponent("cell", true);
  
    for (int c = 0; c < ncells_wghost; ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      double a = TemperatureSource;
      sol_c[0][c] = a / 2 - a / M_PI * atan(20 * (xc[0] - 1.0));
    }

    int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    Epetra_MultiVector& sol_f = *solution_->ViewComponent("face", true);
  
    for (int f = 0; f < nfaces_wghost; ++f) {
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      double a = TemperatureSource;
      sol_f[0][f] = a / 2 - a / M_PI * atan(100 * (xf[0] - 1.0));
    }

    *solution0_ = *solution_;
  }

  void UpdateValues(const CompositeVector& u) { 
    const Epetra_MultiVector& uc = *u.ViewComponent("cell", true); 
    const Epetra_MultiVector& kc = *k->ViewComponent("cell", true); 
    const Epetra_MultiVector& dkdT_c = *dkdT->ViewComponent("cell", true); 

    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    for (int c = 0; c < ncells_wghost; c++) {
      double u = uc[0][c];
      kc[0][c] = Conduction(c, u);
      dkdT_c[0][c] = ConductionDerivative(c, u);
    }

    upwind_->Compute(*flux_, bc_model, bc_value, *k, *k, &Problem::Conduction);
    upwind_->Compute(*flux_, bc_model, bc_value, *dkdT, *dkdT, &Problem::ConductionDerivative);
  }

  // virtual member functions
  void Residual(const Teuchos::RCP<Amanzi::CompositeVector>& u,
                const Teuchos::RCP<Amanzi::CompositeVector>& f) {
    op_->Init();
    UpdateValues(*u);
    op_->UpdateMatrices(Teuchos::null, u);
    op_->AddAccumulationTerm(*solution0_, *phi_, dT, "cell");
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
    UpdateValues(*up);
    op_->UpdateMatrices(flux_, up);
    op_->AddAccumulationTerm(*solution0_, *phi_, dT, "cell");
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
  Teuchos::RCP<Operators::UpwindStandard<Problem> > upwind_;

  double dT;
  Teuchos::RCP<CompositeVector> phi_;

  Teuchos::RCP<CompositeVector> solution_;
  Teuchos::RCP<CompositeVector> solution0_;
  Teuchos::RCP<CompositeVector> flux_;
};

typedef double(Problem::*ModelUpwindFn)(int c, double T) const; 

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
  Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 3.0, 1.0, 60, 10, gm);

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
  Epetra_MultiVector p0(*problem->solution()->ViewComponent("cell"));

  // solve
  newton_picard->Solve(problem->solution());
  Epetra_MultiVector& p1 = *problem->solution()->ViewComponent("cell");

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"energy.gmv");
    GMV::start_data();
    GMV::write_cell_data(p0, 0, "p0");
    GMV::write_cell_data(p1, 0, "p1");
    GMV::close_data_file();
  }
};


