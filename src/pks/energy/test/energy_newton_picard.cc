/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This is the energy component of the Amanzi code.

*/

#include <iostream>
#include "UnitTest++.h"

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "AmanziComm.hh"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "Mesh_Algorithms.hh"
#include "Operator.hh"
#include "PDE_Accumulation.hh"
#include "PDE_DiffusionMFD.hh"
#include "SolverFactory.hh"
#include "SolverFnBase.hh"
#include "UpwindFlux.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Energy
#include "Energy_PK.hh"

const double TemperatureSource = 1.0;
const double TemperatureFloor = 0.02;
const std::string SOLVERS[3] = { "NKA", "Newton-Picard", "JFNK" };

namespace Amanzi {

class HeatConduction : public AmanziSolvers::SolverFnBase<CompositeVector> {
 public:
  HeatConduction() { op_name_ = "diffusion operator Newton-Picard"; }
  HeatConduction(std::string& op_name) { op_name_ = "diffusion operator " + op_name; }
  ~HeatConduction(){};

  // EOS
  double Conduction(int c, double T) const;
  double ConductionDerivative(int c, double T) const;

  // Initialization routines
  void Init(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::ParameterList& plist);
  void InitialGuess();

  // Virtual member functions for base class Solv erFnBase
  void Residual(const Teuchos::RCP<CompositeVector>& u, const Teuchos::RCP<CompositeVector>& f);

  void UpdatePreconditioner(const Teuchos::RCP<const Amanzi::CompositeVector>& up);

  int ApplyPreconditioner(const Teuchos::RCP<const CompositeVector>& u,
                          const Teuchos::RCP<CompositeVector>& hu);

  double ErrorNorm(const Teuchos::RCP<const CompositeVector>& u,
                   const Teuchos::RCP<const CompositeVector>& du);

  void ChangedSolution(){};

  // supporting members
  void UpdateValues(const CompositeVector& u);

  // access
  CompositeVectorSpace& cvs() { return *cvs_; }
  Teuchos::RCP<CompositeVector> solution() { return solution_; }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;
  std::string op_name_;

  Teuchos::RCP<CompositeVectorSpace> cvs_;
  Teuchos::RCP<Operators::Operator> op_;
  Teuchos::RCP<Operators::PDE_DiffusionMFD> op_diff_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc_;

  Teuchos::RCP<Operators::BCs> bc_;

  std::vector<WhetStone::Tensor> K;
  Teuchos::RCP<CompositeVector> k, dkdT;
  Teuchos::RCP<Operators::UpwindFlux> upwind_;

  double dT;
  Teuchos::RCP<CompositeVector> phi_;

  Teuchos::RCP<CompositeVector> solution_, solution0_;
  Teuchos::RCP<CompositeVector> flux_;
};


/* ******************************************************************
* Equation of states.
****************************************************************** */
double
HeatConduction::Conduction(int c, double T) const
{
  AMANZI_ASSERT(T > 0.0);
  return T * T * T;
}


double
HeatConduction::ConductionDerivative(int c, double T) const
{
  AMANZI_ASSERT(T > 0.0);
  return 3 * T * T;
}


/* ******************************************************************
* Initialization requires global parameter list.
****************************************************************** */
void
HeatConduction::Init(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::ParameterList& plist)
{
  mesh_ = mesh;
  plist_ = plist;

  // create generic vector
  cvs_ = Teuchos::rcp(new CompositeVectorSpace());
  cvs_->SetMesh(mesh);
  cvs_->SetGhosted(true);
  cvs_->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  cvs_->SetOwned(false);
  cvs_->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  // solutions at T=T0 and T=T0+dT
  solution_ = Teuchos::rcp(new CompositeVector(*cvs_));
  solution0_ = Teuchos::rcp(new CompositeVector(*cvs_));
  flux_ = Teuchos::rcp(new CompositeVector(*cvs_));

  InitialGuess();

  // create BCs
  int ncells_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  bc_ = Teuchos::rcp(
    new Operators::BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc_->bc_model();
  std::vector<double>& bc_value = bc_->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);

    if (fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = 0.0;
    } else if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = TemperatureSource;
    } else if (fabs(xf[0] - 3.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = TemperatureFloor;
    }
  }

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

  // create temperature-dependent data
  k = Teuchos::rcp(new CompositeVector(*cvs_));
  dkdT = Teuchos::rcp(new CompositeVector(*cvs_));

  // Create upwind model
  Teuchos::ParameterList& ulist = plist.sublist("PK operator").sublist("upwind");
  Teuchos::RCP<HeatConduction> problem = Teuchos::rcp(new HeatConduction());
  upwind_ = Teuchos::rcp(new Operators::UpwindFlux(mesh_));
  upwind_->Init(ulist);

  // Update conductivity values
  UpdateValues(*solution_);

  // create the operators
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist(op_name_);
  op_diff_ = Teuchos::rcp(new Operators::PDE_DiffusionMFD(olist, mesh_));
  op_diff_->Init(olist);
  op_diff_->SetBCs(bc_, bc_);

  op_ = op_diff_->global_operator();
  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op_));
  op_->Init();

  // set up the local matrices
  Teuchos::RCP<std::vector<WhetStone::Tensor>> Kptr = Teuchos::rcpFromRef(K);
  op_diff_->Setup(Kptr, k, dkdT);
  op_diff_->UpdateMatrices(flux_.ptr(), solution_.ptr());
  op_diff_->UpdateMatricesNewtonCorrection(flux_.ptr(), solution_.ptr(), 1.0);
  op_acc_->AddAccumulationDelta(*solution0_, *phi_, *phi_, dT, "cell");

  // form the global matrix
  op_diff_->ApplyBCs(true, true, true);
  op_->SymbolicAssembleMatrix();
  op_->AssembleMatrix();

  // create preconditoner
  Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  op_->set_inverse_parameters(slist);
  op_->InitializeInverse();
  op_->ComputeInverse();
}


/* ******************************************************************
* Initialization for solution at time T=T0.
****************************************************************** */
void
HeatConduction::InitialGuess()
{
  int ncells_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  Epetra_MultiVector& sol_c = *solution_->ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost; ++c) {
    const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
    double a = TemperatureSource;
    sol_c[0][c] = a / 2 - a / M_PI * atan(20 * (xc[0] - 1.0));
  }

  int nfaces_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  Epetra_MultiVector& sol_f = *solution_->ViewComponent("face", true);

  for (int f = 0; f < nfaces_wghost; ++f) {
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    double a = TemperatureSource;
    sol_f[0][f] = a / 2 - a / M_PI * atan(100 * (xf[0] - 1.0));
  }

  *solution0_ = *solution_;
}


/* ******************************************************************
* Residual re-calculates the matrix.
****************************************************************** */
void
HeatConduction::Residual(const Teuchos::RCP<CompositeVector>& u,
                         const Teuchos::RCP<CompositeVector>& f)
{
  op_->Init();
  UpdateValues(*u);
  op_diff_->UpdateMatrices(Teuchos::null, u.ptr());
  op_acc_->AddAccumulationDelta(*solution0_, *phi_, *phi_, dT, "cell");
  op_diff_->ApplyBCs(true, true, true);
  op_->ComputeNegativeResidual(*u, *f);
}


/* ******************************************************************
* Precondtioner re-calculates the matrix.
****************************************************************** */
void
HeatConduction::UpdatePreconditioner(const Teuchos::RCP<const CompositeVector>& up)
{
  op_diff_->UpdateFlux(up.ptr(), flux_.ptr());

  // Calculate new matrix.
  UpdateValues(*up);
  op_->Init();
  op_diff_->UpdateMatrices(flux_.ptr(), up.ptr());
  op_diff_->UpdateMatricesNewtonCorrection(flux_.ptr(), up.ptr(), 1.0);
  op_acc_->AddAccumulationDelta(*solution0_, *phi_, *phi_, dT, "cell");
  op_diff_->ApplyBCs(true, true, true);

  // Assemble matrix and calculate preconditioner.
  op_->AssembleMatrix();

  Teuchos::ParameterList pc_list = plist_.sublist("preconditioners").sublist("Hypre AMG");
  op_->set_inverse_parameters(pc_list);
  op_->InitializeInverse();
  op_->ComputeInverse();
}


int
HeatConduction::ApplyPreconditioner(const Teuchos::RCP<const CompositeVector>& u,
                                    const Teuchos::RCP<CompositeVector>& hu)
{
  int ierr = op_->ApplyInverse(*u, *hu);
  return (ierr > 0) ? 0 : 1;
  //return op_->ApplyInverse(*u, *hu);
}


/* ******************************************************************
* Definition of error in the nonlinear solver.
****************************************************************** */
double
HeatConduction::ErrorNorm(const Teuchos::RCP<const CompositeVector>& u,
                          const Teuchos::RCP<const CompositeVector>& du)
{
  double norm_du;
  du->NormInf(&norm_du);
  return norm_du;
}


/* ******************************************************************
* Recalculates conductivity in cells and on faces.
****************************************************************** */
void
HeatConduction::UpdateValues(const CompositeVector& u)
{
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);
  Epetra_MultiVector& kc = *k->ViewComponent("cell", true);
  Epetra_MultiVector& dkdT_c = *dkdT->ViewComponent("cell", true);

  int ncells_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  for (int c = 0; c < ncells_wghost; c++) {
    double temp = uc[0][c];
    kc[0][c] = Conduction(c, temp);
    dkdT_c[0][c] = ConductionDerivative(c, temp);
  }

  std::vector<int>& bc_model = bc_->bc_model();
  std::vector<double>& bc_value = bc_->bc_value();

  // add boundary face component
  Epetra_MultiVector& k_f = *k->ViewComponent("face", true);
  Epetra_MultiVector& dkdT_f = *dkdT->ViewComponent("face", true);

  for (int f = 0; f != bc_model.size(); ++f) {
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);
      k_f[0][f] = Conduction(c, bc_value[f]);
      dkdT_f[0][f] = ConductionDerivative(c, bc_value[f]);
    }
  }

  upwind_->Compute(*flux_, bc_model, *k);
  upwind_->Compute(*flux_, bc_model, *dkdT);
}

} // namespace Amanzi


/* ******************************************************************
* Comparison of various nonlinear solvers.
****************************************************************** */
TEST(NKA)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  using namespace Amanzi::Energy;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  // read parameter list
  std::string xmlFileName = "test/energy_newton_picard.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 3.0, 1.0, 60, 10);

  // create problem
  for (int i = 0; i < 3; ++i) {
    if (MyPID == 0) std::cout << "Test: nonlinear diffusion: " << SOLVERS[i] << std::endl;
    std::string solver_name(SOLVERS[i]);
    Teuchos::RCP<HeatConduction> problem = Teuchos::rcp(new HeatConduction(solver_name));
    problem->Init(mesh, plist);

    // create the Solver
    AmanziSolvers::SolverFactory<CompositeVector, CompositeVectorSpace> factory;
    Teuchos::RCP<AmanziSolvers::Solver<CompositeVector, CompositeVectorSpace>> solver =
      factory.Create(SOLVERS[i], plist);
    solver->Init(problem, problem->cvs());

    // solve
    problem->InitialGuess();
    solver->Solve(problem->solution());

    // checks
    CHECK(solver->num_itrs() < 9);
    if (i == 1) CHECK_CLOSE(0., solver->residual(), 1e-8);

    Epetra_MultiVector p0(*problem->solution()->ViewComponent("cell"));
    Epetra_MultiVector& p1 = *problem->solution()->ViewComponent("cell");
    if (MyPID == 0) {
      GMV::open_data_file(*mesh, (std::string) "energy.gmv");
      GMV::start_data();
      GMV::write_cell_data(p0, 0, "p0");
      GMV::write_cell_data(p1, 0, "p1");
      GMV::close_data_file();
    }
  }
}
