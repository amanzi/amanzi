/*
  The transport component of the Amanzi code, serial unit tests.
  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "GMVMesh.hh"

#include "SolverNewton.hh"
#include "SolverFnNLFV.hh"
#include "LinearOperatorFactory.hh"

#include "Transport_State.hh"
#include "Transport_PK.hh"
#include "DispersionMatrixFactory.hh"

using namespace Teuchos;
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziTransport;
using namespace Amanzi::AmanziGeometry;

SUITE(DISPERSION_MATRIX) {
  class Problem {
   public:
    Problem() {};
    ~Problem() { delete comm; }

    void Init() {
      comm = new Epetra_MpiComm(MPI_COMM_WORLD);

      /* read parameter list */
      string xmlFileName = "test/transport_matrix.xml";

      ParameterXMLFileReader xmlreader(xmlFileName);
      plist = xmlreader.getParameters();

      /* create an MSTK mesh framework */
      ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
      GeometricModelPtr gm = new GeometricModel(2, region_list, (Epetra_MpiComm *)comm);

      FrameworkPreference pref;
      pref.clear();
      pref.push_back(MSTK);

      MeshFactory meshfactory(comm);
      meshfactory.preference(pref);
      mesh = meshfactory(0.0,0.0, 1.0,1.0, 10,10, gm);

      /* update vebose object */
      Amanzi::VerboseObject::hide_line_prefix = true;

      /* create a transport state from the MPC state and populate it */
      TS = rcp(new Transport_State(mesh, 1));
      TS->Initialize();

      AmanziGeometry::Point velocity(3.0, 4.0);
      TS->set_darcy_flux(velocity);
      TS->set_porosity(0.5);
      TS->set_water_saturation(1.0);
      TS->set_prev_water_saturation(1.0);
      TS->set_water_density(1000.0);
      TS->set_total_component_concentration(0.0, 0);

      /* initialize a transport process kernel from a transport state */
      TPK = rcp(new Transport_PK(plist, TS));
      TPK->InitPK();
    }
    
    void InitSOL(Epetra_Vector& u) {
      double x, y, tmp;
      for (int c = 0; c < u.MyLength(); c++) {
        const AmanziGeometry::Point& xm = mesh->cell_centroid(c);
        x = xm[0];  
        y = xm[1];  

        tmp = x * y * (1.0 - x) * (1.0 - y);
        u[c] = tmp * tmp;
      }
    }

    void InitRHS(Epetra_Vector& b, double r) {
      double x, y, x1, y1, tmp, a;
      for (int c = 0; c < b.MyLength(); c++) {
        const AmanziGeometry::Point& xm = mesh->cell_centroid(c);
        x = xm[0];  
        y = xm[1];  

        x1 = 1.0 - x;
        y1 = 1.0 - y;
        tmp = -28 * y * y * y1 * y1 * (1 - 6 * x + 6 * x * x)
              -96 * x * x1 * (1 - 2 * x) * y * y1 * (1 - 2 * y)
              -42 * x * x * x1 * x1 * (1 - 6 * y + 6 * y * y);

        a = x * y * x1 * y1;
        b[c] = (tmp + r * a * a) * mesh->cell_volume(c);
      }
    }

   public:
    Epetra_MpiComm* comm;
    RCP<Mesh> mesh;
    ParameterList plist;

    RCP<Transport_State> TS;
    RCP<Transport_PK> TPK;
  };

  /* **************************************************************** */
  TEST_FIXTURE(Problem, MATRIX_TPFA) {
    cout << endl
         << "Test: bi-quartic solution with TPFA (no convergence)" << endl;
    Init();

    /* generate a dispersion matrix */
    const Epetra_Vector& phi = TS->ref_porosity();
    const Epetra_Vector& flux = TS->ref_darcy_flux();
    const Epetra_Vector& ws = TS->ref_water_saturation();

    Teuchos::RCP<Dispersion_TPFA> matrix = Teuchos::rcp(new Dispersion_TPFA(&(TPK->dispersion_models()), mesh, TS));
    matrix->Init();
    matrix->InitPreconditioner("Hypre AMG", plist.sublist("Preconditioners"));
    matrix->SymbolicAssembleGlobalMatrix();
    matrix->CalculateDispersionTensor(flux, phi, ws);
    matrix->AssembleGlobalMatrixTPFA(TS);

    /* populate right-hand side and solution */
    Epetra_Vector u(phi), r(phi), b(phi);
    InitSOL(u);
    InitRHS(b, 0.0);

    matrix->Apply(u, r);
    r.Update(-1.0, b, 1.0);

    double bnorm, residual; 
    b.Norm2(&bnorm);
    r.Norm2(&residual);
    residual /= bnorm;
    cout << "Relative residual: " << residual << endl;
    CHECK(residual < 0.5);
  }

  /* **************************************************************** */
  TEST_FIXTURE(Problem, MATRIX_NLFV) {
    cout << endl
         << "Test: bi-quartic solution with NLFV (spartial convergence)" << endl;
    Init();

    /* generate a dispersion matrix */
    const Epetra_Vector& phi = TS->ref_porosity();
    const Epetra_Vector& flux = TS->ref_darcy_flux();
    const Epetra_Vector& ws = TS->ref_water_saturation();

    /* create matrix */
    Teuchos::RCP<Dispersion_TPFA> matrix = Teuchos::rcp(new Dispersion_TPFA(&(TPK->dispersion_models()), mesh, TS));
    matrix->Init();
    matrix->InitPreconditioner("Hypre AMG", plist.sublist("Preconditioners"));

    matrix->SymbolicAssembleGlobalMatrix();
    matrix->CalculateDispersionTensor(flux, phi, ws);
    matrix->InitNLFV();
    matrix->CreateFluxStencils();

    /* populate solution and right-hand side */
    Epetra_Vector u(phi), b(phi);
    InitSOL(u);
    InitRHS(b, 1.0);

    /* update matrix */
    matrix->AssembleGlobalMatrixNLFV(u);
    matrix->AddTimeDerivative(1.0, phi, ws);

    /* Compute residual */
    Epetra_Vector r(phi);
    matrix->Apply(u, r);
    r.Update(-1.0, b, 1.0);

    double bnorm, residual; 
    b.Norm2(&bnorm);
    r.Norm2(&residual);
    residual /= bnorm;
    cout << "Relative residual: " << residual << endl;
    CHECK(residual < 0.2);
  }

  /* **************************************************************** */
  TEST_FIXTURE(Problem, MATRIX_FACTORY) {
    cout << endl
         << "Test: bi-quartic solution with TPFA (using factory)" << endl;
    Init();

    /* generate a dispersion matrix */
    const Epetra_Vector& phi = TS->ref_porosity();
    const Epetra_Vector& flux = TS->ref_darcy_flux();
    const Epetra_Vector& ws = TS->ref_water_saturation();

    DispersionMatrixFactory factory;
    Teuchos::RCP<Dispersion> matrix = factory.Create("tpfa", &(TPK->dispersion_models()), mesh, TS);
    matrix->Init();
    matrix->InitPreconditioner("Hypre AMG", plist.sublist("Preconditioners"));
    matrix->SymbolicAssembleGlobalMatrix();
    matrix->CalculateDispersionTensor(flux, phi, ws);
    matrix->AssembleGlobalMatrix(phi);

    /* populate right-hand side and solution */
    Epetra_Vector u(phi), r(phi), b(phi);
    InitSOL(u);
    InitRHS(b, 0.0);

    /* calcualte residual */
    matrix->Apply(u, r);
    r.Update(-1.0, b, 1.0);

    double bnorm, residual; 
    b.Norm2(&bnorm);
    r.Norm2(&residual);
    residual /= bnorm;
    cout << "Relative residual: " << residual << endl;
    CHECK(residual < 0.5);
  }

  /* **************************************************************** */
  TEST_FIXTURE(Problem, MATRIX_NLFV_SOLVER) {
    cout << endl 
         << "Test: bi-quartic solution with NLFV (nonlinear convergence)" << endl;
    Init();

    /* generate a dispersion matrix */
    const Epetra_Vector& phi = TS->ref_porosity();
    const Epetra_Vector& flux = TS->ref_darcy_flux();
    const Epetra_Vector& ws = TS->ref_water_saturation();

    Teuchos::RCP<Dispersion_TPFA> matrix = Teuchos::rcp(new Dispersion_TPFA(&(TPK->dispersion_models()), mesh, TS));
    matrix->Init();
    matrix->InitPreconditioner("Hypre AMG", plist.sublist("Preconditioners"));

    matrix->SymbolicAssembleGlobalMatrix();
    matrix->CalculateDispersionTensor(flux, phi, ws);
    matrix->InitNLFV();
    matrix->CreateFluxStencils();

    /* create matrix */
    Epetra_Vector u(phi);
    u.PutScalar(1.0);
    matrix->AssembleGlobalMatrixNLFV(u);
    matrix->AddTimeDerivative(1.0, phi, ws);

    /* populate right-hand side */
    Epetra_Vector r(phi), b(phi);
    InitRHS(b, 1.0);

    double bnorm; 
    b.Norm2(&bnorm);

    matrix->UpdatePreconditioner();

    AmanziSolvers::LinearOperatorFactory<Dispersion_TPFA, Epetra_Vector, Epetra_BlockMap> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Dispersion_TPFA, Epetra_Vector, Epetra_BlockMap> >
       solver = factory.Create("Dispersion Solver", plist.sublist("Solvers"), matrix);

    double snorm, residual; 
    for (int n = 0; n < 10; n++) {
      matrix->AssembleGlobalMatrixNLFV(u);
      matrix->AddTimeDerivative(1.0, phi, ws);

      u.Norm2(&snorm);
      residual = solver->TrueResidual(b, u);
      cout << "||r||=" << residual <<  "  ||u||=" << snorm << endl;
      CHECK(residual < 0.12 / pow(1.7, double(n)));

      solver->ApplyInverse(b, u);
    }

    // check spartial convergence
    Epetra_Vector err(u);
    for (int c = 0; c < u.MyLength(); c++) {
      const AmanziGeometry::Point& xm = mesh->cell_centroid(c);
      double x = xm[0];  
      double y = xm[1];  

      double tmp = x * y * (1.0 - x) * (1.0 - y);
      err[c] -= tmp * tmp;
    }

    double enorm;
    err.Norm2(&enorm);
    // CHECK(enorm / bnorm < 0.3);
  }

  /* **************************************************************** */
  TEST_FIXTURE(Problem, MATRIX_NLFV_PICARD) {
    cout << endl 
         << "Test: bi-quartic solution with NLFV (Picard)" << endl;
    Init();

    /* populate the solution guess and right-hand side */
    const Epetra_Vector& phi = TS->ref_porosity();
    const Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(phi));
    const Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(phi));

    u->PutScalar(1.0);
    InitRHS(*b, 1.0);

    /* create the function class */
    Teuchos::RCP<SolverFnNLFV<Epetra_Vector> > fn = 
        Teuchos::rcp(new SolverFnNLFV<Epetra_Vector>(mesh, TPK, b));

    // create the Solver
    const Epetra_BlockMap& map = b->Map();
    Teuchos::ParameterList& nlist = TPK->nonlin_solvers_list;

    Teuchos::RCP<AmanziSolvers::SolverNewton<Epetra_Vector, Epetra_BlockMap> > picard =
        Teuchos::rcp(new AmanziSolvers::SolverNewton<Epetra_Vector, Epetra_BlockMap>(nlist, fn, map));

    picard->Solve(u);
  }
}

