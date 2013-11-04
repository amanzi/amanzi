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

#include "State.hh"
#include "Transport_PK.hh"
#include "DispersionMatrixFactory.hh"
#include "Dispersion.hh"
#include "Dispersion_TPFA.hh"
#include "Dispersion_MFD.hh"
#include "Dispersion_NLFV.hh"

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

      /* create a simple state and populate it */
      Amanzi::VerboseObject::hide_line_prefix = true;

      std::vector<std::string> component_names;
      component_names.push_back("Component 0");

      S = rcp(new State());
      S->RegisterDomainMesh(mesh);

      TPK = rcp(new Transport_PK(plist, S, component_names));
      TPK->CreateDefaultState(mesh, 1);

      /* create RCP pointers */
      ws = S->GetFieldData("water_saturation")->ViewComponent("cell", false);
      phi = S->GetFieldData("porosity")->ViewComponent("cell", false);
      flux = S->GetFieldData("darcy_flux")->ViewComponent("face", true);

      /* modify the default state for the problem at hand */
      std::string passwd("state"); 
      Teuchos::RCP<Epetra_MultiVector> 
          darcy_flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

      AmanziGeometry::Point velocity(3.0, 4.0);
      int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
      for (int f = 0; f < nfaces_owned; f++) {
        const AmanziGeometry::Point& normal = mesh->face_normal(f);
        (*darcy_flux)[0][f] = velocity * normal;
      }

      S->GetFieldData("porosity", passwd)->PutScalar(0.5);

      /* initialize a transport process kernel */
      TPK->InitPK();
      TPK->TimeStep(1.0);
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
    RCP<const Mesh> mesh;
    ParameterList plist;

    RCP<State> S;
    RCP<Transport_PK> TPK;
    Teuchos::RCP<const Epetra_MultiVector> ws, phi, flux;
  };

  /* **************************************************************** */
  TEST_FIXTURE(Problem, MATRIX_TPFA) {
    cout << endl << "Test: TPFA matrix (no spatial convergence)" << endl;
    Init();

    /* generate a dispersion matrix */
    Teuchos::RCP<Dispersion_TPFA> 
        matrix = Teuchos::rcp(new Dispersion_TPFA(&(TPK->dispersion_models()), mesh, S));
    matrix->Init();
    matrix->InitPreconditioner("Hypre AMG", plist.sublist("Preconditioners"));
    matrix->SymbolicAssembleMatrix();
    matrix->CalculateDispersionTensor(*flux, *phi, *ws);
    matrix->AssembleMatrix(*phi);

    /* populate right-hand side and solution */
    const Epetra_BlockMap map = phi->Map();
    Epetra_Vector u(map), r(map), b(map);
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
  TEST_FIXTURE(Problem, MATRIX_MFD) {
    cout << endl << "Test: MFD matrix" << endl;
    Init();

    /* generate a dispersion matrix */
    Teuchos::RCP<Dispersion_MFD>
        matrix = Teuchos::rcp(new Dispersion_MFD(&(TPK->dispersion_models()), mesh, S));
    matrix->Init();
    matrix->InitPreconditioner("Hypre AMG", plist.sublist("Preconditioners"));
    matrix->SymbolicAssembleMatrix();
    matrix->CalculateDispersionTensor(*flux, *phi, *ws);
    matrix->AssembleMatrix(*phi);

    /* populate right-hand side and solution */
    const Epetra_Map& map = matrix->super_map();
    Epetra_Vector u(map), r(map), b(map);
    u.PutScalar(1.0);
    b.PutScalar(0.0);

    matrix->Apply(u, r);
    r.Update(-1.0, b, 1.0);

    double bnorm, residual; 
    b.Norm2(&bnorm);
    r.Norm2(&residual);
    cout << "Relative residual: " << residual << endl;
    CHECK(residual < 0.5);
  }

  /* **************************************************************** */
  TEST_FIXTURE(Problem, MATRIX_NLFV) {
    cout << endl << "Test: NLFV matrix (spartial convergence)" << endl;
    Init();

    /* generate a dispersion matrix */
    Teuchos::RCP<Dispersion_NLFV>
        matrix = Teuchos::rcp(new Dispersion_NLFV(&(TPK->dispersion_models()), mesh, S));
    matrix->Init();
    matrix->InitNLFV();
    matrix->CalculateDispersionTensor(*flux, *phi, *ws);

    matrix->SymbolicAssembleMatrix();
    matrix->ModifySymbolicAssemble();

    /* populate solution and right-hand side */
    const Epetra_BlockMap map = phi->Map();
    Epetra_Vector u(map), b(map);
    InitSOL(u);
    InitRHS(b, 1.0);

    /* update matrix */
    matrix->AssembleMatrix(u);
    matrix->AddTimeDerivative(1.0, *phi, *ws);

    /* Compute residual */
    Epetra_Vector r(map);
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
    cout << endl << "Test: Factory of matrices" << endl;
    Init();

    /* generate a dispersion matrix */
    DispersionMatrixFactory factory;
    Teuchos::RCP<Dispersion> matrix = factory.Create("tpfa", &(TPK->dispersion_models()), mesh, S);
    matrix->InitPreconditioner("Hypre AMG", plist.sublist("Preconditioners"));

    matrix->CalculateDispersionTensor(*flux, *phi, *ws);
    matrix->SymbolicAssembleMatrix();
    matrix->ModifySymbolicAssemble();

    /* populate right-hand side and solution */
    const Epetra_BlockMap& map = phi->Map();
    Epetra_Vector u(map), r(map), b(map);
    InitSOL(u);
    InitRHS(b, 0.0);

    matrix->AssembleMatrix(u);

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
  TEST_FIXTURE(Problem, MATRIX_PICARD) {
    cout << endl << "Test: Nonlinear convergence" << endl;
    Init();

    /* generate a dispersion matrix */
    Teuchos::RCP<Dispersion_NLFV>
        matrix = Teuchos::rcp(new Dispersion_NLFV(&(TPK->dispersion_models()), mesh, S));
    matrix->Init();
    matrix->InitNLFV();

    matrix->CalculateDispersionTensor(*flux, *phi, *ws);
    matrix->SymbolicAssembleMatrix();
    matrix->ModifySymbolicAssemble();

    /* create matrix */
    const Epetra_BlockMap& map = phi->Map();
    Epetra_Vector u(map);
    u.PutScalar(1.0);

    matrix->AssembleMatrix(u);
    matrix->AddTimeDerivative(1.0, *phi, *ws);

    /* populate right-hand side */
    Epetra_Vector r(map), b(map);
    InitRHS(b, 1.0);

    double bnorm; 
    b.Norm2(&bnorm);

    matrix->InitPreconditioner("Hypre AMG", plist.sublist("Preconditioners"));
    matrix->UpdatePreconditioner();

    AmanziSolvers::LinearOperatorFactory<Dispersion_NLFV, Epetra_Vector, Epetra_BlockMap> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Dispersion_NLFV, Epetra_Vector, Epetra_BlockMap> >
        solver = factory.Create("Dispersion Solver", plist.sublist("Solvers"), matrix);

    double snorm, residual; 
    for (int n = 0; n < 10; n++) {
      matrix->AssembleMatrix(u);
      matrix->AddTimeDerivative(1.0, *phi, *ws);

      u.Norm2(&snorm);
      residual = solver->TrueResidual(b, u);
      cout << "||r||=" << residual << "  ||u||=" << snorm << endl;
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
  TEST_FIXTURE(Problem, MATRIX_SOLVER) {
    cout << endl << "Test: NLFV coupled with Newton" << endl;
    Init();

    /* populate the solution guess and right-hand side */
    const Epetra_BlockMap& map = phi->Map();
    const Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(map));
    const Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(map));

    u->PutScalar(1.0);
    InitRHS(*b, 1.0);

    /* create matrix sparcity structure */
    DispersionMatrixFactory factory;
    Teuchos::RCP<Dispersion> matrix = factory.Create("nlfv", &(TPK->dispersion_models()), mesh, S);
    TPK->init_dispersion_matrix(matrix);
    matrix->InitPreconditioner("Hypre AMG", plist.sublist("Preconditioners"));

    matrix->CalculateDispersionTensor(*flux, *phi, *ws);
    matrix->SymbolicAssembleMatrix();
    matrix->ModifySymbolicAssemble();

    /* create the function class */
    Teuchos::RCP<SolverFnNLFV<Epetra_Vector> > fn = 
        Teuchos::rcp(new SolverFnNLFV<Epetra_Vector>(mesh, TPK, b));

    /* create nonlinear solver */
    Teuchos::ParameterList& nlist = TPK->nonlin_solvers_list;

    Teuchos::RCP<AmanziSolvers::SolverNewton<Epetra_Vector, Epetra_BlockMap> > picard =
        Teuchos::rcp(new AmanziSolvers::SolverNewton<Epetra_Vector, Epetra_BlockMap>(nlist, fn, map));

    picard->Solve(u);
  }
}

