#include <ios>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "Mesh_STK.hh"
#include "MeshFactory.hh"
#include "CompositeVector.hh"
#include "TreeVector.hh"

#include "ImplicitTIBDF2fnBase.hh"
#include "ImplicitTIBDF2.hh"

using namespace Amanzi;

SUITE(ODEIntegrationTests) {
  // data structures for testing
  struct test_data {
    Epetra_MpiComm *comm;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;

    Teuchos::RCP<CompositeVector> init_vec;
    Teuchos::RCP<TreeVector> init;
    Teuchos::RCP<TreeVector> u;
    Teuchos::RCP<TreeVector> u_dot;
    Teuchos::RCP<TreeVector> u_ex;

    Teuchos::RCP<Teuchos::ParameterList> plist;

    test_data() {
      // comm and mesh for maps
      comm = new Epetra_MpiComm(MPI_COMM_SELF);
      AmanziMesh::MeshFactory mesh_fact(*comm);
      mesh = mesh_fact(0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 2, 1, 1);

      // plist for setting up integrator
      plist = Teuchos::rcp(new Teuchos::ParameterList());
      Teuchos::ParameterList& verblist = plist->sublist("VerboseObject");
      verblist.set("Verbosity Level","high");

      // non-ghosted u
      init_vec = Teuchos::rcp(new CompositeVector(mesh, "cell", AmanziMesh::CELL, 1, false));
      init = Teuchos::rcp(new TreeVector("init"));
      init->set_data(init_vec);

      // u, u_dot, and exact soln
      u = Teuchos::rcp(new TreeVector("u_dot", *init));
      u_dot = Teuchos::rcp(new TreeVector("u_dot", *init));
      u_ex = Teuchos::rcp(new TreeVector("u_ex", *init));
    }
    ~test_data() {
      delete comm;
    }
  };

  // ODE for testing
  class nonlinearODE : public Amanzi::ImplicitTIBDF2fnBase {
  public:
    nonlinearODE(double atol, double rtol) :
      rtol_(rtol), atol_(atol) {
    }

    void fun(double t, Teuchos::RCP<const Amanzi::TreeVector> u,
             Teuchos::RCP<const Amanzi::TreeVector> udot,
             Teuchos::RCP<Amanzi::TreeVector> f) {
      // f = udot - u^2
      // note that the exact solution is
      // uex = u0/(1-u0(t-t0))

      f->Multiply(1.0, *u, *u, 0.0);
      f->Update(1.0, *udot, -1.0);
    }

    void precon(Teuchos::RCP<const Amanzi::TreeVector> u, Teuchos::RCP<Amanzi::TreeVector> Pu) {
      *Pu = *u;
    }

    double enorm(Teuchos::RCP<const Amanzi::TreeVector> u, Teuchos::RCP<const Amanzi::TreeVector> du) {
      double norm_du, norm_u;
      du->NormInf(&norm_du);
      u->NormInf(&norm_u);
      return  norm_du/(atol_+rtol_*norm_u);
    }

    void update_precon(double t, Teuchos::RCP<const Amanzi::TreeVector> up, double h, int& errc) {
      // do nothing since the preconditioner is the identity
    }

    void compute_udot(double t,  Teuchos::RCP<const Amanzi::TreeVector> u,  Teuchos::RCP<const Amanzi::TreeVector> udot) {
    }

    double atol_, rtol_;
  };

  // Tests
  TEST_FIXTURE(test_data, NonlinearODE) {
    std::cout << "Test: NonlinearODE" << std::endl;

    // set the parameter list for BDF2
    plist->set("Nonlinear solver max iterations", 10);
    plist->set("Nonlinear solver tolerance", 0.01);
    plist->set("NKA max vectors",5);
    plist->set("NKA drop tolerance",0.01);

    // create the PDE problem
    nonlinearODE NF (1e-6, 1e-6);

    // create the time stepper
    Amanzi::ImplicitTIBDF2 TS(NF, init);
    TS.setParameterList(plist);

    // initial value
    u->PutScalar(-1.0);
    u_dot->PutScalar(1.0);

    // initial time
    double t=0.0;

    // final time
    double tout = 2.0;

    // initial time step
    double h = 1.0e-5;
    double hnext;

    // initialize the state of the time stepper
    TS.set_initial_state(t, u, u_dot);

    // iterate until the final time
    int i=0;
    double tlast = t;
    do {
      if (tlast + h > tout) {
        std::cout << "adjusting h, to hit the final time exactly...\n";
        h = tout - tlast;
      }
      TS.bdf2_step(h, 0.0, u, hnext);

      u->Print(std::cout);

      TS.commit_solution(h, u);

      TS.write_bdf2_stepping_statistics();

      h = hnext;
      i++;

      tlast=TS.most_recent_time();
    } while (tout > tlast);

    // compute the error with the exact solution
    u_ex->PutScalar(-1.0/3.0);
    u->Update(1.0, *u_ex, -1.0);

    double norm;
    u->NormInf(&norm);

    CHECK_CLOSE(0.0,norm,1e-4);
  }
}
