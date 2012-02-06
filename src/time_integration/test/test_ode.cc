#include <ios>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"
#include "Epetra_BlockMap.h"

#include "TreeVector.hh"

#include "ImplicitTIBDF2fnBase.hh"
#include "ImplicitTIBDF2.hh"


TEST(NonlinearODE) {

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


  std::cout << "Test: NonlinearODE" << std::endl;

  // create the parameter list for BDF2
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::rcp(new Teuchos::ParameterList());
  plist->set("Nonlinear solver max iterations", 10);
  plist->set("Nonlinear solver tolerance", 0.01);
  plist->set("NKA max vectors",5);
  plist->set("NKA drop tolerance",0.01);
  Teuchos::ParameterList& verblist = plist->sublist("VerboseObject");

  verblist.set("Verbosity Level","high");

  Epetra_Comm* comm = new Epetra_SerialComm();


  Epetra_BlockMap map(1,1,0, *comm);

  Teuchos::RCP<Epetra_MultiVector> x = Teuchos::rcp(new Epetra_MultiVector(map,1,false));
  Teuchos::RCP<Amanzi::TreeVector> initvec = Teuchos::rcp( new Amanzi::TreeVector("initvector"));
  initvec->PushBack(x);

  // create the PDE problem
  nonlinearODE NF (1e-6, 1e-6);

  // create the time stepper
  Amanzi::ImplicitTIBDF2 TS(NF, initvec);
  TS.setParameterList(plist);

  // create the initial condition
  Teuchos::RCP<Amanzi::TreeVector> u_tree = Teuchos::rcp( new Amanzi::TreeVector("u", initvec));

  // initial value
  u_tree->PutScalar(-1.0);

  // initial time
  double t=0.0;

  // final time
  double tout = 2.0;

  // create udot and compute its initial value
  Teuchos::RCP<Amanzi::TreeVector> udot_tree = Teuchos::rcp( new Amanzi::TreeVector("udot", initvec));
  // initial time derivative
  udot_tree->PutScalar(1.0);

  // initial time step
  double h = 1.0e-5;
  double hnext;

  // initialize the state of the time stepper
  TS.set_initial_state(t, u_tree, udot_tree);

  // iterate until the final time
  int i=0;
  double tlast = t;
  do {
    if (tlast + h > tout) {
      std::cout << "adjusting h, to hit the final time exactly...\n";
      h = tout - tlast;
    }
    TS.bdf2_step(h,0.0,u_tree,hnext);

    u_tree->Print(std::cout);

    TS.commit_solution(h,u_tree);

    TS.write_bdf2_stepping_statistics();

    h = hnext;
    i++;

    tlast=TS.most_recent_time();

  } while (tout > tlast);

  delete comm;

  // compute the error with the exact solution
  Teuchos::RCP<Amanzi::TreeVector> uex_tree = Teuchos::rcp( new Amanzi::TreeVector("u", initvec));
  uex_tree->PutScalar(-1.0/3.0);

  u_tree->Update(1.0, *uex_tree, -1.0);

  double norm;
  u_tree->NormInf(&norm);

  CHECK_CLOSE(0.0,norm,1e-4);

}
