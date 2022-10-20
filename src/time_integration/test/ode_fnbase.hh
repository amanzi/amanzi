#ifndef ODE_TEST_FNBASE_HH_
#define ODE_TEST_FNBASE_HH_

#include "BDFFnBase.hh"
#include "FnBaseDefs.hh"

// ODE for testing
class nonlinearODE : public Amanzi::BDFFnBase<Epetra_Vector> {
public:
  nonlinearODE(double atol, double rtol, bool exact_jacobian) :
      exact_jacobian_(exact_jacobian), atol_(atol), rtol_(rtol) {
  }

  void FunctionalResidual(double told, double tnew, Teuchos::RCP<Epetra_Vector> uold,
                  Teuchos::RCP<Epetra_Vector> unew,
                  Teuchos::RCP<Epetra_Vector> f) {
    // f = udot - u^2
    // note that the exact solution is
    // uex = u0/(1-u0(t-t0))

    // compute udot... f <-- (unew-uold)/(tnew-told)
    *f = *unew;
    double hinv(1.0 / (tnew - told));
    f->Update(-hinv, *uold, hinv);

    // f <-- f - unew^2
    f->Multiply(-1.0, *unew, *unew, 1.0);

    std::cout.precision(10);
    for (int c = 0; c != unew->MyLength(); ++c) {
      std::cout << "Res: u_old = " << (*uold)[c] << ", u_new = " << (*unew)[c] << ", f = " << (*f)[c] << std::endl;
    }
  }

  int ApplyPreconditioner(Teuchos::RCP<const Epetra_Vector> u, Teuchos::RCP<Epetra_Vector> Pu) {
    return Pu->ReciprocalMultiply(1., *Pu_, *u, 0.);
  }

  double ErrorNorm(const Teuchos::RCP<const Epetra_Vector>& u,
                   const Teuchos::RCP<const Epetra_Vector>& du,
                   const Teuchos::RCP<const Epetra_Vector>& res,
                   const AmanziSolvers::ConvergenceMonitor& monitor) {
    double norm_du, norm_u;
    du->NormInf(&norm_du);
    u->NormInf(&norm_u);
    return  norm_du/(atol_+rtol_*norm_u);
  }

  void UpdatePreconditioner(double t, Teuchos::RCP<const Epetra_Vector> up, double h) {
    // do nothing since the preconditioner is the identity
    if (Pu_ == Teuchos::null) {
      Pu_ = Teuchos::rcp(new Epetra_Vector(*up));
    }

    if (exact_jacobian_) {
      *Pu_ = *up;
      Teuchos::RCP<Epetra_MultiVector> ones = Teuchos::rcp(new Epetra_Vector(*Pu_));
      ones->PutScalar(1.0);
      Pu_->Update(1.0/h, *ones, -2.0);
    } else {
      Pu_->PutScalar(1.0/h);
    }
  }


  void compute_udot(double t, Teuchos::RCP<const Epetra_Vector> u, Teuchos::RCP<const Epetra_Vector> udot) {};


  bool IsAdmissible(Teuchos::RCP<const Epetra_Vector> up) { return true; }
  bool ModifyPredictor(double h, Teuchos::RCP<const Epetra_Vector> u0, 
                       Teuchos::RCP<Epetra_Vector> u) { return false; }
  Amanzi::AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const Epetra_Vector> res,
                       Teuchos::RCP<const Epetra_Vector> u,
                       Teuchos::RCP<Epetra_Vector> du) { 
    return Amanzi::AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }
  void ChangedSolution() {}


  bool exact_jacobian_;
  double atol_, rtol_;
  Teuchos::RCP<Epetra_Vector> Pu_;
};



#endif
