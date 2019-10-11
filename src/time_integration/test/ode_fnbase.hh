#ifndef ODE_TEST_FNBASE_HH_
#define ODE_TEST_FNBASE_HH_

#include "BDFFnBase.hh"
#include "FnBaseDefs.hh"

using namespace Amanzi;


// ODE for testing
class nonlinearODE : public Amanzi::BDFFnBase<Vector_type> {
 public:
  nonlinearODE(double atol, double rtol, bool exact_jacobian)
    : rtol_(rtol), atol_(atol), exact_jacobian_(exact_jacobian)
  {}

  void
  FunctionalResidual(double told, double tnew, Teuchos::RCP<Vector_type> uold,
                     Teuchos::RCP<Vector_type> unew,
                     Teuchos::RCP<Vector_type> f)
  {
    // f = udot - u^2
    // note that the exact solution is
    // uex = u0/(1-u0(t-t0))

    // compute udot... f <-- (unew-uold)/(tnew-told)
    f->assign(*unew);
    double hinv(1.0 / (tnew - told));
    f->update(-hinv, *uold, hinv);

    // f <-- f - unew^2
    f->elementWiseMultiply(-1.0, *unew, *unew, 1.0);

    {
      auto u0v = uold->getLocalViewHost();
      auto u1v = unew->getLocalViewHost();
      auto fv = f->getLocalViewHost();
      std::cout.precision(10);
      for (int c = 0; c != unew->getLocalLength(); ++c) {
        std::cout << "Res: u_old = " << u0v(c, 0) << ", u_new = " << u1v(c, 0)
                  << ", f = " << fv(c, 0) << std::endl;
      }
    }
  }

  int ApplyPreconditioner(Teuchos::RCP<const Vector_type> u,
                          Teuchos::RCP<Vector_type> Pu)
  {
    Pu->elementWiseMultiply(1., *Pu_, *u, 0.);
    return 0;
  }

  double ErrorNorm(Teuchos::RCP<const Vector_type> u,
                   Teuchos::RCP<const Vector_type> du)
  {
    double norm_du, norm_u;
    norm_du = du->normInf();
    norm_u = u->normInf();
    return norm_du / (atol_ + rtol_ * norm_u);
  }

  void
  UpdatePreconditioner(double t, Teuchos::RCP<const Vector_type> up, double h)
  {
    // do nothing since the preconditioner is the identity
    if (Pu_ == Teuchos::null) {
      Pu_ = Teuchos::rcp(new Vector_type(up->getMap()));
    }

    if (exact_jacobian_) {
      Pu_->assign(*up);
      Vector_type ones(Pu_->getMap());
      ones.putScalar(1.0);
      Pu_->update(1.0 / h, ones, -2.0);
    } else {
      Pu_->putScalar(1.0 / h);
    }
    Pu_->reciprocal(*Pu_);
  }

  void compute_udot(double t, Teuchos::RCP<const Vector_type> u,
                    Teuchos::RCP<const Vector_type> udot){};

  bool IsAdmissible(Teuchos::RCP<const Vector_type> up) { return true; }
  bool ModifyPredictor(double h, Teuchos::RCP<const Vector_type> u0,
                       Teuchos::RCP<Vector_type> u)
  {
    return false;
  }
  Amanzi::AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const Vector_type> res,
                   Teuchos::RCP<const Vector_type> u,
                   Teuchos::RCP<Vector_type> du)
  {
    return Amanzi::AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }
  void ChangedSolution() {}


  bool exact_jacobian_;
  double atol_, rtol_;
  Teuchos::RCP<Vector_type> Pu_;
};


#endif
