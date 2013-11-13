#ifndef _BDFFNBASE_HH_
#define _BDFFNBASE_HH_


namespace Amanzi {

// this is the interface definition for the BDF class
// the nonlinear functional, preconditioner, and error
// functions must be derived from this class to be
// usable with BDF2::Dae or BDF1::Dae

template<class Vec>
class BDFFnBase {

public:

  // computes the non-linear functional f = f(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<Vec> u_old,
                   Teuchos::RCP<Vec> u_new, Teuchos::RCP<Vec> f) = 0;

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const Vec> u, Teuchos::RCP<Vec> Pu) = 0;

  // computes a norm on u-du and returns the result
  virtual double enorm(Teuchos::RCP<const Vec> u, Teuchos::RCP<const Vec> du) = 0;

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const Vec> up, double h) = 0;

  // check the admissibility of a solution
  // override with the actual admissibility check
  virtual bool is_admissible(Teuchos::RCP<const Vec> up) = 0;

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool modify_predictor(double h, Teuchos::RCP<Vec> up) = 0;

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual bool modify_correction(double h, Teuchos::RCP<const Vec> res,
          Teuchos::RCP<const Vec> u, Teuchos::RCP<Vec> du) = 0;

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void changed_solution() = 0;

};

}

#endif  // _BDFFNBASE_HH_
