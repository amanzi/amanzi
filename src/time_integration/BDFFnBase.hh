#ifndef _BDFFNBASE_HH_
#define _BDFFNBASE_HH_

#include "TreeVector.hh"

namespace Amanzi {

  // this is the interface definition for the BDF2 class
  // the nonlinear functional, preconditioner, and error
  // functions must be derived from this class to be
  // usable with BDF2::Dae

  class BDFFnBase {

  public:

    // computes the non-linear functional f = f(t,u,udot)
    virtual void fun(double told, double tnew, Teuchos::RCP<Amanzi::TreeVector> uold,
                     Teuchos::RCP<Amanzi::TreeVector> unew, Teuchos::RCP<Amanzi::TreeVector> f) = 0;

    // applies preconditioner to u and returns the result in Pu
    virtual void precon(Teuchos::RCP<const Amanzi::TreeVector> u, Teuchos::RCP<Amanzi::TreeVector> Pu) = 0;

    // computes a norm on u-du and returns the result
    virtual double enorm(Teuchos::RCP<const Amanzi::TreeVector> u, Teuchos::RCP<const Amanzi::TreeVector> du) = 0;

    // updates the preconditioner
    virtual void update_precon(double t, Teuchos::RCP<const Amanzi::TreeVector> up, double h) = 0;

    // check the admissibility of a solution
    // override with the actual admissibility check
    virtual bool is_admissible(Teuchos::RCP<const Amanzi::TreeVector> up) {
      return true;
    }

  };

}

#endif  // _BDFFNBASE_HH_
