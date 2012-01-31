#ifndef _IMPLICITTIBDF2FNBASE_HH_
#define _IMPLICITTIBDF2FNBASE_HH_

#include "TreeVector.hh"

namespace Amanzi {

  // this is the interface definition for the BDF2 class
  // the nonlinear functional, preconditioner, and error 
  // functions must be derived from this class to be 
  // usable with BDF2::Dae

  class ImplicitTIBDF2fnBase {
    
  public:
        
    // computes the non-linear functional f = f(t,u,udot) 
    virtual void fun(const double t, const Teuchos::RCP<Amanzi::TreeVector> u, 
		     const Teuchos::RCP<Amanzi::TreeVector> udot, Teuchos::RCP<Amanzi::TreeVector> f) = 0;

    // applies preconditioner to u and returns the result in Pu
    virtual void precon(const Teuchos::RCP<Amanzi::TreeVector> u, Teuchos::RCP<Amanzi::TreeVector> Pu) = 0;

    // computes a norm on u-du and returns the result 
    virtual double enorm(const Teuchos::RCP<Amanzi::TreeVector> u, Teuchos::RCP<Amanzi::TreeVector> du) = 0;

    // updates the preconditioner
    virtual void update_precon(const double t, const Teuchos::RCP<Amanzi::TreeVector> up, const double h, int& errc) = 0;
   
    // check the admissibility of a solution
    // override with the actual admissibility check
    virtual bool is_admissible(Teuchos::RCP<Amanzi::TreeVector> up) { return true; }

  };

}

#endif  // _IMPLICITTIBDF2FNBASE_HH_
