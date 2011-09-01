#ifndef _EXPLICIT_FNBASE_HPP_
#define _EXPLICIT_FNBASE_HPP_

#include "Epetra_MultiVector.h"

namespace Explicit_TI {

  // this is the interface definition for the explicit 
  // Runge Kutta time integration class

  class fnBase {
    
  public:
        
    // computes the  functional f = f(t,u) 
    virtual void fun(const double t, const Epetra_MultiVector& u, Epetra_MultiVector& f) = 0;
  };

}

#endif  // _EXPLICIT_FNBASE_HPP_
