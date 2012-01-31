#ifndef _EXPLICITTIFNBASE_HH_
#define _EXPLICITTIFNBASE_HH_


#include "TreeVector.hh"

namespace Amanzi {

// this is the interface definition for the explicit 
// Runge Kutta time integration class
class ExplicitTIfnBase {
 public:
  // computes the  functional f = f(t,u) 
  virtual void fun(const double t, const Amanzi::TreeVector& u, Amanzi::TreeVector& f) = 0;
};

}  // namespace Amanzi


#endif 
