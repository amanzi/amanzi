#ifndef _EXPLICIT_FNBASE_HPP_
#define _EXPLICIT_FNBASE_HPP_

#include "Epetra_Vector.h"

namespace Amanzi {
namespace Explicit_TI {

// this is the interface definition for the explicit 
// Runge Kutta time integration class
class fnBase {
 public:
  // computes the  functional f = f(t,u) 
  virtual void fun(const double t, const Epetra_Vector& u, Epetra_Vector& f) = 0;
};

}  // namespace Explicit_TI
}  // namespace Amanzi


#endif 
