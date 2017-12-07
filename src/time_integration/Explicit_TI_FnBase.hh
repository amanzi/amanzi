#ifndef AMANZI_EXPLICIT_FNBASE_HH_
#define AMANZI_EXPLICIT_FNBASE_HH_

namespace Amanzi {
namespace Explicit_TI {

// this is the interface definition for the explicit 
// Runge Kutta time integration class
template<class Vector>
class fnBase {
 public:
  // computes the  functional f = f(t,u) 
  virtual void Dudt(double t, const Vector& u, Vector& f) = 0;
};

}  // namespace Explicit_TI
}  // namespace Amanzi


#endif 
