/*
  Time Integration 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Markus Berndt (berndt@lanl.gov)
*/

#ifndef AMANZI_EXPLICIT_FNBASE_HH_
#define AMANZI_EXPLICIT_FNBASE_HH_

namespace Amanzi {
namespace Explicit_TI {

// this is the interface definition for the explicit
// Runge Kutta time integration class
template <class Vector>
class fnBase {
 public:
  // modifies solution before each call of functional f(t, u). Since no modifications
  // is often made for a low-order scheme, we provide the empty body. A DG scheme may
  // overload this function as the way to limit intermediate solutions.
  virtual void ModifySolution(const double t, Vector& u){};

  // computes functional f = f(t,u)
  virtual void FunctionalTimeDerivative(const double t, const Vector& u, Vector& f) = 0;
};

} // namespace Explicit_TI
} // namespace Amanzi


#endif
