/*
  This is the Nonlinear Solver component of the Amanzi code.

  Example of a Vector base class used in our templates. The 
  routines below are mandatory for any implementation of a 
  vector class to be compatible with Amanzi.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#ifndef AMANZI_VECTOR_BASE_
#define AMANZI_VECTOR_BASE_

#include "Teuchos_RCP.hpp"
#include "VectorSpaceBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

class VectorBase {
 public:
  virtial VectorBase(const VectorBase& u) = 0;

  // Underlying vector map.
  virtual VectorSpace& Map() = 0;

  // (*this) = a * (*this) + b * u
  virtual Update(double b, const VectorBase& u, double a);

  // Eucleadian norm of vector
  virtual Norm2(double* norm);

  // dot product: a = (*this) * u
  Dot(const VectorBase& u, double* a);
};

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
