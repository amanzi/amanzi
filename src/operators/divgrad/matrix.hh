/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

  Base Operators::Matrix class.

*/


#ifndef OPERATORS_MATRIX_HH_
#define OPERATORS_MATRIX_HH_


namespace Amanzi {

class TreeVector;

namespace Operators {


enum Matrix_bc {
  MATRIX_BC_NULL = 0,
  MATRIX_BC_DIRICHLET,
  MATRIX_BC_FLUX
};


class Matrix {

public:
  virtual void Apply(const TreeVector& X,
                     const Teuchos::Ptr<TreeVector>& Y) const = 0;
  virtual void ApplyInverse(const TreeVector& X,
                            const Teuchos::Ptr<TreeVector>& Y) const = 0;

};


} // namespace
} // namespace


#endif
