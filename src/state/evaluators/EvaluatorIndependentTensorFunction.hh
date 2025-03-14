/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

/*!

This evaluator is typically used for providing tensors that are functions of space
and time.  The evaluator consists of a list of (region, function) pairs, and the
functions are evaluated across that region at each timestep.  If the problem is
time-independent, the `"constant in time`" option results in a performance
improvement (as the functions need only be evaluated once).  This leverages the
exaustive functional format capability provided in Amanzi's Functions_ library.

This evaluator is used by providing the option:

`"evaluator type`" == `"independent variable tensor`"

.. _independent-variable-tensor-function-evaluator-spec:
.. admonition:: independent-variable-tensor-function-evaluator-spec

   * `"tensor type`" ``[string]`` One of:
     - `"scalar`" A single scalar, isotropic (1 DoF)
     - `"horizontal and vertical`" Diagonal with identical x-y entries and
       different z entries. (2 DoF)
     - `"diagonal`" Diagonal (3 DoF for 3D, 2 for 2D)
     - `"full symmetric`" Full (symmetric) tensor, (6 DoF for 3D, 3 for 2D)
     - `"full`" Full (nonsymmetric) tensor, (9 DoF for 3D, 6 for 2D)
   * `"constant in time`" ``[bool]`` **false** If true, only evaluate the
     functions once as they are time-independent.
   * `"rescaling`" ``[double]`` **1** Scales the tensors by a uniform, constant
     scalar.
   * `"function`" ``[composite-vector-function-spec-list]`` describes the
     function.  The number of degrees of freedom required in this function are
     given by the tensor type.

*/


#pragma once

#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"
#include "TensorVector.hh"

namespace Amanzi {

namespace Functions {
class CompositeVectorFunction;
}

class EvaluatorIndependentTensorFunction
  : public EvaluatorIndependent<TensorVector, TensorVector_Factory> {
 public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit EvaluatorIndependentTensorFunction(Teuchos::ParameterList& plist);
  EvaluatorIndependentTensorFunction(const EvaluatorIndependentTensorFunction& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void EnsureCompatibility(State& S) override;

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override;

 protected:
  Teuchos::RCP<Functions::CompositeVectorFunction> func_;
  int num_funcs_;
  int dimension_;
  std::string tensor_type_;
  int rank_;
  double rescaling_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentTensorFunction> fac_;
};


namespace Impl {
void CopyVectorToTensorVector(const Epetra_MultiVector& v, int j, TensorVector& tv);
}


} // namespace Amanzi
