/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//! EvaluatorAlgebraic is algebraic in data from the same time tag.

/*!

Algebraic evaluators are secondary evaluators that read only dependencies of
the same type as they calculate.  This allows requirements placed on the
calculated variable to be pushed down to the dependencies, checking
consistency, and also allows derivatives to be calculated automatically.

Algebraic variable evaluators, such as equations of state, water retention
evaluators, internal energy evaluators, etc should inherit this class,
implementing the missing Update_() and UpdateFieldDerivative_() methods.

*/

#ifndef STATE_EVALUATOR_ALGEBRAIC_HH_
#define STATE_EVALUATOR_ALGEBRAIC_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "exceptions.hh"

#include "EvaluatorSecondary.hh"

namespace Amanzi {

// By default, this class adds nothing on top of EvaluatorSecondary.
// Specializations can do useful things though.
template <typename Data_t, typename DataFactory_t = NullFactory>
class EvaluatorAlgebraic : public EvaluatorSecondary<Data_t, DataFactory_t> {
public:
  using EvaluatorSecondary<Data_t, DataFactory_t>::EvaluatorSecondary;

  virtual void EnsureCompatibility(State &S) override {
    EvaluatorSecondary<Data_t, DataFactory_t>::EnsureCompatibility(S);
  }

protected:
  virtual void UpdateDerivative_(State &S, const Key &wrt_key,
                                 const Key &wrt_tag) override {
    EvaluatorSecondary<Data_t, DataFactory_t>::UpdateDerivative_(S, wrt_key,
                                                                 wrt_tag);
  }
};

template <> void EvaluatorAlgebraic<double>::EnsureCompatibility(State &S);

template <>
void EvaluatorAlgebraic<double>::UpdateDerivative_(State &S, const Key &wrt_key,
                                                   const Key &wrt_tag);

template <>
void EvaluatorAlgebraic<CompositeVector,
                        CompositeVectorSpace>::EnsureCompatibility(State &S);

template <>
void EvaluatorAlgebraic<CompositeVector, CompositeVectorSpace>::
    UpdateDerivative_(State &S, const Key &wrt_key, const Key &wrt_tag);

} // namespace Amanzi

#endif
