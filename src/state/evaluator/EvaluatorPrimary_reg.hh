/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "evaluator/EvaluatorPrimary.hh"

namespace Amanzi {

// registry of method
template <>
Utils::RegisteredFactory<Evaluator, EvaluatorPrimary<double>>
  EvaluatorPrimary<double>::fac_("primary variable double");

template <>
Utils::RegisteredFactory<
  Evaluator, EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>
  EvaluatorPrimary<CompositeVector, CompositeVectorSpace>::fac_(
    "primary variable");

} // namespace Amanzi
