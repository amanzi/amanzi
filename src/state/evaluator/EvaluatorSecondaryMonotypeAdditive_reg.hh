/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/EvaluatorSecondaryMonotypeAdditive.hh"
namespace Amanzi {

template <>
Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeAdditive<
                                      CompositeVector, CompositeVectorSpace>>
  EvaluatorSecondaryMonotypeAdditive<CompositeVector,
                                     CompositeVectorSpace>::fac_("additive");

} // namespace Amanzi
