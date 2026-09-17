/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/
/*!

Algebraic evaluator with detached dependencies.

*/

#ifndef AMANZI_STATE_EVALUATOR_SECONDARY_MONOTYPE_DETACHED_HH_
#define AMANZI_STATE_EVALUATOR_SECONDARY_MONOTYPE_DETACHED_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "exceptions.hh"
#include "Debugger.hh"

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class EvaluatorSecondaryMonotypeDetached
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit EvaluatorSecondaryMonotypeDetached(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist) {};

 protected:
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override;

 protected:
  KeyTagSet detached_dependencies_;
};

} // namespace Amanzi

#endif
