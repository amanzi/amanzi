/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Multiphase PK

  Field evaluator for noninear complimentary problem, function F.
*/

#ifndef AMANZI_MULTIPHASE_NCP_F_HH_
#define AMANZI_MULTIPHASE_NCP_F_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiphaseBaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class NCP_F : public MultiphaseBaseEvaluator {
 public:
  NCP_F(Teuchos::ParameterList& plist);
  NCP_F(const NCP_F& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key saturation_liquid_key_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
