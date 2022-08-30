/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Field evaluator for noninear complimentary problem, function G.
*/

#ifndef AMANZI_MULTIPHASE_NCP_MOLE_FRACTIONS_HH_
#define AMANZI_MULTIPHASE_NCP_MOLE_FRACTIONS_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiphaseBaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class NCP_MoleFractions : public MultiphaseBaseEvaluator {
 public:
  NCP_MoleFractions(Teuchos::ParameterList& plist);
  NCP_MoleFractions(const NCP_MoleFractions& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key x_vapor_key_, x_gas_key_;

  static Utils::RegisteredFactory<Evaluator, NCP_MoleFractions> fac_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
