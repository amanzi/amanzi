/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Multiphase PK

  Field evaluator for a total component storage (water, hydrogen,
  etc) storage, the conserved quantity:

    TCS = phi * (rho_l * s_l + rho_g * s_g)
*/


#ifndef AMANZI_MULTIPHASE_TOTAL_COMPONENT_STORAGE_JAFFRE_HH_
#define AMANZI_MULTIPHASE_TOTAL_COMPONENT_STORAGE_JAFFRE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiphaseBaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class TotalComponentStorage_Jaffre : public MultiphaseBaseEvaluator {
 public:
  TotalComponentStorage_Jaffre(Teuchos::ParameterList& plist);
  TotalComponentStorage_Jaffre(const TotalComponentStorage_Jaffre& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key saturation_liquid_key_, porosity_key_;
  Key mol_density_liquid_key_, mol_density_gas_key_;

  static Utils::RegisteredFactory<Evaluator, TotalComponentStorage_Jaffre> fac_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
