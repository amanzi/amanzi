/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Multiphase PK

  Field evaluator for a total component stirage (water, hydrogen,
  etc) storage, the conserved quantity:

    TCS = phi * (rho_l * s_l * x_l + rho_g * s_g * x_g)

  where x_p is the mole fraction of a component in phase p.
*/


#ifndef AMANZI_MULTIPHASE_TOTAL_COMPONENT_STORAGE_HH_
#define AMANZI_MULTIPHASE_TOTAL_COMPONENT_STORAGE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiphaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class TotalComponentStorage : public MultiphaseEvaluator {
 public:
  TotalComponentStorage(Teuchos::ParameterList& plist);
  TotalComponentStorage(const TotalComponentStorage& other);

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
  Key x_liquid_key_, x_gas_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, TotalComponentStorage> fac_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
