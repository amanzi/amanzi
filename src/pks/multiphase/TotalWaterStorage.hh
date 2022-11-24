/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Field evaluator for total water storage, the conserved quantity:

    TCS = phi * rho_l * s_l.
*/


#ifndef AMANZI_MULTIPHASE_TOTAL_WATER_STORAGE_HH_
#define AMANZI_MULTIPHASE_TOTAL_WATER_STORAGE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiphaseBaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class TotalWaterStorage : public MultiphaseBaseEvaluator {
 public:
  TotalWaterStorage(Teuchos::ParameterList& plist);
  TotalWaterStorage(const TotalWaterStorage& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void Init_();

 protected:
  Key saturation_liquid_key_, porosity_key_, x_vapor_key_;
  Key mol_density_liquid_key_, mol_density_gas_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, TotalWaterStorage> fac_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
