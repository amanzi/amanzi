/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Field evaluator for conserved quantity: the total component storage

    TCS = phi * (s_l * C_l + s_g * C_g)

  where C_p is the concentration a component in phase p.
*/


#ifndef AMANZI_MULTIPHASE_TOTAL_COMPONENT_STORAGE_TCC_HH_
#define AMANZI_MULTIPHASE_TOTAL_COMPONENT_STORAGE_TCC_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiphaseBaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class TotalComponentStorage_Tcc : public MultiphaseBaseEvaluator {
 public:
  TotalComponentStorage_Tcc(Teuchos::ParameterList& plist);
  TotalComponentStorage_Tcc(const TotalComponentStorage_Tcc& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key saturation_liquid_key_, porosity_key_;
  Key tcc_liquid_key_, tcc_gas_key_;
  
 private:
  static Utils::RegisteredFactory<Evaluator, TotalComponentStorage_Tcc> fac_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
