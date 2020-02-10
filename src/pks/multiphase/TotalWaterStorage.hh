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

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void Init_();

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result);

  virtual void EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key saturation_liquid_key_, porosity_key_, x_vapor_key_;
  Key molar_density_liquid_key_, molar_density_gas_key_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator, TotalWaterStorage> reg_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
