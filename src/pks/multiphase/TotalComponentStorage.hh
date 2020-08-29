/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for a total component stirage (water, hydrogen,
  etc) storage, the conserved quantity:

    TCS = phi * (rho_l * s_l * X_l + rho_g * s_g * X_g)

  where x_p is the mole fraction of a component in phase p.
*/


#ifndef AMANZI_MULTIPHASE_TOTAL_COMPONENT_STORAGE_HH_
#define AMANZI_MULTIPHASE_TOTAL_COMPONENT_STORAGE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiphaseBaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class TotalComponentStorage : public MultiphaseBaseEvaluator {
 public:
  TotalComponentStorage(Teuchos::ParameterList& plist);
  TotalComponentStorage(const TotalComponentStorage& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void Init_();

  // interface to FieldEvaluator
  virtual void EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result);

  virtual void EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key saturation_liquid_key_, porosity_key_;
  Key molar_density_liquid_key_, molar_density_gas_key_;
  Key mole_fraction_liquid_key_, mole_fraction_gas_key_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator, TotalComponentStorage> reg_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
