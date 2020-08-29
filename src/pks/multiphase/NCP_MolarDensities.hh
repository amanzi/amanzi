/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Field evaluator for noninear complimentary problem, function G.
*/

#ifndef AMANZI_MULTIPHASE_NCP_MOLAR_DENSITIES_HH_
#define AMANZI_MULTIPHASE_NCP_MOLAR_DENSITIES_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiphaseBaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class NCP_MolarDensities : public MultiphaseBaseEvaluator {
 public:
  NCP_MolarDensities(Teuchos::ParameterList& plist);
  NCP_MolarDensities(const NCP_MolarDensities& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result);

  virtual void EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key x_vapor_key_, tcc_gas_key_, molar_density_gas_key_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
