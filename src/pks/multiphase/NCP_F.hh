/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

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

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result);

  virtual void EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key saturation_liquid_key_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
