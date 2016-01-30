/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for total volumetric water content which is the 
  conserved quantity density in Richards's equation.

  Density is constant.
*/


#ifndef AMANZI_FLOW_VWCONTENT_EVALUATOR_CONST_DENSITY_HH_
#define AMANZI_FLOW_VWCONTENT_EVALUATOR_CONST_DENSITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "VWContentEvaluator.hh"

namespace Amanzi {
namespace Flow {

class VWContentEvaluator_ConstDensity : public VWContentEvaluator {
 public:
  explicit VWContentEvaluator_ConstDensity(Teuchos::ParameterList& plist);

  virtual void Init_();

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);
};

}  // namespace Flow
}  // namespace Amanzi

#endif
