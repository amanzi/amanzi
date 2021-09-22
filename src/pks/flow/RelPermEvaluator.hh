/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Relative permeability as a function of capillary pressure, k=k(pc).
*/

#ifndef AMANZI_FLOW_REL_PERM_EVALUATOR_HH_
#define AMANZI_FLOW_REL_PERM_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"
#include "RelPerm.hh"
#include "WRM.hh"
#include "WRMPartition.hh"

namespace Amanzi {
namespace Flow {

class RelPermEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  RelPermEvaluator(Teuchos::ParameterList& plist,
                   const Teuchos::Ptr<State>& S,
                   const Teuchos::RCP<WRMPartition>& wrm);
  RelPermEvaluator(const RelPermEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
      Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  void InitializeFromPlist_(const Teuchos::Ptr<State>& S);

 protected:
  Teuchos::RCP<RelPerm> relperm_;

  Teuchos::RCP<WRMPartition> wrm_;
  Key pressure_key_;

  double patm_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
