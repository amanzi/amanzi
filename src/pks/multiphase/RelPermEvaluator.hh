/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Relative permeability as a function of liquid saturation.
*/

#ifndef AMANZI_MULTIPHASE_REL_PERM_EVALUATOR_HH_
#define AMANZI_MULTIPHASE_REL_PERM_EVALUATOR_HH_

// Amanzi
#include "secondary_variable_field_evaluator.hh"

// Multiphase
#include "MultiphaseTypeDefs.hh"
#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

class RelPermEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  RelPermEvaluator(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<WRMmpPartition>& wrm);
  RelPermEvaluator(const RelPermEvaluator& other);

  // inteface functions to FieldEvaluator
  virtual Teuchos::RCP<FieldEvaluator> Clone() const override;

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) override;

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
      Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 private:
  Teuchos::RCP<WRMmpPartition> wrm_;
 
  std::string phase_name_;
  Key saturation_liquid_key_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
