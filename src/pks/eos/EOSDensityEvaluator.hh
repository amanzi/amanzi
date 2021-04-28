/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOSFieldEvaluator is the interface between state/data and the model, an EOS.
*/

#ifndef AMANZI_EOS_EVALUATOR_HH_
#define AMANZI_EOS_EVALUATOR_HH_

#include "EOS_Density.hh"
#include "Factory.hh"
#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace AmanziEOS {

class EOSDensityEvaluator : public SecondaryVariablesFieldEvaluator {
 public:
  enum EOSMode { EOS_MODE_MASS, EOS_MODE_MOLAR, EOS_MODE_BOTH };

  explicit EOSDensityEvaluator(Teuchos::ParameterList& plist);

  EOSDensityEvaluator(const EOSDensityEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results);

  Teuchos::RCP<EOS_Density> get_EOS() { return eos_; }

 protected:
  // the actual model
  Teuchos::RCP<EOS_Density> eos_;
  EOSMode mode_;

  // Keys for field dependencies
  Key temp_key_;
  Key pres_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator, EOSDensityEvaluator> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
