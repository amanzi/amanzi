/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  The IEM Evaluator simply calls the IEM with the correct arguments.
*/

#ifndef AMANZI_ENERGY_RELATIONS_IEM_EVALUATOR_
#define AMANZI_ENERGY_RELATIONS_IEM_EVALUATOR_

#include "factory.hh"
#include "iem.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {

class IEMEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  // constructor format for all derived classes
  explicit
  IEMEvaluator(Teuchos::ParameterList& plist);
  IEMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEM>& iem);
  IEMEvaluator(const IEMEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& results);

  Teuchos::RCP<IEM> get_IEM() { return iem_; }

 protected:
  void InitializeFromPlist_();

  Key temp_key_;
  Teuchos::RCP<IEM> iem_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IEMEvaluator> factory_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
