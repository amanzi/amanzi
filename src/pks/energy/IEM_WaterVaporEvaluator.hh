/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  The internal energu model evaluator simply calls the IEM 
  with the correct arguments.
*/

#ifndef AMANZI_ENERGY_IEM_WATER_VAPOR_EVALUATOR_HH_
#define AMANZI_ENERGY_IEM_WATER_VAPOR_EVALUATOR_HH_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "IEM_WaterVapor.hh"

namespace Amanzi {
namespace Energy {

class IEM_WaterVaporEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  // constructor format for all derived classes
  explicit
  IEM_WaterVaporEvaluator(Teuchos::ParameterList& plist);
  IEM_WaterVaporEvaluator(Teuchos::ParameterList& plist,
                          const Teuchos::RCP<IEM_WaterVapor>& iem);
  IEM_WaterVaporEvaluator(const IEM_WaterVaporEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& results);

  Teuchos::RCP<IEM_WaterVapor> iem() { return iem_; }

 protected:
  void InitializeFromPlist_();

  Key temp_key_;
  Key mol_frac_key_;
  Teuchos::RCP<IEM_WaterVapor> iem_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IEM_WaterVaporEvaluator> factory_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
