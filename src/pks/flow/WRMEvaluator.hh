/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The WRM Evaluator simply calls the WRM with the correct arguments.
*/

#ifndef AMANZI_FLOW_WRM_EVALUATOR_HH_
#define AMANZI_FLOW_WRM_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variables_field_evaluator.hh"
#include "WRM.hh"
#include "WRMPartition.hh"

namespace Amanzi {
namespace Flow {

class WRMEvaluator : public SecondaryVariablesFieldEvaluator {
 public:
  // constructor format for all derived classes
  explicit
  WRMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrm);
  WRMEvaluator(const WRMEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  void InitializeFromPlist_();

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<WRMPartition> wrm_;
  Key pressure_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator, WRMEvaluator> factory_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
