/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The fracture permeability evaluator simply calls the permability 
  model with the correct arguments.
*/

#ifndef AMANZI_FLOW_POROSITY_EVALUATOR_HH_
#define AMANZI_FLOW_POROSITY_EVALUATOR_HH_

#include "factory.hh"
#include "secondary_variables_field_evaluator.hh"

#include "FracturePermModel.hh"
#include "FracturePermModelPartition.hh"

namespace Amanzi {
namespace Flow {

class FracturePermModelEvaluator : public SecondaryVariablesFieldEvaluator {
 public:
  // constructor format for all derived classes
  explicit
  FracturePermModelEvaluator(Teuchos::ParameterList& plist,
                             Teuchos::RCP<FracturePermModelPartition> fpm);
  FracturePermModelEvaluator(const FracturePermModelEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<FracturePermModelPartition> fpm_;
  Key my_key_;
  Key aperture_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator, FracturePermModelEvaluator> factory_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
