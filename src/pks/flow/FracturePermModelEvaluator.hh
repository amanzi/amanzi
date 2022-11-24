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

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Tag.hh"

#include "FracturePermModel.hh"
#include "FracturePermModelPartition.hh"

namespace Amanzi {
namespace Flow {

class FracturePermModelEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  explicit FracturePermModelEvaluator(Teuchos::ParameterList& plist,
                                      Teuchos::RCP<FracturePermModelPartition> fpm);
  FracturePermModelEvaluator(const FracturePermModelEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<FracturePermModelPartition> fpm_;
  Key aperture_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, FracturePermModelEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
