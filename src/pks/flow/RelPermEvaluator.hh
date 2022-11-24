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

#include "EvaluatorSecondaryMonotype.hh"
#include "Tag.hh"

#include "WRM.hh"
#include "WRMPartition.hh"

namespace Amanzi {
namespace Flow {

class RelPermEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  RelPermEvaluator(Teuchos::ParameterList& plist,
                   const Teuchos::Ptr<State>& S,
                   const Teuchos::RCP<WRMPartition>& wrm);
  RelPermEvaluator(const RelPermEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  void InitializeFromPlist_(const Teuchos::Ptr<State>& S);

 protected:
  Teuchos::RCP<WRMPartition> wrm_;
  Key pressure_key_;

  double patm_;
};

} // namespace Flow
} // namespace Amanzi

#endif
