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
#include "EvaluatorSecondaryMonotype.hh"

// Multiphase
#include "MultiphaseTypeDefs.hh"
#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

class RelPermEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  RelPermEvaluator(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<WRMmpPartition>& wrm);
  RelPermEvaluator(const RelPermEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  Teuchos::RCP<WRMmpPartition> wrm_;
 
  int phase_;
  Key saturation_liquid_key_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
