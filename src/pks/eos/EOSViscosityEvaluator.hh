/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOSFieldEvaluator is the interface between state/data and the model, an EOS.
*/

#ifndef AMANZI_EOS_VISCOSITY_EVALUATOR_HH_
#define AMANZI_EOS_VISCOSITY_EVALUATOR_HH_

#include "EOS_Viscosity.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace AmanziEOS {

class EOSViscosityEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit EOSViscosityEvaluator(Teuchos::ParameterList& plist);
  EOSViscosityEvaluator(const EOSViscosityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  // the actual model
  Teuchos::RCP<EOS_Viscosity> visc_;

  // Keys for fields
  // dependencies
  Key temp_key_, pres_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, EOSViscosityEvaluator> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
