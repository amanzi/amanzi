/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  The interface between state/data and the model, an EOS.
*/

#ifndef AMANZI_EOS_ISOBARIC_EOS_EVALUATOR_HH_
#define AMANZI_EOS_ISOBARIC_EOS_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

#include "EOS_Density.hh"

namespace Amanzi {
namespace AmanziEOS {

class IsobaricEOSEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  enum EOSMode { EOS_MODE_MASS, EOS_MODE_MOLAR, EOS_MODE_BOTH };

  explicit IsobaricEOSEvaluator(Teuchos::ParameterList& plist);
  IsobaricEOSEvaluator(const IsobaricEOSEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  Teuchos::RCP<EOS_Density> get_EOS() { return eos_; }

 protected:
  // the actual model
  Teuchos::RCP<EOS_Density> eos_;
  EOSMode mode_;

  Key temp_key_, pres_key_;
  Tag tag_;

 private:
  static Utils::RegisteredFactory<Evaluator, IsobaricEOSEvaluator> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
