/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  The internal energu model evaluator simply calls the IEM
  with the correct arguments.
*/

#ifndef AMANZI_ENERGY_IEM_WATER_VAPOR_EVALUATOR_HH_
#define AMANZI_ENERGY_IEM_WATER_VAPOR_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "IEM_WaterVapor.hh"

namespace Amanzi {
namespace Energy {

class IEM_WaterVaporEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  explicit IEM_WaterVaporEvaluator(Teuchos::ParameterList& plist);
  IEM_WaterVaporEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEM_WaterVapor>& iem);
  IEM_WaterVaporEvaluator(const IEM_WaterVaporEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  Teuchos::RCP<IEM_WaterVapor> iem() { return iem_; }

 protected:
  void InitializeFromPlist_();

  Key temp_key_;
  Key mol_frac_key_;
  Teuchos::RCP<IEM_WaterVapor> iem_;

 private:
  static Utils::RegisteredFactory<Evaluator, IEM_WaterVaporEvaluator> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
