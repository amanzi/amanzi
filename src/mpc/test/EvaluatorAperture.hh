/*
  MpC

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_APERTURE_EVALUATOR_HH_
#define AMANZI_FLOW_APERTURE_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"
#include "Tag.hh"

namespace Amanzi {

class EvaluatorAperture : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  EvaluatorAperture(Teuchos::ParameterList& plist);
  EvaluatorAperture(const EvaluatorAperture& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
   std::string pressure_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorAperture> reg_;
};


/* ******************************************************************
* Constructor.
****************************************************************** */
EvaluatorAperture::EvaluatorAperture(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("aperture key"), Tags::DEFAULT));
  }
  pressure_key_ = plist_.get<std::string>("pressure key");
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
EvaluatorAperture::EvaluatorAperture(const EvaluatorAperture& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_) {};


Teuchos::RCP<Evaluator> EvaluatorAperture::Clone() const {
  return Teuchos::rcp(new EvaluatorAperture(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void EvaluatorAperture::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = 1e-5 + (p_c[0][c] - 11e+6) / 1e+11;
  } 
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void EvaluatorAperture::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = 1e-11;
    }
  }
}

}  // namespace Amanzi

#endif
