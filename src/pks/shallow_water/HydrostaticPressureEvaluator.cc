/*
  Shallow Water PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include "errors.hh"

#include "HydrostaticPressureEvaluator.hh"

namespace Amanzi {
namespace ShallowWater {

/* ******************************************************************
* Simple constructor
****************************************************************** */
HydrostaticPressureEvaluator::HydrostaticPressureEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key", "surface-ponded_pressure"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  shallow_water_model_ = plist.get<int>("use shallow water model", 1);

  if (shallow_water_model_){
     primary_variable_key_ = plist_.get<std::string>("ponded depth key", Keys::getKey(domain, "ponded_depth"));
  } else {
     primary_variable_key_ = plist_.get<std::string>("water depth key", Keys::getKey(domain, "water_depth"));
  }

  dependencies_.insert(std::make_pair(primary_variable_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator> HydrostaticPressureEvaluator::Clone() const {
  return Teuchos::rcp(new HydrostaticPressureEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void HydrostaticPressureEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& h_c = *S.Get<CompositeVector>(primary_variable_key_).ViewComponent("cell");
  const double rho = S.Get<double>("const_fluid_density");
  const double patm = S.Get<double>("atmospheric_pressure");
  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));
  
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = rho * g * h_c[0][c] + patm;
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void HydrostaticPressureEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const double rho = S.Get<double>("const_fluid_density");
  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));
  
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == primary_variable_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = rho * g;
    }
  } else {
    AMANZI_ASSERT(false);
  }
}

}  // namespace ShallowWater
}  // namespace Amanzi

