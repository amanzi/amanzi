/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The porosity model evaluator simply calls the porosity model with 
  the correct arguments.
*/

#include "FlowDefs.hh"
#include "PorosityModel_Compressible.hh"
#include "PorosityModel_Constant.hh"
#include "PorosityModelEvaluator.hh"
#include "PorosityModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
PorosityModelEvaluator::PorosityModelEvaluator(
    Teuchos::ParameterList& plist, Teuchos::RCP<PorosityModelPartition> pom)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
      pom_(pom)
{
  InitializeFromPlist_();
}


PorosityModelEvaluator::PorosityModelEvaluator(const PorosityModelEvaluator& other)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
      pressure_key_(other.pressure_key_),
      pom_(other.pom_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator> PorosityModelEvaluator::Clone() const {
  return Teuchos::rcp(new PorosityModelEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void PorosityModelEvaluator::InitializeFromPlist_()
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("porosity key"), Tags::DEFAULT));
  }

  // my dependency is pressure.
  pressure_key_ = plist_.get<std::string>("pressure key");
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void PorosityModelEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  auto& phi_c = *results[0]->ViewComponent("cell", false);
  const auto& pres_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  int ncells = phi_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    phi_c[0][c] = pom_->second[(*pom_->first)[c]]->Porosity(pres_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void PorosityModelEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  auto& phi_c = *results[0]->ViewComponent("cell");
  const auto& pres_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  int ncells = phi_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    phi_c[0][c] = pom_->second[(*pom_->first)[c]]->dPorositydPressure(pres_c[0][c]);
  }
}

}  // namespace Flow
}  // namespace Amanzi
