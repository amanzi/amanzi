/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Evaluators

  The porosity model evaluator simply calls the porosity model with
  the correct arguments.
*/

// Amanzi
#include "PDE_HelperDiscretization.hh"

// Amanzi::Evaluators
#include "Porosity_Compressible.hh"
#include "Porosity_Constant.hh"
#include "PorosityEvaluator.hh"
#include "PorosityModelPartition.hh"

namespace Amanzi {
namespace Evaluators {

/* ******************************************************************
* Two constructors.
****************************************************************** */
PorosityEvaluator::PorosityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
    initialized_(false),
    poroelasticity_(false)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("porosity key"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  poroelasticity_ = plist_.get<bool>("poroelasticity", false);
  thermoelasticity_ = plist_.get<bool>("thermoelasticity", false);

  // primary dependencre is pressure.
  pressure_key_ = plist_.get<std::string>("pressure key");
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));

  if (poroelasticity_) {
    biot_key_ = Keys::getKey(domain, "biot_coefficient");
    strain_key_ = plist_.get<std::string>("volumetric strain key");
    dependencies_.insert(std::make_pair(strain_key_, Tags::DEFAULT));
  }

  if (thermoelasticity_) {
    biot_key_ = Keys::getKey(domain, "biot_coefficient");
    temperature_key_ = plist_.get<std::string>("temperature key");
    dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));
  }
}


PorosityEvaluator::PorosityEvaluator(const PorosityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pom_(other.pom_),
    initialized_(other.initialized_),
    pressure_key_(other.pressure_key_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
PorosityEvaluator::Clone() const
{
  return Teuchos::rcp(new PorosityEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
PorosityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  if (!initialized_) {
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = S.GetMesh(Keys::getDomain(my_keys_.front().first));
    pom_ = CreatePorosityModelPartition(mesh, Teuchos::rcpFromRef(plist_.sublist("parameters")));
    initialized_ = true;
  }

  auto& phi_c = *results[0]->ViewComponent("cell", false);
  const auto& pres_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  int ncells = phi_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    phi_c[0][c] = pom_->second[(*pom_->first)[c]]->PorosityValue(pres_c[0][c]);
  }

  if (poroelasticity_) {
    const auto& e_c = *S.Get<CompositeVector>(strain_key_).ViewComponent("cell");
    const auto& b_c = *S.Get<CompositeVector>(biot_key_).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) {
      phi_c[0][c] += b_c[0][c] * e_c[0][c]; // e0 = 0.0
    }
  }

  if (thermoelasticity_) {
    const auto& temp_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");
    const auto& b_c = *S.Get<CompositeVector>(biot_key_).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) {
      double phi0 = pom_->second[(*pom_->first)[c]]->PorosityValueReference();
      double as = pom_->second[(*pom_->first)[c]]->getThermalCoefficients().second;
      phi_c[0][c] -= (b_c[0][c] - phi0) * as * (temp_c[0][c] - 273.15);
    }
  }

  // optional copy of cell data to boundary
  Operators::CellToBoundaryFaces(*results[0]);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
PorosityEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& results)
{
  auto& phi_c = *results[0]->ViewComponent("cell");
  const auto& pres_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  int ncells = phi_c.MyLength();
  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      phi_c[0][c] = pom_->second[(*pom_->first)[c]]->dPorositydPressure(pres_c[0][c]);
    }
  } else if (wrt_key == strain_key_) {
    const auto& b_c = *S.Get<CompositeVector>(biot_key_).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) phi_c[0][c] = b_c[0][c];
  }

  Operators::CellToBoundaryFaces(*results[0]);
}

} // namespace Evaluators
} // namespace Amanzi
