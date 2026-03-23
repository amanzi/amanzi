/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  Field evaluator for the total internal energy:

    IE = phi * h_liquid * n_liquid + (1 - phi) * rho_rock * u_rock
*/

#include "CommonDefs.hh"

#include "TotalEnergyEvaluatorPH.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor from ParameterList
****************************************************************** */
TotalEnergyEvaluatorPH::TotalEnergyEvaluatorPH(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("energy key"), Tags::DEFAULT));
  }
  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);

  enthalpy_key_ = prefix + "enthalpy";
  mol_density_liquid_key_ = prefix + "molar_density_liquid";

  ie_rock_key_ = prefix + "internal_energy_rock";
  particle_density_key_ = prefix + "particle_density";
  porosity_key_ = prefix + "porosity";

  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));

  dependencies_.insert(std::make_pair(ie_rock_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(particle_density_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TotalEnergyEvaluatorPH::TotalEnergyEvaluatorPH(const TotalEnergyEvaluatorPH& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    enthalpy_key_(other.enthalpy_key_),
    mol_density_liquid_key_(other.mol_density_liquid_key_),
    ie_rock_key_(other.ie_rock_key_),
    particle_density_key_(other.particle_density_key_),
    porosity_key_(other.porosity_key_)
{}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
TotalEnergyEvaluatorPH::Clone() const
{
  return Teuchos::rcp(new TotalEnergyEvaluatorPH(*this));
}


/* ******************************************************************
* Field evaluator.
****************************************************************** */
void
TotalEnergyEvaluatorPH::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& h = *S.Get<CompositeVector>(enthalpy_key_).ViewComponent("cell");
  const auto& eta = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& u_rock = *S.Get<CompositeVector>(ie_rock_key_).ViewComponent("cell");
  const auto& rho_rock = *S.Get<CompositeVector>(particle_density_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = phi[0][c] * eta[0][c] * h[0][c] +
                     (1.0 - phi[0][c]) * u_rock[0][c] * rho_rock[0][c];
  }
}


/* ******************************************************************
* Field derivative evaluator.
****************************************************************** */
void
TotalEnergyEvaluatorPH::EvaluatePartialDerivative_(const State& S,
                                                   const Key& wrt_key,
                                                   const Tag& wrt_tag,
                                                   const std::vector<CompositeVector*>& results)
{
  const auto& h = *S.Get<CompositeVector>(enthalpy_key_).ViewComponent("cell");
  const auto& eta = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& u_rock = *S.Get<CompositeVector>(ie_rock_key_).ViewComponent("cell");
  const auto& rho_rock = *S.Get<CompositeVector>(particle_density_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = eta[0][c] * h[0][c] - rho_rock[0][c] * u_rock[0][c];
    }
  } else if (wrt_key == mol_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * h[0][c];
    }
  } else if (wrt_key == enthalpy_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * eta[0][c];
    }
  } else if (wrt_key == ie_rock_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = (1.0 - phi[0][c]) * rho_rock[0][c];
    }
  } else if (wrt_key == particle_density_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = (1.0 - phi[0][c]) * u_rock[0][c];
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace Energy
} // namespace Amanzi
