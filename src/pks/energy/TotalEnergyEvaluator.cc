/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for the total internal energy:

    IE = phi * s_liquid * n_liquid * u_liquid 
       + phi * s_gas * n_gas * u_gas
       + (1 - phi) * rho_rock * u_rock
*/

#include "TotalEnergyEvaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor from ParameterList
****************************************************************** */
TotalEnergyEvaluator::TotalEnergyEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), aperture_(false)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("energy key"), Tags::DEFAULT));
  }
  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);

  vapor_diffusion_ = plist_.get<bool>("vapor diffusion");
  ie_rock_key_ = plist_.get<std::string>("internal energy rock key");
  ie_liquid_key_ = prefix + "internal_energy_liquid";
  ie_gas_key_ = prefix + "internal_energy_gas";

  mol_density_liquid_key_ = prefix + "molar_density_liquid";
  mol_density_gas_key_ = prefix + "molar_density_gas";

  particle_density_key_ = plist_.get<std::string>("particle density key");
  porosity_key_ = prefix + "porosity";
  sat_liquid_key_ = prefix + "saturation_liquid";

  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(sat_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(ie_liquid_key_, Tags::DEFAULT));

  if (vapor_diffusion_) {
    dependencies_.insert(std::make_pair(mol_density_gas_key_, Tags::DEFAULT));
    dependencies_.insert(std::make_pair(ie_gas_key_, Tags::DEFAULT));
  }

  dependencies_.insert(std::make_pair(ie_rock_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(particle_density_key_, Tags::DEFAULT));

  if (plist_.isParameter("aperture key")) {
    aperture_key_ = plist_.get<std::string>("aperture key");
    aperture_ = true;
  }
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TotalEnergyEvaluator::TotalEnergyEvaluator(const TotalEnergyEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    vapor_diffusion_(other.vapor_diffusion_),
    particle_density_key_(other.particle_density_key_),
    porosity_key_(other.porosity_key_),
    sat_liquid_key_(other.sat_liquid_key_),
    ie_rock_key_(other.ie_rock_key_),
    ie_liquid_key_(other.ie_liquid_key_),
    ie_gas_key_(other.ie_gas_key_),
    mol_density_liquid_key_(other.mol_density_liquid_key_),
    mol_density_gas_key_(other.mol_density_gas_key_){};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
TotalEnergyEvaluator::Clone() const
{
  return Teuchos::rcp(new TotalEnergyEvaluator(*this));
}


/* ******************************************************************
* Field evaluator.
****************************************************************** */
void
TotalEnergyEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& s_l = *S.Get<CompositeVector>(sat_liquid_key_).ViewComponent("cell");
  const auto& n_l = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& u_l = *S.Get<CompositeVector>(ie_liquid_key_).ViewComponent("cell");

  Teuchos::RCP<const Epetra_MultiVector> n_g, u_g;
  if (vapor_diffusion_) {
    n_g = S.Get<CompositeVector>(mol_density_gas_key_).ViewComponent("cell");
    u_g = S.Get<CompositeVector>(ie_gas_key_).ViewComponent("cell");
  }

  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& u_rock = *S.Get<CompositeVector>(ie_rock_key_).ViewComponent("cell");
  const auto& rho_rock = *S.Get<CompositeVector>(particle_density_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c] * u_l[0][c] +
                     (1.0 - phi[0][c]) * u_rock[0][c] * rho_rock[0][c];
    if (vapor_diffusion_) {
      double s_g = 1.0 - s_l[0][c];
      result_v[0][c] += phi[0][c] * s_g * (*n_g)[0][c] * (*u_g)[0][c];
    }
  }

  if (aperture_) {
    const auto& aperture = *S.Get<CompositeVector>(aperture_key_).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) { result_v[0][c] *= aperture[0][c]; }
  }
}


/* ******************************************************************
* Field derivative evaluator.
****************************************************************** */
void
TotalEnergyEvaluator::EvaluatePartialDerivative_(const State& S,
                                                 const Key& wrt_key,
                                                 const Tag& wrt_tag,
                                                 const std::vector<CompositeVector*>& results)
{
  const auto& s_l = *S.Get<CompositeVector>(sat_liquid_key_).ViewComponent("cell");
  const auto& n_l = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& u_l = *S.Get<CompositeVector>(ie_liquid_key_).ViewComponent("cell");

  Teuchos::RCP<const Epetra_MultiVector> n_g, u_g;
  if (vapor_diffusion_) {
    n_g = S.Get<CompositeVector>(mol_density_gas_key_).ViewComponent("cell");
    u_g = S.Get<CompositeVector>(ie_gas_key_).ViewComponent("cell");
  }

  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& u_rock = *S.Get<CompositeVector>(ie_rock_key_).ViewComponent("cell");
  const auto& rho_rock = *S.Get<CompositeVector>(particle_density_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = s_l[0][c] * n_l[0][c] * u_l[0][c] - rho_rock[0][c] * u_rock[0][c];
      if (vapor_diffusion_) {
        double s_g = 1.0 - s_l[0][c];
        result_v[0][c] += s_g * (*n_g)[0][c] * (*u_g)[0][c];
      }
    }
  } else if (wrt_key == sat_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_l[0][c] * u_l[0][c];
      if (vapor_diffusion_) { result_v[0][c] -= (*n_g)[0][c] * (*u_g)[0][c]; }
    }
  } else if (wrt_key == mol_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) { result_v[0][c] = phi[0][c] * s_l[0][c] * u_l[0][c]; }
  } else if (wrt_key == ie_liquid_key_) {
    for (int c = 0; c != ncells; ++c) { result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c]; }

  } else if (wrt_key == mol_density_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      double s_g = 1.0 - s_l[0][c];
      result_v[0][c] = phi[0][c] * s_g * (*u_g)[0][c];
    }
  } else if (wrt_key == ie_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      double s_g = 1.0 - s_l[0][c];
      result_v[0][c] = phi[0][c] * s_g * (*n_g)[0][c];
    }

  } else if (wrt_key == ie_rock_key_) {
    for (int c = 0; c != ncells; ++c) { result_v[0][c] = (1.0 - phi[0][c]) * rho_rock[0][c]; }
  } else if (wrt_key == particle_density_key_) {
    for (int c = 0; c != ncells; ++c) { result_v[0][c] = (1.0 - phi[0][c]) * u_rock[0][c]; }
  } else {
    AMANZI_ASSERT(0);
  }

  if (aperture_) {
    const auto& aperture = *S.Get<CompositeVector>(aperture_key_).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) { result_v[0][c] *= aperture[0][c]; }
  }
}

} // namespace Energy
} // namespace Amanzi
