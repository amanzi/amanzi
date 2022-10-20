/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for water storage which is the conserved quantity 
  in Richards' equation.
*/

#include "CommonDefs.hh"
#include "WaterStorage.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Constructor.
****************************************************************** */
WaterStorage::WaterStorage(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
    aperture_(false) {
  Init_();
};


/* ******************************************************************
* Initialization.
****************************************************************** */
void WaterStorage::Init_()
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("water storage key"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  saturation_key_ = plist_.get<std::string>("saturation key");
  porosity_key_ = plist_.get<std::string>("porosity key");
  mol_density_liquid_key_ = Keys::getKey(domain, "molar_density_liquid");

  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(saturation_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_liquid_key_, Tags::DEFAULT));

  water_vapor_ = plist_.get<bool>("water vapor", false);

  if (water_vapor_) {
    dependencies_.insert(std::make_pair(Keys::getKey(domain, "molar_density_gas"), Tags::DEFAULT));
    dependencies_.insert(std::make_pair(Keys::getKey(domain, "molar_fraction_gas"), Tags::DEFAULT));
  }

  if (plist_.isParameter("aperture key")) {
    aperture_ = true;
    aperture_key_ = plist_.get<std::string>("aperture key");
    dependencies_.insert(std::make_pair(Keys::getKey(domain, "aperture"), Tags::DEFAULT));
  }
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
WaterStorage::WaterStorage(const WaterStorage& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    water_vapor_(other.water_vapor_),
    aperture_(other.aperture_) {};


Teuchos::RCP<Evaluator> WaterStorage::Clone() const {
  return Teuchos::rcp(new WaterStorage(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void WaterStorage::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& s_l = *S.Get<CompositeVector>(saturation_key_).ViewComponent("cell");
  const auto& n_l = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");

  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (water_vapor_) {
    const auto& n_g = *S.Get<CompositeVector>("molar_density_gas").ViewComponent("cell");
    const auto& x_g = *S.Get<CompositeVector>("molar_fraction_gas").ViewComponent("cell");
    
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * (s_l[0][c] * n_l[0][c]
                                  + (1.0 - s_l[0][c]) * n_g[0][c] * x_g[0][c]);
    }
  } else {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c];
    }
  } 
     
  if (aperture_) {
    const auto& aperture = *S.Get<CompositeVector>(aperture_key_).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] *= aperture[0][c];
    }
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void WaterStorage::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& s_l = *S.Get<CompositeVector>(saturation_key_).ViewComponent("cell");
  const auto& n_l = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");

  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (water_vapor_) {
    const auto& n_g = *S.Get<CompositeVector>("molar_density_gas").ViewComponent("cell");
    const auto& x_g = *S.Get<CompositeVector>("molar_fraction_gas").ViewComponent("cell");

    if (wrt_key == porosity_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = (s_l[0][c] * n_l[0][c] + (1.0 - s_l[0][c]) * n_g[0][c] * x_g[0][c]);
      }
    } else if (wrt_key == saturation_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * n_l[0][c];
      }
    } else if (wrt_key == mol_density_liquid_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * s_l[0][c];
      }
    } else if (wrt_key == "molar_density_gas") {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * (1.0 - s_l[0][c]) * x_g[0][c];
      }
    } else if (wrt_key == "molar_fraction_gas") {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * (1.0 - s_l[0][c]) * n_g[0][c];
      }
    } else if (wrt_key == aperture_key_ && aperture_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * (s_l[0][c] * n_l[0][c] + (1.0 - s_l[0][c]) * n_g[0][c] * x_g[0][c]);
      }
    } else {
      AMANZI_ASSERT(0);
    }

  } else {
    if (wrt_key == porosity_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = s_l[0][c] * n_l[0][c];
      }
    } else if (wrt_key == saturation_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * n_l[0][c];
      }
    } else if (wrt_key == mol_density_liquid_key_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * s_l[0][c];
      }
    } else if (wrt_key == aperture_key_ && aperture_) {
      for (int c = 0; c != ncells; ++c) {
        result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c];
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }

  if (aperture_ && wrt_key != aperture_key_) {
    const auto& aperture = *S.Get<CompositeVector>(aperture_key_).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] *= aperture[0][c];
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi
