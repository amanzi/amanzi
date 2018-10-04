/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Determine the volumetric ponded and snow depths from ponded depth and snow depth.
/*!

* `"maximum relief key`" ``[string]`` **DOMAIN-maximum_relief**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"ponded depth key`" ``[string]`` **DOMAIN-ponded_depth**
         The true height of the water surface.
* `"snow depth key`" ``[string]`` **DOMAIN-snow_depth**
         The true height of the snow surface.

*/

#include "volumetric_height_subgrid_evaluator.hh"

namespace Amanzi {
namespace Flow {


VolumetricHeightSubgridEvaluator::VolumetricHeightSubgridEvaluator(Teuchos::ParameterList& plist) :
     SecondaryVariablesFieldEvaluator(plist)
{
  Key a_key = Keys::cleanPListName(plist.name());
  Key domain = Keys::getDomain(a_key);
  Key domain_surf, domain_snow;
  if (domain == "surface" || domain == "snow") {
    domain_surf = "surface";
    domain_snow = "snow";
  } else if (boost::starts_with(domain, "surface_")) {
    domain_surf = domain;
    domain_snow = std::string("snow_") + domain.substr(8,domain.size());
  } else if (boost::starts_with(domain, "snow_")) {
    domain_snow = domain;
    domain_surf = std::string("surface_") + domain.substr(5,domain.size());
  } else {
    Errors::Message msg("VolumetricHeightSubgridEvaluator: not sure how to interpret domain.");
    Exceptions::amanzi_throw(msg);
  }

  // my keys
  vol_pd_key_ = Keys::readKey(plist, domain_surf, "volumetric ponded depth", "volumetric_ponded_depth");
  my_keys_.push_back(vol_pd_key_);
  vol_sd_key_ = Keys::readKey(plist, domain_snow, "volumetric snow depth", "volumetric_snow_depth");
  my_keys_.push_back(vol_sd_key_);
  
  // dependencies
  pd_key_ = Keys::readKey(plist_, domain, "ponded depth key", "ponded_depth");
  dependencies_.insert(pd_key_);
  sd_key_ = Keys::readKey(plist_, domain, "snow depth key", "snow_depth");
  dependencies_.insert(sd_key_);

  delta_max_key_ = Keys::readKey(plist_, domain, "microtopographic relief", "microtopographic_relief"); 
  dependencies_.insert(delta_max_key_);
  
  delta_ex_key_ = Keys::readKey(plist_, domain, "excluded volume", "excluded_volume");
  dependencies_.insert(delta_ex_key_);
}


void
VolumetricHeightSubgridEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                             const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  auto& vpd = *results[0]->ViewComponent("cell",false);
  auto& vsd = *results[1]->ViewComponent("cell",false);
  const auto& pd = *S->GetFieldData(pd_key_)->ViewComponent("cell",false);
  const auto& sd = *S->GetFieldData(sd_key_)->ViewComponent("cell",false);
  const auto& del_max = *S->GetFieldData(delta_max_key_)->ViewComponent("cell",false);
  const auto& del_ex = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell",false);

  for (int c=0; c!=vpd.MyLength(); ++c){
    double vol_tot = f_(pd[0][c] + sd[0][c], del_max[0][c], del_ex[0][c]);
    vpd[0][c] = f_(pd[0][c], del_max[0][c], del_ex[0][c]);
    vsd[0][c] = vol_tot - vpd[0][c];
  }
}


void
VolumetricHeightSubgridEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  auto& vpd = *results[0]->ViewComponent("cell",false);
  auto& vsd = *results[1]->ViewComponent("cell",false);
  const auto& pd = *S->GetFieldData(pd_key_)->ViewComponent("cell",false);
  const auto& sd = *S->GetFieldData(sd_key_)->ViewComponent("cell",false);
  const auto& del_max = *S->GetFieldData(delta_max_key_)->ViewComponent("cell",false);
  const auto& del_ex = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell",false);

  if (wrt_key == pd_key_) {
    for (int c=0; c!=vpd.MyLength(); ++c){
      vpd[0][c] = f_prime_(pd[0][c], del_max[0][c], del_ex[0][c]);
      vsd[0][c] = f_prime_(pd[0][c] + sd[0][c], del_max[0][c], del_ex[0][c]);
    }
  } else {
    Errors::Message msg("VolumetricHeightSubgridEvaluator: Not Implemented: no derivatives implemented other than ponded depth.");
    Exceptions::amanzi_throw(msg);
  }
}

} //namespace
} //namespace
