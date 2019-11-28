/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Author: Ahmad Jan (jana@ornl.gov)
*/

//! Evaluates potential evapotranpiration (PET) 
//! Models are provide in the PERM-IV, Version 4, see pages 90-93, Equations 1-57 to 1-60
/*!

Requires the following dependencies:

* `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
* `"relative humicity key`" ``[string]`` **DOMAIN-relative_humidity**
* `"elevation key`" ``[string]`` **DOMAIN-elevation**
* `"shortwave radiation key`" ``[string]`` **DOMAIN-shortwave_radiation**

*/

#include "Key.hh"
#include "pet_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


PETEvaluator::PETEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  pt_alpha_ = plist.get<double>("adjustment factor alpha [-]",1.26);
  auto domain = Keys::getDomain(my_key_);
  at_key_ = Keys::readKey(plist, domain, "air temperature", "air_temperature");
  dependencies_.insert(at_key_);
  at_inter_key_ = Keys::readKey(plist, domain, "air temperature inter", "air_temperature_inter");
  dependencies_.insert(at_inter_key_);
  rel_hum_key_ = Keys::readKey(plist, domain, "relative humidity", "relative_humidity");
  dependencies_.insert(rel_hum_key_);
  elev_key_ = Keys::readKey(plist, domain, "elevation", "elevation");
  dependencies_.insert(elev_key_);
  swr_key_ = Keys::readKey(plist, domain, "shortwave radiation", "shortwave_radiation");
  dependencies_.insert(swr_key_);
  
}

// Required methods from SecondaryVariableFieldEvaluator
void
PETEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& air_temp = *S->GetFieldData(at_key_)->ViewComponent("cell", false);
  const auto& air_temp_inter = *S->GetFieldData(at_inter_key_)->ViewComponent("cell", false);
  const auto& rel_hum = *S->GetFieldData(rel_hum_key_)->ViewComponent("cell", false);
  const auto& elev = *S->GetFieldData(elev_key_)->ViewComponent("cell", false);
  const auto& swr = *S->GetFieldData(swr_key_)->ViewComponent("cell",false);

  auto& res = *result->ViewComponent("cell", false);

  for (int c=0; c!=res.MyLength(); ++c) {
    double lh_vap = LatentHeatVaporization(air_temp[0][c]);
    double ps_const = PsychrometricConstant(lh_vap, elev[0][c]);
    double vp_slope = VPSlope(air_temp[0][c]);

    double hf_ground = HeatFluxDensity(air_temp_inter[0][c], air_temp[0][c]);

    double s1 = vp_slope / (vp_slope + ps_const);
    double sw_rad = 2.064 * swr[0][c];
    double s2 = sw_rad / 23.88 - hf_ground;

    res[0][c] = 0.24* pt_alpha_ * (1./lh_vap) * s1 * s2 / 86400.; // convert mm to m, and per day to per second
    res[0][c] = std::max(res[0][c],0.0);
  }

}

double  
PETEvaluator::LatentHeatVaporization(double temp) 
{
  // convert temperature to Fahrenheit
  double temp_f = 1.8 * (temp - 273.15) + 32;
  return 597.3 - (0.5653 * temp_f);
}

double  
PETEvaluator::PsychrometricConstant(double lh_vap, double elev) 
{
  // convert elevation [m] to elevation [ft]
  double elev_ft = elev * 3.281;
  return 1.6286 * (101.3 - (0.003215 * elev_ft)) / lh_vap;
}

double
PETEvaluator::VPSlope(double temp) 
{
  // temperature conversion from K to C
  double temp_c = (temp - 273.15);

  double x = 17.26939*temp_c / (temp_c + 237.3);
  return 4098 * (0.6108 * std::exp(x)) / (std::pow(temp_c+237.3,2.0));
}

double
PETEvaluator::HeatFluxDensity(double temp_inter, double temp_next) 
{
  // temperature conversion from K to C
  double temp_prev = temp_inter - 273.15;
  double temp_curr = temp_next - 273.15;

  return -4.2 * (temp_prev - temp_curr);
}

} //namespace
} //namespace
} //namespace

