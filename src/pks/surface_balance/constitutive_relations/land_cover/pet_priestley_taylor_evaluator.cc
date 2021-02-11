/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

#include "Key.hh"
#include "pet_priestley_taylor_evaluator.hh"

namespace Amanzi {
namespace LandCover {
namespace Relations {

namespace PriestleyTaylor {

double
latentHeatVaporization_water(double temp_air)
{
  // convert temperature to Fahrenheit
  double temp_f = 1.8 * (temp_air - 273.15) + 32;
  // note 1/4184 converts cal/gm to J/kg
  return 597.3 - (0.5653 * temp_f) / 4184;
}

double
latentHeatVaporization_snow(double temp_air)
{
  // convert temperature to Fahrenheit
  double temp_f = 1.8 * (temp_air - 273.15) + 32;
  // note 1/4184 converts cal/gm to J/kg
  return 597.3 - (0.5653 * temp_f) / 4184;
}


double
psychrometricConstant(double lh_vap, double elev)
{
  // convert elevation [m] to elevation [ft]
  double elev_ft = elev * 3.281;
  return 1.6286 * (101.3 - (0.003215 * elev_ft)) / lh_vap;
}

double
vaporPressureSlope(double temp_air)
{
  // temperature conversion from K to C
  double temp_c = (temp_air - 273.15);
  double x = 17.26939*temp_c / (temp_c + 237.3);
  return 4098 * (0.6108 * std::exp(x)) / std::pow(temp_c+237.3,2);
}

double
groundHeatFlux(double temp_ground, double temp_air)
{
  return -4.2 * (temp_ground - temp_air);
}

} // namespace PriestleyTaylor


PETPriestleyTaylorEvaluator::PETPriestleyTaylorEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  pt_alpha_ = plist.get<double>("adjustment factor alpha [-]",1.26);

  auto domain = Keys::getDomain(my_key_);

  air_temp_key_ = Keys::readKey(plist, domain, "air temperature", "air_temperature");
  dependencies_.insert(air_temp_key_);

  ground_temp_ = Keys::readKey(plist, domain, "ground temperature", "temperature");
  dependencies_.insert(ground_temp_);

  rel_hum_key_ = Keys::readKey(plist, domain, "relative humidity", "relative_humidity");
  dependencies_.insert(rel_hum_key_);

  elev_key_ = Keys::readKey(plist, domain, "elevation", "elevation");
  dependencies_.insert(elev_key_);

  rad_key_ = Keys::readKey(plist, domain, "net radiation", "net_radiation");
  dependencies_.insert(rad_key_);

  limiter_ = plist.get<bool>("include limiter", false);
  if (limiter_) {
    limiter_key_ = Keys::readKey(plist, domain, "limiter");
    dependencies_.insert(limiter_key_);
  }

  one_minus_limiter_ = plist.get<bool>("include 1 - limiter", false);
  if (one_minus_limiter_) {
    one_minus_limiter_key_ = Keys::readKey(plist, domain, "1 - limiter");
    dependencies_.insert(one_minus_limiter_key_);
  }

  is_snow_ = plist.get<bool>("sublimate snow", false);
}

// Required methods from SecondaryVariableFieldEvaluator
void
PETPriestleyTaylorEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& air_temp = *S->GetFieldData(air_temp_key_)->ViewComponent("cell", false);
  const auto& air_temp_inter = *S->GetFieldData(ground_temp_)->ViewComponent("cell", false);
  const auto& rel_hum = *S->GetFieldData(rel_hum_key_)->ViewComponent("cell", false);
  const auto& elev = *S->GetFieldData(elev_key_)->ViewComponent("cell", false);
  const auto& rad = *S->GetFieldData(rad_key_)->ViewComponent("cell",false);

  auto& res = *result->ViewComponent("cell", false);

  for (int c=0; c!=res.MyLength(); ++c) {
    double lh_vap;
    if (is_snow_)
      lh_vap = PriestleyTaylor::latentHeatVaporization_snow(air_temp[0][c]);
    else
      lh_vap = PriestleyTaylor::latentHeatVaporization_water(air_temp[0][c]);
    double ps_const = PriestleyTaylor::psychrometricConstant(lh_vap, elev[0][c]);
    double vp_slope = PriestleyTaylor::vaporPressureSlope(air_temp[0][c]);
    double hf_ground = PriestleyTaylor::groundHeatFlux(air_temp_inter[0][c], air_temp[0][c]);

    double s1 = vp_slope / (vp_slope + ps_const);
    double s2 = rad[0][c] * 86400/1e6 - hf_ground; // converts net radiation in [W m^-2] to [MJ m^-2 d^-1]

    // note: conversions from mm to m (1e-3), MJ to J (1e6), and per day to per
    // second (/86400)
    double unit_conversion = 1e-3 * 1e6 / 86400.;
    res[0][c] = unit_conversion * pt_alpha_ / lh_vap * s1 * s2;
    res[0][c] = std::max(res[0][c],0.0);
  }

  // apply a limiter if requested
  if (limiter_) {
    const auto& limiter = *S->GetFieldData(limiter_key_)->ViewComponent("cell", false);
    res.Multiply(1, limiter, res, 0);
  }
  if (one_minus_limiter_) {
    const auto& one_minus_limiter = *S->GetFieldData(one_minus_limiter_key_)->ViewComponent("cell", false);
    res.Multiply(-1, limiter, res, 1);
  }
}


} //namespace
} //namespace
} //namespace

