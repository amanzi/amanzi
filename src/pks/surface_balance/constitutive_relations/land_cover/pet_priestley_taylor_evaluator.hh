/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (coonet@ornl.gov)
*/
//! Evaluates potential evapotranpiration (PET) using Priestley & Taylor formulation.
/*!

This implementation is based on models provided in the PRMS-IV, Version 4, see
pages 90-93, Equations 1-57 to 1-60

Requires the following dependencies:

.. _pet-priestley-taylor-evaluator-spec:
.. admonition:: pet-priestley-taylor-evaluator-spec:

   * `"include limiter`" ``[bool]`` If true, multiply potential ET by a limiter
     to get an actual ET.
   * `"include 1 - limiter`" ``[bool]`` If true, multiply potential ET by
     1 - a limiter (e.g. a limiter that partitions between two pools) to get
     actual ET.
   * `"sublimate snow`" ``[bool]`` **false** If true, use latent heat of
      vaporization of snow, not water.

   KEYS:
   * `"air temperature`" Air temp, in [K]
   * `"ground temperature`" Ground or leaf temp, in [K].  Note this may be the
      same as air temperature.
   * `"relative humidity`" [-]
   * `"elevation`" Elevation [m]
   * `"net radiation`" [W m^-2] Net radiation balance, positive to the ground.
   * `"limiter`" [-] See `"include limiter`" above.
   * `"1 - limiter`" [-] See `"include 1 - limiter`" above.

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

namespace PriestleyTaylor {

//
// PRMS-IV eqn 1-59, calculates ground heat flux density in units of
// [MJ m^-2 d^-1] given a daily-averaged ground and air temperature
// (in C or K).  We convert to W/m^2
//
double groundHeatFlux(double temp_ground, double temp_air);

//
// PRMS-IV eqn 1-58, calculates the slope of vapor pressure as a function of
// daily averaged air temperature [K], in [KPa C^-1]
//
double vaporPressureSlope(double temp_air);

//
// PRMS-IV eqn 1-57, calculates the psychrometric constant in [KPa C^-1] as a
// function of an elevation (lapse rate fixed) and a latent heat of
// vaporization in [cal gm^-1].
//
double psychrometricConstant(double lh_vap, double evel);

//
// PRMS-IV eqn 1-51, calculates the latent heat of vaporization [cal g^-1] as a
// function of the daily averaged air temperature [K] for liquid water
//
double latentHeatVaporization_water(double temp_air);

//
// PRMS-IV eqn 1-51, calculates the latent heat of vaporization [cal g^-1] as a
// function of the daily averaged air temperature [K] for snow -- note this is
// currently the same as the water value, but should get modified for snow!
//
double latentHeatVaporization_snow(double temp_air);


} // namespace PriestleyTaylor


class PETPriestleyTaylorEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  PETPriestleyTaylorEvaluator(Teuchos::ParameterList& plist);
  PETPriestleyTaylorEvaluator(const PETPriestleyTaylorEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new PETPriestleyTaylorEvaluator(*this));
  }

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Exceptions::amanzi_throw("NotImplemented: PETPriestleyTaylorEvaluator currently does not provide derivatives.");
  }

 protected:

  Key air_temp_key_;
  Key ground_temp_key_;
  Key rel_hum_key_;
  Key elev_key_;
  Key rad_key_;
  Key limiter_key_;
  Key one_minus_limiter_key_;

  double pt_alpha_;
  bool limiter_;
  bool one_minus_limiter_;
  bool is_snow_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,PETPriestleyTaylorEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

