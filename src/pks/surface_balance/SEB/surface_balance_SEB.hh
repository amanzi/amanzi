/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Adam Atchley, Satish Karra

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   0=(1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qc(Ts) + Qe(Ts)
   Qc = the energy delived to the subsurface
   The rate of water entering the surface skin occurs only when Ts > 0
   In which case Ts is set = 0 and the excess energy = Qm and the melt rate (Mr) is
   delivered to the surface skin via:
   Mr = Qm/(ROWw*Hf)
   ROWw = density of water
   Hf = latent heat of fusion
   The temperature of water in assumed to be 0 C

   In cases without snow the energy balance equations is:

   Qex = 0=(1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qe(Ts)
   Qex is the energy derived to the subsurface skin
   All water entering the surface skin is assumed to be precipitated
   or condensed on the surface and has a temperature of Air.
   ------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_SEB_HH_
#define PKS_ENERGY_SEB_HH_

#include "pk_factory.hh"
#include "pk_physical_bdf_base.hh"
#include "primary_variable_field_evaluator.hh"


namespace Amanzi {
namespace SurfaceBalance {

class SurfaceBalanceSEB : public PKPhysicalBase {

 public:

  SurfaceBalanceSEB(Teuchos::Ptr<State> S, const Teuchos::RCP<Teuchos::ParameterList>& plist,
                    Teuchos::ParameterList& FElist,
                    const Teuchos::RCP<TreeVector>& solution);

  // SurfaceBalanceSEB is a PK
  // -- Setup data
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- provide a timestep size
  virtual double get_dt() {
    return dt_;
  }

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // -- advance via one of a few methods
  virtual bool advance(double dt);

 protected:
  // A few options for advance
  void CalculateSEB_();

 protected:

  // multiple primary variables
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_esource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_wsource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_wtemp_;

  double dt_;
  double min_wind_speed_;
  double albedo_trans_;
  double snow_ground_trans_;
  double no_snow_trans_;

 private:
  // factory registration
  static RegisteredPKFactory<SurfaceBalanceSEB> reg_;
};

} // namespace
} // namespace

#endif
