/*
  SurfaceBalanceEvaluator evaluates SEB as a nonlinear source instead of as a PK.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Debugger.hh"
#include "surface_balance_evaluator.hh"
#include "SnowEnergyBalance.hh"

namespace Amanzi {
namespace SurfaceBalance {

SurfaceBalanceEvaluator::SurfaceBalanceEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist) {

  // my keys
  my_keys_.push_back("surface_mass_source");
  my_keys_.push_back("surface_conducted_energy_source");
  my_keys_.push_back("surface_mass_source_temperature");
  my_keys_.push_back("snow_depth");
  my_keys_.push_back("snow_density");
  my_keys_.push_back("days_no_snow");

  // parameters
  min_wind_speed_ = plist_.get<double>("minimum wind speed", 1.0);
  snow_ground_trans_ = plist_.get<double>("minimum snow depth", 0.02);
  albedo_trans_ = plist_.get<double>("albedo transition depth", 0.02);

  // dependencies
  dependencies_.insert("surface_temperature");
  dependencies_.insert("ponded_depth");
  dependencies_.insert("surface_porosity");
  dependencies_.insert("surface_vapor_pressure");
  dependencies_.insert("air_temperature");
  dependencies_.insert("incoming_shortwave_radiation");
  dependencies_.insert("relative_humidity");
  dependencies_.insert("wind_speed");
  dependencies_.insert("precipitation_rain");
  dependencies_.insert("precipitation_snow");

};


SurfaceBalanceEvaluator::SurfaceBalanceEvaluator(const SurfaceBalanceEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    min_wind_speed_(other.min_wind_speed_),
    snow_ground_trans_(other.snow_ground_trans_),
    albedo_trans_(other.albedo_trans_),
    db_(other.db_) {}


Teuchos::RCP<FieldEvaluator>
SurfaceBalanceEvaluator::Clone() const {
  return Teuchos::rcp(new SurfaceBalanceEvaluator(*this));
}


void
SurfaceBalanceEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  S->RequireScalar("T_prev", my_keys_[0]);
  SecondaryVariablesFieldEvaluator::EnsureCompatibility(S);

  // Some of the provided vectors are not used, so we'll Require them here
  // just in case.
  S->RequireField("surface_mass_source", "surface_mass_source")
      ->SetMesh(S->GetMesh("surface"))->AddComponent("cell",AmanziMesh::CELL,1);
  S->RequireField("surface_conducted_energy_source",
                  "surface_conducted_energy_source")
      ->SetMesh(S->GetMesh("surface"))->AddComponent("cell",AmanziMesh::CELL,1);
  S->RequireField("surface_mass_source_temperature",
                  "surface_mass_source_temperature")
      ->SetMesh(S->GetMesh("surface"))->AddComponent("cell",AmanziMesh::CELL,1);
  S->RequireField("snow_depth", "snow_depth")
      ->SetMesh(S->GetMesh("surface"))->AddComponent("cell",AmanziMesh::CELL,1);
  S->RequireField("snow_density", "snow_density")
      ->SetMesh(S->GetMesh("surface"))->AddComponent("cell",AmanziMesh::CELL,1);
  S->RequireField("days_no_snow", "days_no_snow")
      ->SetMesh(S->GetMesh("surface"))->AddComponent("cell",AmanziMesh::CELL,1);

  // And additionally require this evaluator, again because not all are
  // required.  They will simply be aliased in state.
  S->RequireFieldEvaluator("surface_mass_source");
  S->RequireFieldEvaluator("surface_mass_source_temperature");
  S->RequireFieldEvaluator("surface_conducted_energy_source");
  S->RequireFieldEvaluator("snow_depth");
  S->RequireFieldEvaluator("snow_density");
  S->RequireFieldEvaluator("days_no_snow");

}

void
SurfaceBalanceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  if (db_ == Teuchos::null) {
    // Debugger
    db_ = Teuchos::rcp(new Debugger(S->GetMesh("surface"), my_keys_[0], plist_));
  }

  // Pull dependencies out of state.
  double dt = S->time() - *S->GetScalarData("T_prev");
  *S->GetScalarData("T_prev", my_keys_[0]) = S->time();

  const Epetra_MultiVector& temp = *S->GetFieldData("surface_temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& pd = *S->GetFieldData("ponded_depth")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& poro = *S->GetFieldData("surface_porosity")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& vp = *S->GetFieldData("surface_vapor_pressure")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& air_temp = *S->GetFieldData("air_temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& rad = *S->GetFieldData("incoming_shortwave_radiation")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& humid = *S->GetFieldData("relative_humidity")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& wind = *S->GetFieldData("wind_speed")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& precip_rain = *S->GetFieldData("precipitation_rain")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& precip_snow = *S->GetFieldData("precipitation_snow")
      ->ViewComponent("cell",false);


  Epetra_MultiVector& Q_mass = *results[0]->ViewComponent("cell",false);
  Epetra_MultiVector& Q_energy = *results[1]->ViewComponent("cell",false);
  Epetra_MultiVector& T_Qmass = *results[2]->ViewComponent("cell",false);
  Epetra_MultiVector& snow_depth = *results[3]->ViewComponent("cell",false);
  Epetra_MultiVector& snow_dens = *results[4]->ViewComponent("cell",false);
  Epetra_MultiVector& days_no_snow = *results[5]->ViewComponent("cell",false);

  // Create the SEB data structure
  SurfaceEnergyBalance::LocalData data;
  // Populate the SEB data structure with constants
  data.st_energy.stephB = 0.00000005670373;// Stephan-boltzmann Constant ------- [W/m^2 K^4]
  data.st_energy.Hf = 333500.0;            // Heat of fusion for melting snow -- [J/kg]
  data.st_energy.Ls = 2834000.0;           // Latent heat of sublimation ------- [J/kg]
  data.st_energy.Le = 2497848.;            // Latent heat of vaporization ------ [J/kg]
  data.st_energy.SEs = 0.98;               // Surface Emissivity for snow  ----- [-] ** From P. ReVelle (Thesis)
  data.st_energy.SEtun = 0.92;             // Surface Emissivity for tundray --- [-] ** From P. ReVelle (Thesis)
  data.st_energy.Zr = 2.0;                 // Referance ht of wind speed ------- [m]
  data.st_energy.Zo = 0.005;               // Roughness length  ---------------- [m]
  data.st_energy.VKc = 0.41;               // Von Karman Constant -------------- [-]
  double Cp = 1004.0;                      // Specific heat of air ------------- [J/K kg]
  data.st_energy.Apa = 100;                // Atmospheric Pressure ------------- [KPa]

  data.st_energy.density_w = 1000;         // Density of Water ----------------- [kg/m^3]
  double density_air = 1.275;              // Density of Air ------------------- [kg/m^3]
  data.st_energy.density_frost = 200;      // Density of Frost (condensation) -- [kg/m^3]
  data.st_energy.density_freshsnow = 100;  // Density of Freshly fallebn snow -- [kg/m^3]

  data.st_energy.gZr = 9.807*data.st_energy.Zr;
  data.st_energy.rowaCp = density_air*Cp;
  data.st_energy.rowaLs = density_air*data.st_energy.Ls;
  data.st_energy.rowaLe = density_air*data.st_energy.Le;

  data.vp_snow.relative_humidity = 1.;
  data.vp_ground.relative_humidity = 1.;

  data.st_energy.Dt = dt;
  data.st_energy.AlbedoTrans = albedo_trans_;
  data.st_energy.snow_groundTrans = snow_ground_trans_;

  // Calculating which Day frost density matched snow Defermation fucntion from (Martinec, 1977)
  data.st_energy.NDSfrost = std::pow((data.st_energy.density_frost
          / data.st_energy.density_freshsnow), (1/0.3)) - 1;

  int count = Q_mass.MyLength();
  for (unsigned int c=0; c!=count; ++c) {
    // ATS Calcualted Data
    data.st_energy.Tb = temp[0][c];
    data.st_energy.water_depth = pd[0][c];
    data.st_energy.por = poro[0][c];
    data.vp_ground.temp = temp[0][c];
    data.vp_ground.actual_vaporpressure = vp[0][c];
    data.st_energy.porrowaLe = data.st_energy.por
        * density_air * data.st_energy.Le;

    // MET station data
    data.st_energy.air_temp = air_temp[0][c];
    data.st_energy.QswIn = rad[0][c];
    data.st_energy.Us = std::max(wind[0][c], min_wind_speed_);
    data.st_energy.Pr = precip_rain[0][c] * dt; // SEB expects total precip, not rate
    data.st_energy.Ps = precip_snow[0][c] * dt; // SEB expects total precip, not rate
    data.vp_air.temp = air_temp[0][c];
    data.vp_air.relative_humidity = humid[0][c];

    // STORED INFO FOR SnowEnergyBalanc Model
    data.st_energy.ht_snow = snow_depth[0][c];
    data.st_energy.density_snow = snow_dens[0][c];
    data.st_energy.nosnowdays = days_no_snow[0][c];

    // Snow-ground Smoothing 
    if ((data.st_energy.ht_snow > data.st_energy.snow_groundTrans)
        || (data.st_energy.ht_snow <= 0)) {
      // Run the Snow Energy Balance Model as normal.
      SurfaceEnergyBalance::SnowEnergyBalance(data);
    } else {
      // Transition between Snow and bare ground
      double theta=0.0;
      double QcSnow=0.0, QcGrnd=0.0, MrSnow=0.0, MrGrnd=0.0, TrwSnow=0.0, TrwGrnd=0.0;
      double ZsHold=0.0, ROWsHold=0.0, nosnowdayHold=0.0;

      // Fraction of snow covered ground
      theta = std::pow(data.st_energy.ht_snow / data.st_energy.snow_groundTrans, 2);
      SurfaceEnergyBalance::SnowEnergyBalance(data);

      // Snow covered data
      QcSnow = theta * data.st_energy.fQc;
      MrSnow = data.st_energy.Mr - (data.st_energy.Pr / data.st_energy.Dt) + theta*(data.st_energy.Pr / data.st_energy.Dt);
      TrwSnow = data.st_energy.Trw * MrSnow;

      // Saving Snow data for next timestep
      ZsHold = data.st_energy.ht_snow;
      ROWsHold = data.st_energy.density_snow;
      nosnowdayHold = data.st_energy.nosnowdays;

      // Running the Snow Energy Balance Model for groundcover portion
      data.st_energy.ht_snow = 0;
      SurfaceEnergyBalance::SnowEnergyBalance(data);

      // Ground Cover Data
      QcGrnd = (1-theta)*data.st_energy.fQc;
      MrGrnd = (1-theta)*data.st_energy.Mr;

      if (MrGrnd <= 0) {
        TrwGrnd = 0; // No water
      } else {
        TrwGrnd = data.st_energy.Trw * MrGrnd;
      }

      // Calculating Data for ATS
      data.st_energy.fQc =  QcSnow + QcGrnd;
      data.st_energy.Mr = MrSnow + MrGrnd;
      if ((MrSnow + MrGrnd) <= 0) {
        // No Water delivered to ATS, therefore temperature of water is
        // irrelevant.
        data.st_energy.Trw = 0.0;
      } else {
        // Temperature of water delivered to ATS is averaged by rate/mass of
        // water from snow and no snow areas
        data.st_energy.Trw = (TrwSnow + TrwGrnd) / (MrSnow + MrGrnd);
      }

      // Put saved snow data back in place.
      data.st_energy.ht_snow = ZsHold;
      data.st_energy.density_snow = ROWsHold;
      data.st_energy.nosnowdays = nosnowdayHold;
    }

    // Store results
    Q_energy[0][c] = data.st_energy.fQc;
    Q_mass[0][c] = data.st_energy.Mr;
    T_Qmass[0][c] = data.st_energy.Trw;
    snow_depth[0][c] = data.st_energy.ht_snow;
    snow_dens[0][c] = data.st_energy.density_snow;
    days_no_snow[0][c] = data.st_energy.nosnowdays;

    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Snow depth, snowtemp = " << data.st_energy.ht_snow
                 << ", " << data.st_energy.Ts << std::endl;
      *vo_->os() << "Latent heat = " << data.st_energy.fQe << std::endl;
      *vo_->os() << "Sensible heat = " << data.st_energy.fQh << std::endl;
    }
  }

  // debug
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("air_temp"); vecs.push_back(S->GetFieldData("air_temperature").ptr());
    vnames.push_back("air_temp2"); vecs.push_back(S->GetFieldData("air_temperature").ptr());
    vnames.push_back("Qsw_in"); vecs.push_back(S->GetFieldData("incoming_shortwave_radiation").ptr());
    vnames.push_back("precip_rain"); vecs.push_back(S->GetFieldData("precipitation_rain").ptr());
    vnames.push_back("precip_snow"); vecs.push_back(S->GetFieldData("precipitation_snow").ptr());
    vnames.push_back("soil vapor pressure"); vecs.push_back(S->GetFieldData("surface_vapor_pressure").ptr());
    vnames.push_back("T_ground"); vecs.push_back(S->GetFieldData("surface_temperature").ptr());
    vnames.push_back("water_source"); vecs.push_back(S->GetFieldData("surface_mass_source").ptr());
    vnames.push_back("e_source"); vecs.push_back(S->GetFieldData("surface_conducted_energy_source").ptr());
    db_->WriteVectors(vnames, vecs, false);
  }

}


void SurfaceBalanceEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  ASSERT(0);
}

} // namespace
} // namespace
