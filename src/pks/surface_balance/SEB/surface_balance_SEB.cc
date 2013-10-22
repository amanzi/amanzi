/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 * Author: Ethan Coon, Adam Atchley, Satish Karra
 *
 * DOCUMENT ME:
 * Surface Energy Balance for Snow Surface and Ground Surface
 * Calculates Energy flux, rate or water, and water temperature
 * entering through the surface skin.  Snow surface energy balance
 * is calculated at equilibrium with ground/surface water and Air.
 *
 * 0=(1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qc(Ts) + Qe(Ts)
 * Qc = the energy delived to the subsurface
 * The rate of water entering the surface skin occurs only when Ts > 0
 * In which case Ts is set = 0 and the excess energy = Qm and the melt rate (Mr) is
 * delivered to the surface skin via:
 * Mr = Qm/(ROWw*Hf)
 * ROWw = density of water
 * Hf = latent heat of fusion
 * The temperature of water in assumed to be 0 C
 *
 * In cases without snow the energy balance equations is:
 *
 * Qex = 0=(1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qe(Ts)
 * Qex is the energy derived to the subsurface skin
 * All water entering the surface skin is assumed to be precipitated
 * or condensed on the surface and has a temperature of Air.
 *
 *
 * ------------------------------------------------------------------------- */


#include "surface_top_cells_evaluator.hh"

#include "surface_balance_SEB.hh"
#include "SnowEnergyBalance.hh"

namespace Amanzi {
namespace SurfaceBalance {

RegisteredPKFactory<SurfaceBalanceSEB> SurfaceBalanceSEB::reg_("surface balance SEB");


SurfaceBalanceSEB::SurfaceBalanceSEB(Teuchos::ParameterList& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& solution)  :
    PKPhysicalBase(plist,FElist,solution),
    PKDefaultBase(plist,FElist,solution) {

  // set up additional primary variables
  // -- surface energy source
  Teuchos::ParameterList& esource_sublist = FElist.sublist("surface_conducted_energy_source");
  esource_sublist.set("evaluator name", "surface_conducted_energy_source");
  esource_sublist.set("field evaluator type", "primary variable");

  // -- surface mass source
  Teuchos::ParameterList& wsource_sublist = FElist.sublist("surface_mass_source");
  wsource_sublist.set("evaluator name", "surface_mass_source");
  wsource_sublist.set("field evaluator type", "primary variable");

  // -- surface energy temperature
  Teuchos::ParameterList& wtemp_sublist = FElist.sublist("surface_mass_source_temperature");
  wtemp_sublist.set("evaluator name", "surface_mass_source_temperature");
  wtemp_sublist.set("field evaluator type", "primary variable");


  // timestep size
  dt_ = plist.get<double>("max time step", 1.e99);

  // min wind speed
  min_wind_speed_ = plist.get<double>("minimum wind speed", 0.5);
}


void SurfaceBalanceSEB::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)->
      SetComponent("cell", AmanziMesh::CELL, 1);

  // requirements: other primary variables
  S->RequireField("surface_conducted_energy_source", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_conducted_energy_source");
  Teuchos::RCP<FieldEvaluator> fm = S->GetFieldEvaluator("surface_conducted_energy_source");
  pvfe_esource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_esource_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField("surface_mass_source", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_mass_source");
  fm = S->GetFieldEvaluator("surface_mass_source");
  pvfe_wsource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_wsource_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField("surface_mass_source_temperature", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_mass_source_temperature");
  fm = S->GetFieldEvaluator("surface_mass_source_temperature");
  pvfe_wtemp_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_wtemp_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  // requirements: independent variables (data from MET)
  S->RequireFieldEvaluator("incoming_shortwave_radiation");
  S->RequireField("incoming_shortwave_radiation")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("air_temperature");
  S->RequireField("air_temperature")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("relative_humidity");
  S->RequireField("relative_humidity")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("wind_speed");
  S->RequireField("wind_speed")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("precipitation_rain");
  S->RequireField("precipitation_rain")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("precipitation_snow");
  S->RequireField("precipitation_snow")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirements: stored secondary variables
  S->RequireField("snow_density", name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField("days_of_nosnow", name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // information from ATS data
  S->RequireFieldEvaluator("surface_temperature");
  S->RequireField("surface_temperature")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("ponded_depth");
  S->RequireField("ponded_depth")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface_porosity");
  S->RequireField("surface_porosity")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface_vapor_pressure");
  S->RequireField("surface_vapor_pressure")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
};

// initialize ICs
void SurfaceBalanceSEB::initialize(const Teuchos::Ptr<State>& S) {
  // this call specifies snow depth
  PKPhysicalBase::initialize(S);

  // initialize snow density
  S->GetFieldData("snow_density",name_)->PutScalar(100.);
  S->GetField("snow_density", name_)->set_initialized();

  // initialize days of no snow
  S->GetFieldData("days_of_nosnow",name_)->PutScalar(0.);
  S->GetField("days_of_nosnow", name_)->set_initialized();

};


bool SurfaceBalanceSEB::advance(double dt) {
  // Create the SEB data structure
  SurfaceEnergyBalance::LocalData data;

  // Populate the SEB data structure with constants
  data.st_energy.stephB = 0.00000005670373;// Stephan-boltzmann Constant ------- [W/m^2 K^4]
  data.st_energy.Hf = 333500.0;            // Heat of fusion for melting snow -- [J/kg]
  data.st_energy.Ls = 2834000.0;           // Latent heat of sublimation ------- [J/kg]
  data.st_energy.Le = 2497848.;            // Latent heat of vaporization ------ [J/kg]
  data.st_energy.SEs = 0.98;               // Surface Emissivity for snow  ----- [-] ** From P. ReVelle (Thesis)
  data.st_energy.Zr = 2.0;                 // Referance ht of wind speed ------- [m]
  data.st_energy.Zo = 0.005;               // Roughness length  ---------------- [m]
  data.st_energy.VKc = 0.41;               // Von Karman Constant -------------- [-]
  double Cp = 1004.0;               // Specific heat of air ------------- [J/K kg]
  data.st_energy.Apa = 100;                // Atmospheric Pressure ------------- [KPa]

  data.st_energy.density_w = 1000;         // Density of Water ----------------- [kg/m^3]
  double density_air = 1.275;       // Density of Air ------------------- [kg/m^3]
  data.st_energy.density_frost = 800;      // Density of Frost (condensation) -- [kg/m^3]
  data.st_energy.density_freshsnow = 100;  // Density of Freshly fallebn snow -- [kg/m^3]

  data.st_energy.gZr = 9.807*data.st_energy.Zr;
  data.st_energy.rowaCp = density_air*Cp;
  data.st_energy.rowaLs = density_air*data.st_energy.Ls;
  data.st_energy.rowaLe = density_air*data.st_energy.Le;

  data.vp_snow.relative_humidity = 1.;
  data.vp_ground.relative_humidity = 1.;

  data.st_energy.Dt = dt;

  // Get all data
  // ATS CALCULATED
  S_inter_->GetFieldEvaluator("surface_temperature")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& surf_temp =
      *S_inter_->GetFieldData("surface_temperature")->ViewComponent("cell", false);

  S_inter_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& ponded_depth =
      *S_inter_->GetFieldData("ponded_depth")->ViewComponent("cell", false);

  S_inter_->GetFieldEvaluator("surface_porosity")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& surf_porosity =
      *S_inter_->GetFieldData("surface_porosity")->ViewComponent("cell", false);

  S_inter_->GetFieldEvaluator("surface_vapor_pressure")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& soil_vapor_pressure =
      *S_inter_->GetFieldData("surface_vapor_pressure")->ViewComponent("cell", false);

  // MET DATA
  S_inter_->GetFieldEvaluator("air_temperature")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& air_temp =
      *S_inter_->GetFieldData("air_temperature")->ViewComponent("cell", false);

  S_inter_->GetFieldEvaluator("incoming_shortwave_radiation")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& incomming_shortwave =
      *S_inter_->GetFieldData("incoming_shortwave_radiation")->ViewComponent("cell", false);

  S_inter_->GetFieldEvaluator("relative_humidity")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& relative_humidity =
      *S_inter_->GetFieldData("relative_humidity")->ViewComponent("cell", false);

  S_inter_->GetFieldEvaluator("wind_speed")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& wind_speed =
      *S_inter_->GetFieldData("wind_speed")->ViewComponent("cell", false);

  S_inter_->GetFieldEvaluator("precipitation_rain")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& precip_rain =
      *S_inter_->GetFieldData("precipitation_rain")->ViewComponent("cell", false);

  S_inter_->GetFieldEvaluator("precipitation_snow")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& precip_snow =
      *S_inter_->GetFieldData("precipitation_snow")->ViewComponent("cell", false);


  // Get output data
  Epetra_MultiVector& surf_energy_flux =
      *S_next_->GetFieldData("surface_conducted_energy_source", name_)->ViewComponent("cell", false);

  Epetra_MultiVector& surface_water_flux =
      *S_next_->GetFieldData("surface_mass_source", name_)->ViewComponent("cell", false);

  Epetra_MultiVector& surf_water_temp =
      *S_next_->GetFieldData("surface_mass_source_temperature", name_)->ViewComponent("cell", false);

  Epetra_MultiVector& snow_depth =
      *S_next_->GetFieldData("snow_depth", name_)->ViewComponent("cell", false);

  Epetra_MultiVector& snow_density =
      *S_next_->GetFieldData("snow_density", name_)->ViewComponent("cell", false);

  Epetra_MultiVector& days_of_nosnow =
      *S_next_->GetFieldData("days_of_nosnow", name_)->ViewComponent("cell", false);


  // loop over all cells and call CalculateSEB_
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    // ATS Calcualted Data
    data.st_energy.Tb = surf_temp[0][c];
    data.st_energy.water_depth = ponded_depth[0][c];
    data.st_energy.por = surf_porosity[0][c];
    data.vp_ground.temp = surf_temp[0][c];
    data.vp_ground.actual_vaporpressure = soil_vapor_pressure[0][c]; // FIX THIS   FIX THIS FIX THIS    FIX THIS 11111111
    data.st_energy.porrowaLe = data.st_energy.por*density_air*data.st_energy.Le;
    // MET station data
    data.st_energy.air_temp = air_temp[0][c];
    data.st_energy.QswIn = incomming_shortwave[0][c];
    data.st_energy.Us = std::max(wind_speed[0][c], min_wind_speed_);
    data.st_energy.Pr = precip_rain[0][c] * dt; // SEB expects total precip, not rate
    data.st_energy.Ps = precip_snow[0][c] * dt; // SEB expects total precip, not rate
    data.vp_air.temp = air_temp[0][c];
    data.vp_air.relative_humidity = relative_humidity[0][c];

    // STORED INFO FOR SnowEnergyBalanc Model
    data.st_energy.ht_snow = snow_depth[0][c];
    data.st_energy.density_snow = snow_density[0][c];
    data.st_energy.nosnowdays = days_of_nosnow[0][c];

    SurfaceEnergyBalance::SnowEnergyBalance(data);//.........);

    // STUFF ATS WANTS
    // *** Ask Ethan about "surf_water_flux" naming convention ***
    surf_energy_flux[0][c] = data.st_energy.fQc;
    surface_water_flux[0][c]  = data.st_energy.Mr;
    surf_water_temp[0][c]  = data.st_energy.Trw;

    // STUFF SnowEnergyBalance NEEDS STORED FOR NEXT TIME STEP
    snow_depth[0][c]=data.st_energy.ht_snow;
    snow_density[0][c]=data.st_energy.density_snow;
    days_of_nosnow[0][c]=data.st_energy.nosnowdays;


  }

  // Mark primary variables as changed.
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
  pvfe_esource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wsource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wtemp_->SetFieldAsChanged(S_next_.ptr());

};


} // namespace
} // namespace
