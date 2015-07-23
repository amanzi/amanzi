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

SurfaceBalanceSEB::SurfaceBalanceSEB(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& solution)  :
    PKPhysicalBase(plist,FElist,solution),
    PKDefaultBase(plist,FElist,solution) {

  // set up additional primary variables
  // -- surface energy source
  // Teuchos::ParameterList& esource_sublist =
  //     FElist.sublist("surface_conducted_energy_source");
  // esource_sublist.set("evaluator name", "surface_conducted_energy_source");
  // esource_sublist.set("field evaluator type", "primary variable");

  // -- surface mass source
  Teuchos::ParameterList& wsource_sublist =
      FElist.sublist("surface_mass_source");
  wsource_sublist.set("evaluator name", "surface_mass_source");
  wsource_sublist.set("field evaluator type", "primary variable");

  // -- surface energy temperature
  Teuchos::ParameterList& wtemp_sublist =
      FElist.sublist("surface_mass_source_temperature");
  wtemp_sublist.set("evaluator name", "surface_mass_source_temperature");
  wtemp_sublist.set("field evaluator type", "primary variable");


  // timestep size
  dt_ = plist_->get<double>("max time step", 1.e99);

  // min wind speed
  min_wind_speed_ = plist_->get<double>("minimum wind speed", 1.0);

  // transition snow depth
  snow_ground_trans_ = plist_->get<double>("minimum snow depth", 0.02);
  no_snow_trans_ = plist_->get<double>("zero snow depth", 1.e-5);

  // albedo transition depth
  albedo_trans_ = plist_->get<double>("albedo transition depth", 0.02);

}


void SurfaceBalanceSEB::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)->
      SetComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::RCP<FieldEvaluator> fm;  
  // // requirements: other primary variables
  // S->RequireField("surface_conducted_energy_source", name_)->SetMesh(mesh_)
  //     ->SetComponent("cell", AmanziMesh::CELL, 1);
  // S->RequireFieldEvaluator("surface_conducted_energy_source");
  // fm = S->GetFieldEvaluator("surface_conducted_energy_source");
  // pvfe_esource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  // if (pvfe_esource_ == Teuchos::null) {
  //   Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
  //   Exceptions::amanzi_throw(message);
  // }

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

  S->RequireField("snow_temperature", name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // information from ATS data
  S->RequireFieldEvaluator("surface_temperature");
  S->RequireField("surface_temperature")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("ponded_depth");
  S->RequireField("ponded_depth")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("unfrozen_fraction");
  S->RequireField("unfrozen_fraction")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface_porosity");
  S->RequireField("surface_porosity")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface_vapor_pressure"); // atually mole_fraction not vapor pressure ~ AA
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

  // initialize snow temp
  S->GetFieldData("snow_temperature",name_)->PutScalar(0.);
  S->GetField("snow_temperature", name_)->set_initialized();

  // initialize sources, temps
  //  S->GetFieldData("surface_conducted_energy_source",name_)->PutScalar(0.);
  //  S->GetField("surface_conducted_energy_source",name_)->set_initialized();
  S->GetFieldData("surface_mass_source",name_)->PutScalar(0.);
  S->GetField("surface_mass_source",name_)->set_initialized();
  S->GetFieldData("surface_mass_source_temperature",name_)->PutScalar(273.15);
  S->GetField("surface_mass_source_temperature",name_)->set_initialized();

};


bool SurfaceBalanceSEB::advance(double dt) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_inter_->time()
               << " t1 = " << S_next_->time() << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;
  // Get all data
  // ATS CALCULATED
  S_next_->GetFieldEvaluator("surface_temperature")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_temp =
      *S_next_->GetFieldData("surface_temperature")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& ponded_depth =
      *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("unfrozen_fraction")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& unfrozen_fraction =
      *S_next_->GetFieldData("unfrozen_fraction")->ViewComponent("cell", false);


  S_next_->GetFieldEvaluator("surface_porosity")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_porosity =
      *S_next_->GetFieldData("surface_porosity")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("surface_vapor_pressure")->HasFieldChanged(S_next_.ptr(), name_);
  //  THIS IS MOLE FRACTION OF GAS NEEDS TO BE CONVERTEDT TO VAPOR PRESSURE!
  // Actually mole_fraction not pressure ~AA
  const Epetra_MultiVector& soil_vapor_mol_fraction =
      *S_next_->GetFieldData("surface_vapor_pressure")->ViewComponent("cell", false);

  // MET DATA
  S_next_->GetFieldEvaluator("air_temperature")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& air_temp =
      *S_next_->GetFieldData("air_temperature")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("incoming_shortwave_radiation")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& incoming_shortwave =
      *S_next_->GetFieldData("incoming_shortwave_radiation")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("relative_humidity")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& relative_humidity =
      *S_next_->GetFieldData("relative_humidity")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("wind_speed")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& wind_speed =
      *S_next_->GetFieldData("wind_speed")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("precipitation_rain")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& precip_rain =
      *S_next_->GetFieldData("precipitation_rain")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("precipitation_snow")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& precip_snow =
      *S_next_->GetFieldData("precipitation_snow")->ViewComponent("cell", false);


  // Get output data
  // Epetra_MultiVector& surf_energy_flux =
  //     *S_next_->GetFieldData("surface_conducted_energy_source", name_)->ViewComponent("cell", false);

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

  Epetra_MultiVector& snow_temp =
      *S_next_->GetFieldData("snow_temperature", name_)->ViewComponent("cell", false);

  // Create the SEB data structure
  SurfaceEnergyBalance::LocalData data;
  data.st_energy.dt = dt;
  data.st_energy.AlbedoTrans = albedo_trans_;

  SurfaceEnergyBalance::LocalData data_bare;
  data_bare.st_energy.dt = dt;
  data_bare.st_energy.AlbedoTrans = albedo_trans_;

   data.vp_ground.relative_humidity=1;
   data_bare.vp_ground.relative_humidity=1;


   data.st_energy.Zo=0.005;
   if (air_temp[0][0] > 270){// Little ditty I wrote for the roughness lenght ~ AA 1/10/14
      double Zsmooth = 0.005;
      double Zrough = 0.04;
      double Zfraction = -0.1*air_temp[0][0] + 28;
      if (air_temp[0][0]>=280){
       Zfraction = 0;
       }
     data.st_energy.Zo=(Zsmooth*Zfraction) + (Zrough*(1-Zfraction));
    }

  // loop over all cells and call CalculateSEB_
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    // ATS Calcualted Data
    double density_air = 1.275;       // Density of Air ------------------- [kg/m^3]
    data.st_energy.water_depth = ponded_depth[0][c];
    data.st_energy.water_fraction = unfrozen_fraction[0][c];
    data.st_energy.temp_ground = surf_temp[0][c];
    data.vp_ground.temp = surf_temp[0][c];
    // Convert mol fraction to vapor pressure [moleFraction/atmosphericPressure]
    data.vp_ground.actual_vaporpressure = soil_vapor_mol_fraction[0][c] * data.st_energy.Apa;
    data.st_energy.porrowaLe = surf_porosity[0][c] * density_air * data.st_energy.Le;
    // MET station data
    data.st_energy.temp_air = air_temp[0][c];
    data.st_energy.QswIn = incoming_shortwave[0][c];
    data.st_energy.Us = std::max(wind_speed[0][c], min_wind_speed_);
    data.st_energy.Pr = precip_rain[0][c] * dt; // SEB expects total precip, not rate
    data.st_energy.Ps = precip_snow[0][c] * dt; // SEB expects total precip, not rate
    data.vp_air.temp = air_temp[0][c];
    data.vp_air.relative_humidity = relative_humidity[0][c];
    // STORED INFO FOR SnowEnergyBalanc Model
    data.st_energy.ht_snow = snow_depth[0][c];
    data.st_energy.density_snow = snow_density[0][c];
    data.st_energy.age_snow = days_of_nosnow[0][c];

    // Snow-ground Smoothing
    // -- zero out if just small
    if (data.st_energy.ht_snow < no_snow_trans_)
      data.st_energy.ht_snow = 0.;

    if ((data.st_energy.ht_snow > snow_ground_trans_) ||
        (data.st_energy.ht_snow <= 0)) {
      SurfaceEnergyBalance::SnowEnergyBalance(data);
    } else { // Transition between Snow and bare ground
      double theta=0.0;

      // Fraction of snow covered ground
      //      theta = pow ((data.st_energy.ht_snow / snow_ground_trans_),2);
      theta = data.st_energy.ht_snow / snow_ground_trans_;

      // Calculate as if ht_snow is the min value.
      data.st_energy.ht_snow = snow_ground_trans_;
      SurfaceEnergyBalance::SnowEnergyBalance(data);

      // Calculate as if bare ground
      // ATS Calcualted Data
      data_bare.st_energy.water_depth = ponded_depth[0][c];
      data_bare.st_energy.water_fraction = unfrozen_fraction[0][c];
      data_bare.st_energy.temp_ground = surf_temp[0][c];
      data_bare.vp_ground.temp = surf_temp[0][c];
      //Convert Mol fraction to vapor pressure [molFraction/atmosphericPressure]
      data_bare.vp_ground.actual_vaporpressure = soil_vapor_mol_fraction[0][c] * data_bare.st_energy.Apa; 

      data_bare.st_energy.porrowaLe = surf_porosity[0][c] * density_air * data_bare.st_energy.Le;
      // MET station data
      data_bare.st_energy.temp_air = air_temp[0][c];
      data_bare.st_energy.QswIn = incoming_shortwave[0][c];
      data_bare.st_energy.Us = std::max(wind_speed[0][c], min_wind_speed_);
      data_bare.st_energy.Pr = precip_rain[0][c] * dt; // SEB expects total precip, not rate
      data_bare.st_energy.Ps = precip_snow[0][c] * dt; // SEB expects total precip, not rate
      data_bare.vp_air.temp = air_temp[0][c];
      data_bare.vp_air.relative_humidity = relative_humidity[0][c];
      // STORED INFO FOR SnowEnergyBalanc Model
      data_bare.st_energy.ht_snow = snow_depth[0][c];
      data_bare.st_energy.density_snow = snow_density[0][c];
      data_bare.st_energy.age_snow = days_of_nosnow[0][c];
      data_bare.st_energy.ht_snow = 0.;
      SurfaceEnergyBalance::SnowEnergyBalance(data_bare);

      // Calculating Data for ATS
      data.st_energy.fQc = data.st_energy.fQc * theta + data_bare.st_energy.fQc * (1.-theta);
      data.st_energy.Mr = data.st_energy.Mr * theta + data_bare.st_energy.Mr * (1.-theta);
      data.st_energy.Trw = data.st_energy.Trw * theta + data_bare.st_energy.Trw * (1.-theta);

      // slightly different averaging
      double avg_dens = 0.;
      double avg_age = 0.;
      double ht = 0.;
      if (data_bare.st_energy.ht_snow > 0) {
        avg_dens += data_bare.st_energy.density_snow * data_bare.st_energy.ht_snow;
        avg_age += data_bare.st_energy.age_snow * data_bare.st_energy.ht_snow;
        ht += data_bare.st_energy.ht_snow;
      }
      if (data.st_energy.ht_snow > 0) {
        avg_dens += data.st_energy.density_snow * data.st_energy.ht_snow;
        avg_age += data.st_energy.age_snow * data.st_energy.ht_snow;
        ht += data.st_energy.ht_snow;
      }
      if (ht > 0.) {
        data.st_energy.density_snow = avg_dens / ht;
        data.st_energy.age_snow = avg_age / ht;
      } else {
        data.st_energy.density_snow = data.st_energy.density_frost;
        data.st_energy.age_snow = 0.;
      }

      data.st_energy.ht_snow = data.st_energy.ht_snow * theta + data_bare.st_energy.ht_snow * (1.-theta);
    }

    // STUFF ATS WANTS
    //    surf_energy_flux[0][c] = data.st_energy.fQc;
    surface_water_flux[0][c] = data.st_energy.Mr;
    surf_water_temp[0][c] = data.st_energy.Trw;

    // STUFF SnowEnergyBalance NEEDS STORED FOR NEXT TIME STEP
    snow_depth[0][c] = data.st_energy.ht_snow;
    snow_density[0][c] = data.st_energy.density_snow;
    days_of_nosnow[0][c] = data.st_energy.age_snow;
    snow_temp[0][c] = data.st_energy.temp_snow;
    
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      int rank = mesh_->get_comm()->MyPID();
      Teuchos::RCP<VerboseObject> dcvo = db_->GetVerboseObject(c, rank);
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_HIGH)) {
        *dcvo->os() << "Surface Cell " << c << " SEB:" << std::endl
                    << "  Snow depth, snowtemp = " << data.st_energy.ht_snow << ", " << data.st_energy.temp_snow << std::endl
                    << "  Melt heat = " << data.st_energy.Qm << std::endl
                    << "  ShortWave = " << data.st_energy.fQswIn << std::endl
                    << "  LongWave IN = " << data.st_energy.fQlwIn << std::endl
                    << "  LongWave OUT = " << data.st_energy.fQlwOut << std::endl
                    << "  Latent heat = " << data.st_energy.fQe << std::endl
                    << "  Sensible heat = " << data.st_energy.fQh << std::endl
                    << "GROUND HEAT Qex = " << data.st_energy.fQc << std::endl
                    << "  Ice condensation rate = " << data.st_energy.MIr << std::endl
                    << "WATER SOURCE Mr = " << data.st_energy.Mr << std::endl;
      }
    }
  }

  // Mark primary variables as changed.
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
  //  pvfe_esource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wsource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wtemp_->SetFieldAsChanged(S_next_.ptr());

  // // debug
  // if (vo_->os_OK(Teuchos::VERB_HIGH)) {
  //   std::vector<std::string> vnames;
  //   std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  //   vnames.push_back("air_temp"); vecs.push_back(S_next_->GetFieldData("air_temperature").ptr());
  //   vnames.push_back("Qsw_in"); vecs.push_back(S_next_->GetFieldData("incoming_shortwave_radiation").ptr());
  //   vnames.push_back("precip_rain"); vecs.push_back(S_next_->GetFieldData("precipitation_rain").ptr());
  //   vnames.push_back("precip_snow"); vecs.push_back(S_next_->GetFieldData("precipitation_snow").ptr());
  //   vnames.push_back("soil vapor mol fraction"); vecs.push_back(S_next_->GetFieldData("surface_vapor_pressure").ptr());
  //   vnames.push_back("T_ground"); vecs.push_back(S_next_->GetFieldData("surface_temperature").ptr());
  //   vnames.push_back("water_source"); vecs.push_back(S_next_->GetFieldData("surface_mass_source").ptr());
  //   //    vnames.push_back("e_source"); vecs.push_back(S_next_->GetFieldData("surface_conducted_energy_source").ptr());
  //   db_->WriteVectors(vnames, vecs, true);
  // }

  return false;
};


} // namespace
} // namespace
