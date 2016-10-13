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

#include "surface_balance_SEB_VPL.hh"
#include "SnowEnergyBalance_VPL.hh"

namespace Amanzi {
namespace SurfaceBalance {

SurfaceBalanceSEBVPL::SurfaceBalanceSEBVPL(Teuchos::ParameterList& pk_tree,
                                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                     const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& solution):
  PK(pk_tree, global_list,  S, solution),
  PK_Physical_Default(pk_tree, global_list,  S, solution)
 {

   Teuchos::ParameterList& FElist = global_list->sublist("field evaluators");

  // -- surface mass source
  Teuchos::ParameterList& wsource_sublist =
      FElist.sublist("surface_mass_source");
  wsource_sublist.set("evaluator name", "surface_mass_source");
  wsource_sublist.set("field evaluator type", "primary variable");

  // -- VaporMassSource for VaporFlux at cell center 
  Teuchos::ParameterList& w_v_source_sublist =
    FElist.sublist("mass_source");
  w_v_source_sublist.set("evaluator name", "mass_source");
  w_v_source_sublist.set("field evaluator type", "primary variable");

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


void SurfaceBalanceSEBVPL::Setup(const Teuchos::Ptr<State>& S) {
  PK_Physical_Default::Setup(S);
  subsurf_mesh_ = S->GetMesh(); // needed for VPL, which is treated as subsurface source

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
  //   Errors::Message message("SurfaceBalanceSEBVPL: error, failure to initialize primary variable");
  //   Exceptions::amanzi_throw(message);
  // }

  S->RequireField("surface_mass_source", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_mass_source");
  fm = S->GetFieldEvaluator("surface_mass_source");
  pvfe_wsource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_wsource_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEBVPL: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField("mass_source", name_)->SetMesh(subsurf_mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("mass_source");
  fm = S->GetFieldEvaluator("mass_source");
  pvfe_w_v_source_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_w_v_source_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField("surface_mass_source_temperature", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_mass_source_temperature");
  fm = S->GetFieldEvaluator("surface_mass_source_temperature");
  pvfe_wtemp_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_wtemp_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEBVPL: error, failure to initialize primary variable");
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

  S->RequireField("stored_surface_pressure", name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("stored_SWE", name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // information from ATS data
  S->RequireFieldEvaluator("surface_temperature");
  S->RequireField("surface_temperature")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("ponded_depth");
  S->RequireField("ponded_depth")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

//   S->RequireFieldEvaluator("saturation_liquid");
   S->RequireField("saturation_liquid")->SetMesh(subsurf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  
    S->RequireFieldEvaluator("surface_pressure");
  S->RequireField("surface_pressure")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);  

  S->RequireFieldEvaluator("unfrozen_fraction");
  S->RequireField("unfrozen_fraction")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface_porosity");
  S->RequireField("surface_porosity")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface_vapor_pressure"); // actually mole_fraction not vapor pressure ~ AA
  S->RequireField("surface_vapor_pressure")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
};

// initialize ICs
void SurfaceBalanceSEBVPL::Initialize(const Teuchos::Ptr<State>& S) {
  // this call specifies snow depth
  PK_Physical_Default::Initialize(S);

  // initialize snow density
  S->GetFieldData("snow_density",name_)->PutScalar(100.);
  S->GetField("snow_density", name_)->set_initialized();

  // initialize days of no snow
  S->GetFieldData("days_of_nosnow",name_)->PutScalar(0.);
  S->GetField("days_of_nosnow", name_)->set_initialized();

  // initialize snow temp
  S->GetFieldData("snow_temperature",name_)->PutScalar(0.);
  S->GetField("snow_temperature", name_)->set_initialized();

  // initialize stored surface pressure
  S->GetFieldData("stored_surface_pressure",name_)->PutScalar(0.);
  S->GetField("stored_surface_pressure", name_)->set_initialized();
  S->GetFieldData("stored_SWE",name_)->PutScalar(0.);
  S->GetField("stored_SWE", name_)->set_initialized();
  // initialize sources, temps
  //  S->GetFieldData("surface_conducted_energy_source",name_)->PutScalar(0.);
  //  S->GetField("surface_conducted_energy_source",name_)->set_initialized();
  S->GetFieldData("surface_mass_source",name_)->PutScalar(0.);
  S->GetField("surface_mass_source",name_)->set_initialized();
    S->GetFieldData("mass_source",name_)->PutScalar(0.);
  S->GetField("mass_source",name_)->set_initialized();
  S->GetFieldData("surface_mass_source_temperature",name_)->PutScalar(273.15);
  S->GetField("surface_mass_source_temperature",name_)->set_initialized();

};


  bool SurfaceBalanceSEBVPL::AdvanceStep(double t_old, double t_new, bool reinit) {

    double dt = t_new - t_old;

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

//  S_next_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& saturation_liquid =
    *S_next_->GetFieldData("saturation_liquid")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("surface_pressure")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surface_pressure =
      *S_next_->GetFieldData("surface_pressure")->ViewComponent("cell", false);

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
  
  Epetra_MultiVector& surface_vapor_flux =
      *S_next_->GetFieldData("mass_source", name_)->ViewComponent("cell", false);
  surface_vapor_flux.PutScalar(0.);

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

  Epetra_MultiVector& stored_surface_pressure =
      *S_next_->GetFieldData("stored_surface_pressure", name_)->ViewComponent("cell", false);
  Epetra_MultiVector& stored_SWE =
      *S_next_->GetFieldData("stored_SWE", name_)->ViewComponent("cell", false);

  //0 Create the SEB data structure
  SurfaceEnergyBalance_VPL::LocalData data;
  data.st_energy.dt = dt;
  data.st_energy.AlbedoTrans = albedo_trans_;

  SurfaceEnergyBalance_VPL::LocalData data_bare;
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
//   data.st_energy.stored_fQe = 9999999;

  // loop over all cells and call CalculateSEB_
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    // ATS Calcualted Data
    double density_air = 1.275;       // Density of Air ------------------- [kg/m^3]
    data.st_energy.water_depth = ponded_depth[0][c]; 
    
    AmanziMesh::Entity_ID subsurf_f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
    AmanziMesh::Entity_ID_List cells;
    subsurf_mesh_->face_get_cells(subsurf_f, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);
    data.st_energy.saturation_liquid = saturation_liquid[0][cells[0]];

    data.st_energy.surface_pressure = surface_pressure[0][c];
    //data.st_energy.stored_pressure = surface_pressure[0][c];
    data.st_energy.water_fraction = unfrozen_fraction[0][c];
    data.st_energy.temp_ground = surf_temp[0][c];
    data.vp_ground.temp = surf_temp[0][c];
    // Convert mol fraction to vapor pressure [moleFraction/atmosphericPressure]
    data.vp_ground.actual_vaporpressure = soil_vapor_mol_fraction[0][c] * data.st_energy.Apa;
    data.st_energy.surface_porosity = surf_porosity[0][c];
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
    data.st_energy.stored_surface_pressure = stored_surface_pressure[0][c];

    // Snow-ground Smoothing
    // -- zero out if just small
    if (data.st_energy.ht_snow < no_snow_trans_)
      data.st_energy.ht_snow = 0.;

    if ((data.st_energy.ht_snow > snow_ground_trans_) ||
        (data.st_energy.ht_snow <= 0)) {
      SurfaceEnergyBalance_VPL::SnowEnergyBalance(data);
    } else { // Transition between Snow and bare ground
      double theta=0.0;

      // Fraction of snow covered ground
      //      theta = pow ((data.st_energy.ht_snow / snow_ground_trans_),2);
      theta = data.st_energy.ht_snow / snow_ground_trans_;

      // Calculate as if ht_snow is the min value.
      data.st_energy.ht_snow = snow_ground_trans_;
      SurfaceEnergyBalance_VPL::SnowEnergyBalance(data);

      // Calculate as if bare ground
      // ATS Calcualted Data
      data_bare.st_energy.water_depth = ponded_depth[0][c];

    AmanziMesh::Entity_ID subsurf_f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
    AmanziMesh::Entity_ID_List cells;
    subsurf_mesh_->face_get_cells(subsurf_f, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);
      data_bare.st_energy.saturation_liquid = saturation_liquid[0][cells[0]];  

      data_bare.st_energy.surface_pressure = surface_pressure[0][c];
      data_bare.st_energy.stored_surface_pressure = stored_surface_pressure[0][c];
      data_bare.st_energy.water_fraction = unfrozen_fraction[0][c];
      data_bare.st_energy.temp_ground = surf_temp[0][c];
      data_bare.vp_ground.temp = surf_temp[0][c];
      //Convert Mol fraction to vapor pressure [molFraction/atmosphericPressure]
      data_bare.vp_ground.actual_vaporpressure = soil_vapor_mol_fraction[0][c] * data_bare.st_energy.Apa; 

      data_bare.st_energy.surface_porosity = surf_porosity[0][c];
      // MET station data
      data_bare.st_energy.temp_air = air_temp[0][c];
      data_bare.st_energy.QswIn = incoming_shortwave[0][c];
      data_bare.st_energy.Us = std::max(wind_speed[0][c], min_wind_speed_);
      data_bare.st_energy.Pr = precip_rain[0][c] * dt; // SEB expects total precip, not rate
      data_bare.st_energy.Ps = precip_snow[0][c] * dt; // SEB expects total precip, not rate
      data_bare.vp_air.temp = air_temp[0][c];
      data_bare.vp_air.relative_humidity = relative_humidity[0][c];
      // STORED INFO FOR SnowEnergyBalanc Model
      //      data_bare.st_energy.ht_snow = snow_depth[0][c];
      data_bare.st_energy.density_snow = snow_density[0][c];
      data_bare.st_energy.age_snow = days_of_nosnow[0][c];
      data_bare.st_energy.ht_snow = 0.;
      data.st_energy.stored_surface_pressure = stored_surface_pressure[0][c];

      SurfaceEnergyBalance_VPL::SnowEnergyBalance(data_bare);

      // Calculating Data for ATS
      data.st_energy.fQc = data.st_energy.fQc * theta + data_bare.st_energy.fQc * (1.-theta);
      data.st_energy.Mr = data.st_energy.Mr * theta + data_bare.st_energy.Mr * (1.-theta);
      data.st_energy.Trw = data.st_energy.Trw * theta + data_bare.st_energy.Trw * (1.-theta);

      // slightly different averaging
      double avg_dens(0.), avg_age(0.), ht(0.);
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
   //     surf_energy_flux[0][c] = data.st_energy.fQc;  
    surface_water_flux[0][c] = data.st_energy.Mr;
    surf_water_temp[0][c] = data.st_energy.Trw;

    // surface vapor flux is treated as a volumetric source for the subsurface.
  //  AmanziMesh::Entity_ID subsurf_f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
  //  AmanziMesh::Entity_ID_List cells;
  //  subsurf_mesh_->face_get_cells(subsurf_f, AmanziMesh::OWNED, &cells);
  //  ASSERT(cells.size() == 1);
    // surface mass sources are in m^3 water / (m^2 s)
    // subsurface mass sources are in mol water / (m^3 s)
    surface_vapor_flux[0][cells[0]] = data.st_energy.SurfaceVaporFlux 
      * mesh_->cell_volume(c) * data.st_energy.density_w / 0.0180153
      / subsurf_mesh_->cell_volume(cells[0]);

    // STUFF SnowEnergyBalance NEEDS STORED FOR NEXT TIME STEP
    snow_depth[0][c] = data.st_energy.ht_snow;
    snow_density[0][c] = data.st_energy.density_snow;
    days_of_nosnow[0][c] = data.st_energy.age_snow;
    snow_temp[0][c] = data.st_energy.temp_snow;
    // This is intended to keep the capillary pressure for internal interations (evaluatur) constant ~AA
   // stored_surface_pressure[0][c] = surface_pressure[0][c];  
    stored_surface_pressure[0][c] = surface_pressure[0][c]; 
    std::cout<<"SEB:: Surface_pressure: "<<surface_pressure[0][c]<<"  Stored Surface Pressure: "<<stored_surface_pressure[0][c]<<std::endl;
    stored_SWE[0][c] = data.st_energy.SWE;    
 
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Snow depth, snowtemp = " << data.st_energy.ht_snow << ", " << data.st_energy.temp_snow << std::endl;
      *vo_->os() << "Melt heat, Melt water temp = " << data.st_energy.Qm <<", " << data.st_energy.Trw << std::endl;
      *vo_->os() << "ShortWave = " << data.st_energy.fQswIn << std::endl;
      *vo_->os() << "LongWave IN = " << data.st_energy.fQlwIn << std::endl;
      *vo_->os() << "LongWave OUT = " << data.st_energy.fQlwOut << std::endl;
      *vo_->os() << "Latent heat = " << data.st_energy.fQe << std::endl;
      *vo_->os() << "Sensible heat = " << data.st_energy.fQh << std::endl;
      *vo_->os() << "GROUND HEAT Qex = " << data.st_energy.fQc << std::endl;
      *vo_->os() << "Ice condensation rate = " << data.st_energy.MIr << std::endl;
      *vo_->os() << "ALBEDO = " << data.st_energy.albedo_value << std::endl;
    }

  }

  // Mark primary variables as changed.
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
  //  pvfe_esource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wsource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_w_v_source_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wtemp_->SetFieldAsChanged(S_next_.ptr());

  // debug
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("air_temp"); vecs.push_back(S_next_->GetFieldData("air_temperature").ptr());
    vnames.push_back("Qsw_in"); vecs.push_back(S_next_->GetFieldData("incoming_shortwave_radiation").ptr());
    vnames.push_back("precip_rain"); vecs.push_back(S_next_->GetFieldData("precipitation_rain").ptr());
    vnames.push_back("precip_snow"); vecs.push_back(S_next_->GetFieldData("precipitation_snow").ptr());
    vnames.push_back("soil vapor mol fraction"); vecs.push_back(S_next_->GetFieldData("surface_vapor_pressure").ptr());
    vnames.push_back("T_ground"); vecs.push_back(S_next_->GetFieldData("surface_temperature").ptr());
    vnames.push_back("water_source"); vecs.push_back(S_next_->GetFieldData("surface_mass_source").ptr());
    vnames.push_back("surface_vapor_source"); vecs.push_back(S_next_->GetFieldData("mass_source").ptr());
    //    vnames.push_back("e_source"); vecs.push_back(S_next_->GetFieldData("surface_conducted_energy_source").ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  return false;
};


} // namespace
} // namespace
