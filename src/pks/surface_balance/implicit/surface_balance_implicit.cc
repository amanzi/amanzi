/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Adam Atchley, Satish Karra

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   ------------------------------------------------------------------------- */

#include "boost/algorithm/string/predicate.hpp"

#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"
#include "surface_balance_implicit.hh"

namespace Amanzi {
namespace SurfaceBalance {

SurfaceBalanceImplicit::SurfaceBalanceImplicit(Teuchos::ParameterList& pk_tree,
                                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                     const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& solution):
  PK(pk_tree, global_list,  S, solution),
  PK_PhysicalBDF_Default(pk_tree, global_list,  S, solution),
  modify_predictor_advance_(false)
{
  if (!plist_->isParameter("conserved quantity suffix"))
    plist_->set("conserved quantity suffix", "snow_depth");

  Teuchos::ParameterList& FElist = S->FEList();

  // set up additional primary variables -- this is very hacky...
  // -- surface energy source
  if (domain_ == "surface") {
    domain_ss_ = plist_->get<std::string>("subsurface domain name", "domain");
  } else if (boost::starts_with(domain_, "surface_")) {
    domain_ss_ = plist_->get<std::string>("subsurface domain name", domain_.substr(8,domain_.size()));
  } else if (boost::ends_with(domain_, "_surface")) {
    domain_ss_ = plist_->get<std::string>("subsurface domain name", domain_.substr(0,domain_.size()-8));
  } else {
    plist_->get<std::string>("subsurface domain name");
  }

  // -- surface mass,energy sources
  Teuchos::ParameterList& wsource_sublist =
    FElist.sublist(Keys::getKey(domain_,"mass_source"));
  wsource_sublist.set("evaluator name", Keys::getKey(domain_,"mass_source"));
  wsource_sublist.set("field evaluator type", "primary variable");
  Teuchos::ParameterList& esource_sublist =
    FElist.sublist(Keys::getKey(domain_,"total_energy_source"));
  esource_sublist.set("evaluator name", Keys::getKey(domain_,"total_energy_source"));
  esource_sublist.set("field evaluator type", "primary variable");

  // -- subsurface mass source for VaporFlux at cell center
  Teuchos::ParameterList& w_v_source_sublist =
    FElist.sublist(Keys::getKey(domain_ss_,"mass_source"));
  w_v_source_sublist.set("evaluator name", Keys::getKey(domain_ss_,"mass_source"));
  w_v_source_sublist.set("field evaluator type", "primary variable");
  Teuchos::ParameterList& e_sub_source_sublist =
    FElist.sublist(Keys::getKey(domain_ss_,"total_energy_source"));
  e_sub_source_sublist.set("evaluator name", Keys::getKey(domain_ss_,"total_energy_source"));
  e_sub_source_sublist.set("field evaluator type", "primary variable");

  // Derivatives for PC
  eval_derivatives_ = plist_->get<bool>("evaluate source derivatives", true);

  // min wind speed
  min_wind_speed_ = plist_->get<double>("minimum wind speed [m/s]?", 1.0);
  wind_speed_ref_ht_ = plist_->get<double>("wind speed reference height [m]", 2.0);

  // roughness parameters
  roughness_bare_ground_ = plist_->get<double>("roughness length of bare ground [m]");
  roughness_snow_covered_ground_ = plist_->get<double>("roughness length of snow-covered ground [m]");
  
  // implicit/explicit snow precip
  implicit_snow_ = plist_->get<bool>("implicit snow precipitation", false);

  // Reading in Longwave Radation
  longwave_input_ = FElist.isParameter(Keys::getKey(domain_,"incoming_longwave_radiation"));
  
  // transition snow depth
  snow_ground_trans_ = plist_->get<double>("snow-ground transitional depth", 0.02);
  min_snow_trans_ = plist_->get<double>("minimum snow transitional depth", 1.e-8);
  if (min_snow_trans_ < 0. || snow_ground_trans_ < min_snow_trans_) {
    Errors::Message message("Invalid parameters: snow-ground transitional depth or minimum snow transitional depth.");
    Exceptions::amanzi_throw(message);
  }

  // shortwave radiation key changes if shaded
  sw_incoming_key_ = Keys::readKey(*plist_, domain_,"incoming shortwave radiation", "incoming_shortwave_radiation");

  // modify predictor by calling advance -- this is cheap and sets up BCs
  // correctly for subsurface's call to ModifyPredictorConsistentFaces()
  modify_predictor_advance_ = plist_->get<bool>("modify predictor by advancing", false);
}

// main methods
// -- Setup data.
void
SurfaceBalanceImplicit::Setup(const Teuchos::Ptr<State>& S) {
  PK_PhysicalBDF_Default::Setup(S);
  subsurf_mesh_ = S->GetMesh(domain_ss_); // needed for VPL, which is treated as subsurface source

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)->
      SetComponent("cell", AmanziMesh::CELL, 1);

  // requirements: other primary variables
  Teuchos::RCP<FieldEvaluator> fm;
  S->RequireField(Keys::getKey(domain_,"mass_source"), name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_,"mass_source"));
  fm = S->GetFieldEvaluator(Keys::getKey(domain_,"mass_source"));
  pvfe_wsource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_wsource_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField(Keys::getKey(domain_,"total_energy_source"), name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_,"total_energy_source"));
  fm = S->GetFieldEvaluator(Keys::getKey(domain_,"total_energy_source"));
  pvfe_esource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_esource_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField(Keys::getKey(domain_ss_,"mass_source"), name_)->SetMesh(subsurf_mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_ss_,"mass_source"));
  fm = S->GetFieldEvaluator(Keys::getKey(domain_ss_,"mass_source"));
  pvfe_w_sub_source_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_w_sub_source_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }
  
  S->RequireField(Keys::getKey(domain_ss_,"total_energy_source"), name_)->SetMesh(subsurf_mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_ss_,"total_energy_source"));
  fm = S->GetFieldEvaluator(Keys::getKey(domain_ss_,"total_energy_source"));
  pvfe_e_sub_source_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_e_sub_source_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  // requirements: source derivatives
  if (eval_derivatives_) {
    Key key_cond_temp = Keys::getDerivKey(Keys::getKey(domain_,"total_energy_source"),
            Keys::getKey(domain_,"temperature"));
    S->RequireField(key_cond_temp, name_)->SetMesh(mesh_)
       ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // requirements: diagnostic variables
  S->RequireField(Keys::getKey(domain_,"albedo"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"albedo"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"evaporative_flux"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"evaporative_flux"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"qE_latent_heat"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"qE_latent_heat"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"qE_sensible_heat"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"qE_sensible_heat"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"qE_lw_out"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"qE_lw_out"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"qE_melt_out"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"qE_melt_out"),name_)->set_io_checkpoint(false);
  
  // requirements: independent variables (data from MET)
  S->RequireFieldEvaluator(sw_incoming_key_);
  S->RequireField(sw_incoming_key_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  if (longwave_input_) {
    S->RequireFieldEvaluator(Keys::getKey(domain_,"incoming_longwave_radiation"));
    S->RequireField(Keys::getKey(domain_,"incoming_longwave_radiation"))->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
  } else {
    // it is a diagnostic
    S->RequireField(Keys::getKey(domain_,"incoming_longwave_radiation"), name_)->SetMesh(mesh_)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  S->RequireFieldEvaluator(Keys::getKey(domain_,"air_temperature"));
  S->RequireField(Keys::getKey(domain_,"air_temperature"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"relative_humidity"));
  S->RequireField(Keys::getKey(domain_,"relative_humidity"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"wind_speed"));
  S->RequireField(Keys::getKey(domain_,"wind_speed"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"precipitation_rain"));
  S->RequireField(Keys::getKey(domain_,"precipitation_rain"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"precipitation_snow"));
  S->RequireField(Keys::getKey(domain_,"precipitation_snow"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirements: stored secondary variables
  S->RequireField(Keys::getKey(domain_,"snow_density"), name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField(Keys::getKey(domain_,"snow_age"), name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField(Keys::getKey(domain_,"snow_temperature"), name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField(Keys::getKey(domain_,"stored_SWE"), name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"temperature"));
  S->RequireField(Keys::getKey(domain_,"temperature"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"pressure"));
  S->RequireField(Keys::getKey(domain_,"pressure"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"ponded_depth"));
  S->RequireField(Keys::getKey(domain_,"ponded_depth"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // cannot work with current implementation because WRMs are hard coded into
  //  flow PKs.  Design fail.  Future versions will need this call.
  //  S->RequireFieldEvaluator(Keys::getKey(domain_ss_, "saturation_gas"));
  S->RequireField(Keys::getKey(domain_ss_,"saturation_gas"))->SetMesh(subsurf_mesh_)
       ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"unfrozen_fraction"));
  S->RequireField(Keys::getKey(domain_,"unfrozen_fraction"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_ss_,"porosity"));
  S->RequireField(Keys::getKey(domain_ss_,"porosity"))->SetMesh(subsurf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

}

// -- Initialize owned (dependent) variables.
void
SurfaceBalanceImplicit::Initialize(const Teuchos::Ptr<State>& S) {
  PK_PhysicalBDF_Default::Initialize(S);

  // initialize snow density, age
  ASSERT(plist_->isSublist("initial condition"));
  Teuchos::ParameterList& ic_list = plist_->sublist("initial condition");

  SEBPhysics::ModelParams seb_params;
  
  if (!S->GetField(Keys::getKey(domain_,"snow_density"))->initialized()) {
    if (ic_list.isParameter("restart file")) {
      // initialize density, age from restart file
      S->GetField(Keys::getKey(domain_,"snow_density"), name_)->Initialize(ic_list);
      S->GetField(Keys::getKey(domain_,"snow_density"), name_)->set_initialized();
      S->GetField(Keys::getKey(domain_,"snow_age"), name_)->Initialize(ic_list);
      S->GetField(Keys::getKey(domain_,"snow_age"), name_)->set_initialized();
    } else {
      // initialize density to fresh powder, age to 0
      S->GetFieldData(Keys::getKey(domain_,"snow_density"),name_)->PutScalar(seb_params.density_freshsnow);
      S->GetField(Keys::getKey(domain_,"snow_density"), name_)->set_initialized();
      S->GetFieldData(Keys::getKey(domain_,"snow_age"),name_)->PutScalar(0.);
      S->GetField(Keys::getKey(domain_,"snow_age"), name_)->set_initialized();

    }
  }

  // initialize swe consistently with snow height and density

  Epetra_MultiVector& swe = *S->GetFieldData(Keys::getKey(domain_,"stored_SWE"),name_)->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_ht = *S->GetFieldData(key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_dens = *S->GetFieldData(Keys::getKey(domain_,"snow_density"))->ViewComponent("cell",false);
  for (int c=0; c!=swe.MyLength(); ++c) {
    swe[0][c] = snow_ht[0][c] * snow_dens[0][c] / seb_params.density_water;
  }
  S->GetField(Keys::getKey(domain_,"stored_SWE"), name_)->set_initialized();

  // initialize snow temp
  S->GetFieldData(Keys::getKey(domain_,"snow_temperature"),name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_,"snow_temperature"), name_)->set_initialized();

  // initialize sources, temps
  S->GetFieldData(Keys::getKey(domain_,"mass_source"),name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_,"mass_source"),name_)->set_initialized();
  S->GetFieldData(Keys::getKey(domain_,"total_energy_source"),name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_,"total_energy_source"),name_)->set_initialized();

  S->GetFieldData(Keys::getKey(domain_ss_,"mass_source"),name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_ss_,"mass_source"),name_)->set_initialized();
  S->GetFieldData(Keys::getKey(domain_ss_,"total_energy_source"),name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_ss_,"total_energy_source"),name_)->set_initialized();

  if (eval_derivatives_) {
    Key key_cond_temp = Keys::getDerivKey(Keys::getKey(domain_,"total_energy_source"), Keys::getKey(domain_,"temperature"));
    S->GetFieldData(key_cond_temp,name_)->PutScalar(0.);
    S->GetField(key_cond_temp,name_)->set_initialized();
  }

  // initialize diagnostics
  S->GetField(Keys::getKey(domain_,"albedo"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"evaporative_flux"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"qE_latent_heat"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"qE_sensible_heat"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"qE_lw_out"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"qE_melt_out"),name_)->set_initialized();
  if (!longwave_input_) {
    S->GetField(Keys::getKey(domain_,"incoming_longwave_radiation"),name_)->set_initialized();
  }
}


// computes the non-linear functional g = g(t,u,udot)
void
SurfaceBalanceImplicit::Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                            Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  Teuchos::OSTab tab = vo_->getOSTab();
  double dt = t_new - t_old;
  double T_eps = 1.e-10;

  bool debug = false;
  Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
  int rank = mesh_->get_comm()->MyPID();

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << dt << std::endl;
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("old snow depth");
    vecs.push_back(u_old->Data().ptr());
    vnames.push_back("solver snow depth");
    vecs.push_back(u_new->Data().ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  // pull residual vector
  Epetra_MultiVector& res = *g->Data()->ViewComponent("cell",false);

  // pull old snow data
  const Epetra_MultiVector& snow_depth_old = *u_old->Data()->ViewComponent("cell",false);

  const Epetra_MultiVector& snow_age_old = *S_inter_->GetFieldData(Keys::getKey(domain_,"snow_age"))
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_dens_old = *S_inter_->GetFieldData(Keys::getKey(domain_,"snow_density"))
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& stored_SWE_old = *S_inter_->GetFieldData(Keys::getKey(domain_,"stored_SWE"))
      ->ViewComponent("cell",false);

  // pull current snow data
  const Epetra_MultiVector& snow_depth_new = *u_new->Data()->ViewComponent("cell",false);
  CompositeVector snow_depth_new_computed(*u_new->Data());
  Epetra_MultiVector& snow_depth_new_computed_c = *snow_depth_new_computed.ViewComponent("cell",false);

  Epetra_MultiVector& snow_temp_new = *S_next_->GetFieldData(Keys::getKey(domain_,"snow_temperature"), name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& snow_age_new = *S_next_->GetFieldData(Keys::getKey(domain_,"snow_age"), name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& snow_dens_new = *S_next_->GetFieldData(Keys::getKey(domain_,"snow_density"), name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& stored_SWE_new = *S_next_->GetFieldData(Keys::getKey(domain_,"stored_SWE"), name_)
      ->ViewComponent("cell",false);

  // pull diagnostics
  Epetra_MultiVector& albedo = *S_next_->GetFieldData(Keys::getKey(domain_,"albedo"),name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& evaporative_flux = *S_next_->GetFieldData(Keys::getKey(domain_,"evaporative_flux"),name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& qE_latent_heat = *S_next_->GetFieldData(Keys::getKey(domain_,"qE_latent_heat"),name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& qE_sensible_heat = *S_next_->GetFieldData(Keys::getKey(domain_,"qE_sensible_heat"),name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& qE_lw_out = *S_next_->GetFieldData(Keys::getKey(domain_,"qE_lw_out"),name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& qE_melt_out = *S_next_->GetFieldData(Keys::getKey(domain_,"qE_melt_out"),name_)
      ->ViewComponent("cell",false);

  // pull ATS data
  S_next_->GetFieldEvaluator(Keys::getKey(domain_,"temperature"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_temp =
    *S_next_->GetFieldData(Keys::getKey(domain_,"temperature"))->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator(Keys::getKey(domain_,"pressure"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_pres =
    *S_next_->GetFieldData(Keys::getKey(domain_,"pressure"))->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator(Keys::getKey(domain_,"ponded_depth"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& ponded_depth =
    *S_next_->GetFieldData(Keys::getKey(domain_,"ponded_depth"))->ViewComponent("cell", false);

  // cannot work with current implementation because WRMs are hard coded into
  //  flow PKs.  Design fail.  Future versions will need this call.
  //  S_next_->GetFieldEvaluator(Keys::getKey(domain_ss_,"saturation_gas"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& saturation_gas =
     *S_next_->GetFieldData(Keys::getKey(domain_ss_,"saturation_gas"))->ViewComponent("cell", false);

   S_next_->GetFieldEvaluator(Keys::getKey(domain_,"unfrozen_fraction"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& unfrozen_fraction =
    *S_next_->GetFieldData(Keys::getKey(domain_,"unfrozen_fraction"))->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator(Keys::getKey(domain_ss_,"porosity"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& porosity =
    *S_next_->GetFieldData(Keys::getKey(domain_ss_,"porosity"))->ViewComponent("cell", false);


  // pull Met data
  S_next_->GetFieldEvaluator(Keys::getKey(domain_,"air_temperature"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& air_temp =
    *S_next_->GetFieldData(Keys::getKey(domain_,"air_temperature"))->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator(sw_incoming_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& incoming_shortwave =
    *S_next_->GetFieldData(sw_incoming_key_)->ViewComponent("cell", false);

 Teuchos::RCP<const Epetra_MultiVector> incoming_longwave = Teuchos::null;
 Teuchos::RCP<Epetra_MultiVector> incoming_longwave_diag = Teuchos::null;
  if (longwave_input_) {
    S_next_->GetFieldEvaluator(Keys::getKey(domain_,"incoming_longwave_radiation"))->HasFieldChanged(S_next_.ptr(), name_);
    incoming_longwave =
        S_next_->GetFieldData(Keys::getKey(domain_,"incoming_longwave_radiation"))->ViewComponent("cell", false);
  } else {
    incoming_longwave_diag =
        S_next_->GetFieldData(Keys::getKey(domain_,"incoming_longwave_radiation"), name_)->ViewComponent("cell", false);
  }

  S_next_->GetFieldEvaluator(Keys::getKey(domain_,"relative_humidity"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& relative_humidity =
    *S_next_->GetFieldData(Keys::getKey(domain_,"relative_humidity"))->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator(Keys::getKey(domain_,"wind_speed"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& wind_speed =
    *S_next_->GetFieldData(Keys::getKey(domain_,"wind_speed"))->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator(Keys::getKey(domain_,"precipitation_rain"))->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& precip_rain =
    *S_next_->GetFieldData(Keys::getKey(domain_,"precipitation_rain"))->ViewComponent("cell", false);

  // snow precip need not be updated each iteration
  bool update = false;
  if (implicit_snow_) {
    S_next_->GetFieldEvaluator(Keys::getKey(domain_,"precipitation_snow"))->HasFieldChanged(S_next_.ptr(), name_);
  } else {
    update = S_inter_->GetFieldEvaluator(Keys::getKey(domain_,"precipitation_snow"))->HasFieldChanged(S_inter_.ptr(), name_);
  }
  const Epetra_MultiVector& precip_snow = implicit_snow_ ?
    *S_next_->GetFieldData(Keys::getKey(domain_,"precipitation_snow"))->ViewComponent("cell", false) :
    *S_inter_->GetFieldData(Keys::getKey(domain_,"precipitation_snow"))->ViewComponent("cell", false);

  // pull additional primary variable data
  Epetra_MultiVector& surf_water_flux =
    *S_next_->GetFieldData(Keys::getKey(domain_,"mass_source"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& surf_energy_flux =
    *S_next_->GetFieldData(Keys::getKey(domain_,"total_energy_source"), name_)->ViewComponent("cell", false);
  Teuchos::RCP<Epetra_MultiVector> dsurf_energy_flux_dT;
  if (eval_derivatives_) {
    Key key_cond_temp = Keys::getDerivKey(Keys::getKey(domain_,"total_energy_source"), Keys::getKey(domain_,"temperature"));
    dsurf_energy_flux_dT = S_next_->GetFieldData(key_cond_temp, name_)
      ->ViewComponent("cell", false);
  }
  Epetra_MultiVector& ss_water_flux =
    *S_next_->GetFieldData(Keys::getKey(domain_ss_,"mass_source"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& ss_energy_flux =
    *S_next_->GetFieldData(Keys::getKey(domain_ss_,"total_energy_source"), name_)->ViewComponent("cell", false);

  // loop over each cell
  unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (unsigned int c=0; c!=ncells; ++c) {
    dcvo = Teuchos::null;
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) dcvo = db_->GetVerboseObject(c, rank);
    Teuchos::OSTab dctab = dcvo == Teuchos::null ? vo_->getOSTab() : dcvo->getOSTab();

    SEBPhysics::ModelParams seb_params;
    SEBPhysics::GroundProperties seb_surf;
    SEBPhysics::SnowProperties seb_snow;
    SEBPhysics::MetData seb_met;
    SEBPhysics::SurfaceParams seb_surf_pars;

    double snow_depth = snow_depth_old[0][c];
    double snow_dens = snow_dens_old[0][c];
    double swe = stored_SWE_old[0][c];
    double snow_age = snow_age_old[0][c];
    if (!implicit_snow_) {
      double Ps = std::max(precip_snow[0][c], 0.);
      swe += dt * Ps;
      snow_depth += dt * Ps * seb_params.density_water / seb_params.density_freshsnow;
      snow_age = swe > 0. ? (snow_age * stored_SWE_old[0][c]) / swe : 0.; // age of new snow is 0 here, as increment happens internally
      snow_dens = snow_depth > 0. ? swe * seb_params.density_water / snow_depth : seb_params.density_freshsnow;
    }

    if (snow_depth >= snow_ground_trans_ ||
        snow_depth < min_snow_trans_) {
      // Evaluate the model as usual
      // Initialize the SEB object

      // -- ground properties
      seb_surf.temp = surf_temp[0][c];
      seb_surf.pressure = surf_pres[0][c];
      seb_surf.roughness = roughness_bare_ground_;

      AmanziMesh::Entity_ID subsurf_f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
      AmanziMesh::Entity_ID_List cells;
      subsurf_mesh_->face_get_cells(subsurf_f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);
      seb_surf.saturation_gas = saturation_gas[0][cells[0]];
      seb_surf.density_w = seb_params.density_water; // NOTE: could update this to use true density! --etc
      seb_surf.dz = subsurf_mesh_->cell_volume(cells[0]) / subsurf_mesh_->face_area(subsurf_f) / 2.0;
      ASSERT(seb_surf.dz > 0.);
      
      SEBPhysics::Partition al_part = SEBPhysics::Partitioner()
          .CalcPartition(0., ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb_surf.albedo = al_part.Interpolate(1., seb_surf_pars.a_water, seb_surf_pars.a_ice, seb_surf_pars.a_tundra);

      SEBPhysics::Partition other_part = SEBPhysics::Partitioner(0.02, 0.02)
          .CalcPartition(0., ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb_surf.emissivity = other_part.Interpolate(1., seb_surf_pars.e_water, seb_surf_pars.e_ice, seb_surf_pars.e_tundra);
      seb_surf.porosity = other_part.Interpolate(1., 1., 1., porosity[0][cells[0]]);
      seb_surf.ponded_depth = ponded_depth[0][c];
      
      // -- snow properties
      seb_snow.height = snow_depth < min_snow_trans_ ? 0. : snow_depth;
      seb_snow.density = snow_dens;
      seb_snow.age = snow_age;
      seb_snow.SWE = swe;
      seb_snow.temp = 273.15;
      seb_snow.albedo = SEBPhysics::CalcAlbedoSnow(snow_dens);
      seb_snow.emissivity = seb_surf_pars.e_snow;
      seb_snow.roughness = roughness_snow_covered_ground_;

      // -- met data
      seb_met.Z_Us = wind_speed_ref_ht_;
      seb_met.Us = std::max(wind_speed[0][c], min_wind_speed_);
      seb_met.QswIn = incoming_shortwave[0][c];
      seb_met.Ps = implicit_snow_ ? std::max(precip_snow[0][c],0.) : 0.; // protect against wayward snow distribution models
      seb_met.Pr = precip_rain[0][c];
      seb_met.air_temp = air_temp[0][c];
      seb_met.relative_humidity = relative_humidity[0][c];
     
      if (longwave_input_) {
        seb_met.QlwIn = (*incoming_longwave)[0][c];
      } else {     
        seb_met.QlwIn = SEBPhysics::CalcIncomingLongwave(seb_met.air_temp, seb_met.relative_humidity, seb_params.stephB);
        (*incoming_longwave_diag)[0][c] = seb_met.QlwIn;
      }

      // Run the model
      auto result = SEBPhysics::CalculateSurfaceBalance(dt, seb_surf, seb_snow, seb_met, seb_params, false, dcvo);
      SEBPhysics::SnowProperties& seb_snow_new = std::get<0>(result);
      SEBPhysics::EnergyBalance& seb_eb = std::get<1>(result);
      SEBPhysics::MassBalance& seb_mb = std::get<2>(result);
      SEBPhysics::FluxBalance& seb_flux = std::get<3>(result);

      // Evaluate the residual
      res[0][c] = snow_depth_new[0][c] - seb_snow_new.height;
      snow_depth_new_computed_c[0][c] = seb_snow_new.height;

      // Pull the output
      // -- fluxes
      //    -- surface mass = Pr + Melt - evap if ponded and no snow
      
      //    -- surface energy = diffusion + change due to evap
      surf_energy_flux[0][c] = seb_flux.E_surf;
      surf_water_flux[0][c] = seb_flux.M_surf;

      // -- vapor flux to cells
      //     surface vapor flux is treated as a volumetric source for the subsurface.
      // surface mass sources are in m^3 water / (m^2 s)
      // subsurface mass sources are in mol water / (m^3 s)
      ss_water_flux[0][cells[0]] = seb_flux.M_sub * mesh_->cell_volume(c)
                                * seb_surf.density_w / 0.0180153
                                / subsurf_mesh_->cell_volume(cells[0]);
      // surface energy sources are in W / m^2 
      // subsurface mass sources are in W / m^3
      ss_energy_flux[0][cells[0]] = seb_flux.E_sub * mesh_->cell_volume(c)
                                 / subsurf_mesh_->cell_volume(cells[0]);

      // -- snow properties
      snow_age_new[0][c] = seb_snow_new.age;
      snow_dens_new[0][c] = seb_snow_new.density;
      snow_temp_new[0][c] = seb_snow_new.temp;
      stored_SWE_new[0][c] = seb_snow_new.SWE;
      
      // -- diagnostics
      albedo[0][c] = seb_snow.height > snow_ground_trans_ ? seb_snow.albedo : seb_surf.albedo;
      evaporative_flux[0][c] = -seb_mb.Me;
      qE_latent_heat[0][c] = seb_eb.fQe;
      qE_sensible_heat[0][c] = seb_eb.fQh;
      qE_lw_out[0][c] = seb_eb.fQlwOut;
      qE_melt_out[0][c] = seb_eb.fQm;

      if (eval_derivatives_) {
        // evaluate FD derivative of energy flux wrt surface temperature
        SEBPhysics::GroundProperties seb_surf2(seb_surf);
        seb_surf2.temp += T_eps;
        // for now ignore the effect on unfrozen fraction, and therefore on albedo and emissivity
        auto result2 = SEBPhysics::CalculateSurfaceBalance(dt, seb_surf2, seb_snow, seb_met, seb_params, false, dcvo);
        SEBPhysics::FluxBalance& seb_flux2 = std::get<3>(result2);
        (*dsurf_energy_flux_dT)[0][c] = (seb_flux2.E_surf - seb_flux2.E_surf) / T_eps;
      }

    } else {
      // Evaluate the model twice -- once as bare ground, once with snow, using
      // an area-averaged subgrid model to smooth between the two end-members.
      // The area-weighting parameter is theta:
      double theta = snow_depth / snow_ground_trans_;

      // Evaluate the model as usual
      // -- ground properties
      seb_surf.temp = surf_temp[0][c];
      seb_surf.pressure = surf_pres[0][c];
      seb_surf.roughness = roughness_bare_ground_;

      AmanziMesh::Entity_ID subsurf_f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
      AmanziMesh::Entity_ID_List cells;
      subsurf_mesh_->face_get_cells(subsurf_f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);
      seb_surf.saturation_gas = saturation_gas[0][cells[0]];
      seb_surf.density_w = seb_params.density_water; // NOTE: could update this to use true density! --etc
      seb_surf.dz = subsurf_mesh_->cell_volume(cells[0]) / subsurf_mesh_->face_area(subsurf_f) / 2.0;
      ASSERT(seb_surf.dz > 0.);
      
      SEBPhysics::Partition al_part = SEBPhysics::Partitioner()
          .CalcPartition(0., ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb_surf.albedo = al_part.Interpolate(1., seb_surf_pars.a_water, seb_surf_pars.a_ice, seb_surf_pars.a_tundra);

      SEBPhysics::Partition other_part = SEBPhysics::Partitioner(0.02, 0.02)
          .CalcPartition(0., ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb_surf.emissivity = other_part.Interpolate(1., seb_surf_pars.e_water, seb_surf_pars.e_ice, seb_surf_pars.e_tundra);
      seb_surf.porosity = other_part.Interpolate(1., 1., 1., porosity[0][cells[0]]);
      seb_surf.ponded_depth = ponded_depth[0][c];
      
      // -- snow properties
      seb_snow.height = snow_ground_trans_;
      seb_snow.density = snow_dens;
      seb_snow.age = snow_age;
      seb_snow.SWE = swe / theta;
      seb_snow.temp = 273.15;
      seb_snow.albedo = SEBPhysics::CalcAlbedoSnow(snow_dens);
      seb_snow.emissivity = seb_surf_pars.e_snow;
      seb_snow.roughness = roughness_snow_covered_ground_;

      SEBPhysics::SnowProperties seb_snow_bare;
      seb_snow_bare.height = 0.;
      seb_snow_bare.density = seb_params.density_freshsnow;
      seb_snow_bare.age = 0.;
      seb_snow_bare.SWE = 0.;
      seb_snow_bare.temp = 273.15;
      seb_snow_bare.albedo = 1.;
      seb_snow_bare.emissivity = 1.;
      seb_snow_bare.roughness = 0.;
      
      // -- met data
      seb_met.Z_Us = wind_speed_ref_ht_;
      seb_met.Us = std::max(wind_speed[0][c], min_wind_speed_);
      seb_met.QswIn = incoming_shortwave[0][c];
      seb_met.Ps = implicit_snow_ ? std::max(precip_snow[0][c],0.) : 0.; // protect against wayward snow distribution models
      seb_met.Pr = precip_rain[0][c];
      seb_met.air_temp = air_temp[0][c];
      seb_met.relative_humidity = relative_humidity[0][c];
     
     if (longwave_input_) {
       seb_met.QlwIn = (*incoming_longwave)[0][c];
      } else {     
       seb_met.QlwIn = SEBPhysics::CalcIncomingLongwave(seb_met.air_temp, seb_met.relative_humidity, seb_params.stephB);
       (*incoming_longwave_diag)[0][c] = seb_met.QlwIn;
     }

     // Run the model
     auto result = SEBPhysics::CalculateSurfaceBalance(dt, seb_surf, seb_snow, seb_met, seb_params, false, dcvo);
     auto result_bare = SEBPhysics::CalculateSurfaceBalance(dt, seb_surf, seb_snow_bare, seb_met, seb_params, false, dcvo);

     SEBPhysics::SnowProperties& seb_snow_new = std::get<0>(result);
     SEBPhysics::EnergyBalance& seb_eb = std::get<1>(result);
     SEBPhysics::MassBalance& seb_mb = std::get<2>(result);
     SEBPhysics::FluxBalance& seb_flux = std::get<3>(result);

     SEBPhysics::SnowProperties& bare_snow_new = std::get<0>(result_bare);
     SEBPhysics::EnergyBalance& bare_eb = std::get<1>(result_bare);
     SEBPhysics::MassBalance& bare_mb = std::get<2>(result_bare);
     SEBPhysics::FluxBalance& bare_flux = std::get<3>(result_bare);

     // Evaluate the residual
     double snow_depth_new_tmp = theta * seb_snow_new.height + (1-theta) * bare_snow_new.height;
     res[0][c] = snow_depth_new[0][c] - snow_depth_new_tmp;
     snow_depth_new_computed_c[0][c] = snow_depth_new_tmp;

     // Pull the output
     // -- fluxes
     surf_energy_flux[0][c] = theta * seb_flux.E_surf + (1-theta) * bare_flux.E_surf;
     surf_water_flux[0][c] = theta * seb_flux.M_surf + (1-theta) * bare_flux.M_surf;
     // -- vapor flux to cells
     //     surface vapor flux is treated as a volumetric source for the subsurface.
     // surface mass sources are in m^3 water / (m^2 s)
     // subsurface mass sources are in mol water / (m^3 s)
     double mean_flux = theta * seb_flux.M_sub + (1-theta) * bare_flux.M_sub;
     ss_water_flux[0][cells[0]] = mean_flux * mesh_->cell_volume(c)
                               * seb_surf.density_w / 0.0180153
                               / subsurf_mesh_->cell_volume(cells[0]);
     // surface energy sources are in W / m^2
     // subsurface energy sources are in W / m^3
     double mean_flux_E = theta * seb_flux.E_sub + (1-theta) * bare_flux.E_sub;
     ss_energy_flux[0][cells[0]] = mean_flux_E * mesh_->cell_volume(c)
                               / subsurf_mesh_->cell_volume(cells[0]);

     // -- snow properties: SWE averaged
     double total_swe = theta * seb_snow_new.height * seb_snow_new.density / seb_surf.density_w
                        + (1-theta) * bare_snow_new.height * bare_snow_new.density / seb_surf.density_w;
     if (snow_depth_new_tmp > 0.) {
       snow_age_new[0][c] = (theta * seb_snow_new.height * seb_snow_new.age
                             + (1-theta) * bare_snow_new.height * bare_snow_new.age) / snow_depth_new_tmp;
       snow_dens_new[0][c] = total_swe * seb_surf.density_w / snow_depth_new_tmp;
       stored_SWE_new[0][c] = total_swe;
     } else {
       snow_age_new[0][c] = 0.;
       snow_dens_new[0][0] = seb_snow_new.density;
       stored_SWE_new[0][c] = total_swe;
     }
     snow_temp_new[0][c] = seb_snow_new.temp;

     // -- diagnostics
     albedo[0][c] = theta * seb_snow_new.albedo + (1-theta) * seb_surf.albedo;
     evaporative_flux[0][c] = -(theta * seb_mb.Me + (1-theta) * bare_mb.Me);
     qE_latent_heat[0][c] = (theta * seb_eb.fQe + (1-theta) * bare_eb.fQe);
     qE_sensible_heat[0][c] = (theta * seb_eb.fQh + (1-theta) * bare_eb.fQh);
     qE_lw_out[0][c] = theta * seb_eb.fQlwOut + (1-theta) * bare_eb.fQlwOut;
     qE_melt_out[0][c] = theta * seb_eb.fQm + (1-theta) * bare_eb.fQm;

     // Evaluate derivatives, if requested
     if (eval_derivatives_) {
       // evaluate FD derivative of energy flux wrt surface temperature
       SEBPhysics::GroundProperties seb_surf2(seb_surf);
       seb_surf2.temp += T_eps;
       // for now ignore the effect on unfrozen fraction, and therefore on albedo and emissivity
       auto result2 = SEBPhysics::CalculateSurfaceBalance(dt, seb_surf2, seb_snow, seb_met, seb_params, false, dcvo);
       SEBPhysics::EnergyBalance& seb_eb2 = std::get<1>(result2);

       auto result_bare2 = SEBPhysics::CalculateSurfaceBalance(dt, seb_surf2, seb_snow_bare, seb_met, seb_params, false, dcvo);
       SEBPhysics::EnergyBalance& bare_eb2 = std::get<1>(result_bare2);

       double eflux2 = (theta * seb_eb2.fQc + (1-theta) * bare_eb2.fQc); // MJ
       (*dsurf_energy_flux_dT)[0][c] = (eflux2 - surf_energy_flux[0][c]) / T_eps;
     }
    }
  }  // END CELL LOOP ###############################

  // debug
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("computed snow depth");
    Teuchos::Ptr<CompositeVector> snow_depth_new_computed_p(&snow_depth_new_computed);
    vecs.push_back(snow_depth_new_computed_p);
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.push_back("air_temp"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"air_temperature")).ptr());
    vnames.push_back("rel_hum"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"relative_humidity")).ptr());
    vnames.push_back("precip_rain"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"precipitation_rain")).ptr());
    vnames.push_back("precip_snow");
    if (implicit_snow_) {
      vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"precipitation_snow")).ptr());
    } else {
      vecs.push_back(S_inter_->GetFieldData(Keys::getKey(domain_,"precipitation_snow")).ptr());
    }      
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();
    
    vnames.push_back("p_ground"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"pressure")).ptr());
    vnames.push_back("ponded_depth"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"ponded_depth")).ptr());
    vnames.push_back("snow_depth"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"snow_depth")).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("T_ground"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"temperature")).ptr());
    vnames.push_back("snow_temp"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"snow_temperature")).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("inc shortwave radiation"); 
    vecs.push_back(S_next_->GetFieldData(sw_incoming_key_).ptr());
    vnames.push_back("inc longwave radiation");
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_, "incoming_longwave_radiation")).ptr());
    vnames.push_back("inc latent heat"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"qE_latent_heat")).ptr());
    vnames.push_back("inc sensible heat"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"qE_sensible_heat")).ptr());
    vnames.push_back("out longwave radiation"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"qE_lw_out")).ptr());
    vnames.push_back("out conducted energy"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"total_energy_source")).ptr());
    vnames.push_back("out melting energy"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"qE_melt_out")).ptr());

    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();
    
    vnames.push_back("water_source"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"mass_source")).ptr());
    vnames.push_back("evap flux"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"evaporative_flux")).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("res(snow_diff)"); 
    vecs.push_back(g->Data().ptr());
    
    db_->WriteVectors(vnames, vecs, true);
  }

  // scale to MJ for ATS units
  surf_energy_flux.Scale(1.e-6);
  if (dsurf_energy_flux_dT.get()) dsurf_energy_flux_dT->Scale(1.e-6);
  ss_energy_flux.Scale(1.e-6);
  
  // mark our other primary variables as changed
  pvfe_esource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wsource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_w_sub_source_->SetFieldAsChanged(S_next_.ptr());
  pvfe_e_sub_source_->SetFieldAsChanged(S_next_.ptr());
  
}


// applies preconditioner to u and returns the result in Pu
int SurfaceBalanceImplicit::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  *Pu = *u;
  
  return 0;
}


// updates the preconditioner
void
SurfaceBalanceImplicit::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {}


// error monitor, inf norm is good, this is relative to 1m snow pack
double
SurfaceBalanceImplicit::ErrorNorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();
  double err;
  du->NormInf(&err);
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << "ENorm (cells) = " << err << std::endl;
  }
  return err;
}


bool
SurfaceBalanceImplicit::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  Epetra_MultiVector& u_vec = *u->Data()->ViewComponent("cell",false);
  unsigned int ncells = u_vec.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    u_vec[0][c] = std::max(0., u_vec[0][c]);
  }

  if (modify_predictor_advance_) {
    Teuchos::RCP<TreeVector> res = Teuchos::rcp(new TreeVector(*u));
    Teuchos::RCP<TreeVector> u0_nc = Teuchos::rcp_const_cast<TreeVector>(u0);
    Functional(S_next_->time()-h, S_next_->time(), u0_nc, u, res);
  }

  return true;
}




} // namespace
} // namespace
