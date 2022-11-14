/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (coonet @ ornl.gov)

 ------------------------------------------------------------------------- */

//! SEBTwoComponentEvaluator: evaluates the Surface Energy Balance model on subgrid units.

/*!

Sets up a collection of patches, for portions of the column covered in snow,
ponded water, and vegetated/bare ground.  The surface energy balance on these
area weighted patches are individually calculated then averaged to form the
total quantities.  All down- and up-scaling of relevant quantities are done
through the area weighting, which is calculated by a minimum threshold in snow
and a depression depth/geometry-based approach for water.  All snow is assumed
to first cover water (likely ice), then cover land, as both water and snow
prefer low-lying depressions due to gravity- and wind-driven redistributions,
respectively.

*

*/

#include "seb_twocomponent_evaluator.hh"
#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

SEBTwoComponentEvaluator::SEBTwoComponentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist),
    plist_(plist),
    compatible_(false),
    model_1p1_(false)
{
  // determine the domain
  Key domain = Keys::getDomain(Keys::cleanPListName(plist_.name()));
  Key dtype = Keys::guessDomainType(domain);
  if (dtype == "surface") {
    domain_ = domain;
    domain_ss_ = Keys::readDomainHint(plist_, domain_, "surface", "subsurface");
    domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");
  } else if (dtype == "domain") {
    domain_ss_ = domain;
    domain_ = Keys::readDomainHint(plist_, domain_ss_, "domain", "surface");
    domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");
  } else if (dtype == "snow") {
    domain_snow_ = domain;
    domain_ = Keys::readDomainHint(plist_, domain_snow_, "snow", "surface");
    domain_ss_ = Keys::readDomainHint(plist_, domain_snow_, "snow", "subsurface");
  } else {
    domain_snow_ = plist.get<std::string>("snow domain name", domain_snow_);
    domain_ = plist.get<std::string>("surface domain name", domain_);
    domain_ss_ = plist.get<std::string>("subsurface domain name", domain_ss_);
  }

  // my keys
  // -- sources
  water_source_key_ = Keys::readKey(plist, domain_, "surface water source", "water_source");
  my_keys_.push_back(water_source_key_);
  energy_source_key_ = Keys::readKey(plist, domain_, "surface energy source", "total_energy_source");
  my_keys_.push_back(energy_source_key_);
  ss_water_source_key_ = Keys::readKey(plist, domain_ss_, "subsurface water source", "water_source");
  my_keys_.push_back(ss_water_source_key_);
  ss_energy_source_key_ = Keys::readKey(plist, domain_ss_, "subsurface energy source", "total_energy_source");
  my_keys_.push_back(ss_energy_source_key_);
  snow_source_key_ = Keys::readKey(plist, domain_snow_, "snow mass source - sink", "source_sink");
  my_keys_.push_back(snow_source_key_);

  new_snow_key_ = Keys::readKey(plist, domain_snow_, "new snow source", "source");
  my_keys_.push_back(new_snow_key_);

  // old model from version 1.1, deprecated
  model_1p1_ = plist.get<bool>("use model from ATS 1.1", false);

  // diagnostics and debugging
  diagnostics_ = plist.get<bool>("save diagnostic data", false);
  if (diagnostics_) {
    // -- diagnostics
    albedo_key_ = Keys::readKey(plist, domain_, "albedo", "albedo");
    my_keys_.push_back(albedo_key_);
    melt_key_ = Keys::readKey(plist, domain_snow_, "snowmelt", "melt");
    my_keys_.push_back(melt_key_);
    evap_key_ = Keys::readKey(plist, domain_, "evaporation", "evaporative_flux");
    my_keys_.push_back(evap_key_);
    snow_temp_key_ = Keys::readKey(plist, domain_snow_, "snow temperature", "temperature");
    my_keys_.push_back(snow_temp_key_);
    qE_sh_key_ = Keys::readKey(plist, domain_, "sensible heat flux", "qE_sensible_heat");
    my_keys_.push_back(qE_sh_key_);
    qE_lh_key_ = Keys::readKey(plist, domain_, "latent heat of evaporation", "qE_latent_heat");
    my_keys_.push_back(qE_lh_key_);
    qE_sm_key_ = Keys::readKey(plist, domain_, "latent heat of snowmelt", "qE_snowmelt");
    my_keys_.push_back(qE_sm_key_);
    qE_lw_out_key_ = Keys::readKey(plist, domain_, "outgoing longwave radiation", "qE_lw_out");
    my_keys_.push_back(qE_lw_out_key_);
    qE_cond_key_ = Keys::readKey(plist, domain_, "conducted energy flux", "qE_conducted");
    my_keys_.push_back(qE_cond_key_);
  }

  // dependencies
  // -- met data
  met_sw_key_ = Keys::readKey(plist, domain_,"incoming shortwave radiation", "incoming_shortwave_radiation");
  dependencies_.insert(met_sw_key_);
  met_lw_key_ = Keys::readKey(plist, domain_,"incoming longwave radiation", "incoming_longwave_radiation");
  dependencies_.insert(met_lw_key_);
  met_air_temp_key_ = Keys::readKey(plist, domain_,"air temperature", "air_temperature");
  dependencies_.insert(met_air_temp_key_);
  met_rel_hum_key_ = Keys::readKey(plist, domain_,"relative humidity", "relative_humidity");
  dependencies_.insert(met_rel_hum_key_);
  met_wind_speed_key_ = Keys::readKey(plist, domain_,"wind speed", "wind_speed");
  dependencies_.insert(met_wind_speed_key_);
  met_prain_key_ = Keys::readKey(plist, domain_,"precipitation rain", "precipitation_rain");
  dependencies_.insert(met_prain_key_);
  met_psnow_key_ = Keys::readKey(plist, domain_snow_,"precipitation snow", "precipitation");
  dependencies_.insert(met_psnow_key_);

  // -- snow properties
  snow_depth_key_ = Keys::readKey(plist, domain_snow_, "snow depth", "depth");
  dependencies_.insert(snow_depth_key_);
  snow_dens_key_ = Keys::readKey(plist, domain_snow_, "snow density", "density");
  dependencies_.insert(snow_dens_key_);
  snow_death_rate_key_ = Keys::readKey(plist, domain_snow_, "snow death rate", "death_rate");
  dependencies_.insert(snow_death_rate_key_);

  // -- skin properties
  mol_dens_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(mol_dens_key_);
  mass_dens_key_ = Keys::readKey(plist, domain_, "mass density liquid", "mass_density_liquid");
  dependencies_.insert(mass_dens_key_);
  ponded_depth_key_ = Keys::readKey(plist, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(ponded_depth_key_);
  unfrozen_fraction_key_ = Keys::readKey(plist, domain_, "unfrozen fraction", "unfrozen_fraction");
  dependencies_.insert(unfrozen_fraction_key_);
  sg_albedo_key_ = Keys::readKey(plist, domain_, "albedos", "albedos");
  dependencies_.insert(sg_albedo_key_);
  sg_emissivity_key_ = Keys::readKey(plist, domain_, "emissivities", "emissivities");
  dependencies_.insert(sg_emissivity_key_);
  area_frac_key_ = Keys::readKey(plist, domain_, "area fractions", "area_fractions");
  // explicitly excluded to allow snow_death algorithm to work, see #8
  dependencies_.insert(area_frac_key_);

  surf_temp_key_ = Keys::readKey(plist, domain_, "temperature", "temperature");
  dependencies_.insert(surf_temp_key_);
  surf_pres_key_ = Keys::readKey(plist, domain_, "pressure", "pressure");
  dependencies_.insert(surf_pres_key_);

  // -- subsurface properties for evaporating bare soil
  sat_gas_key_ = Keys::readKey(plist, domain_ss_, "gas saturation", "saturation_gas");
  dependencies_.insert(sat_gas_key_);
  poro_key_ = Keys::readKey(plist, domain_ss_, "porosity", "porosity");
  dependencies_.insert(poro_key_);
  ss_pres_key_ = Keys::readKey(plist, domain_ss_, "subsurface pressure", "pressure");
  dependencies_.insert(ss_pres_key_);

  // parameters
  min_rel_hum_ = plist.get<double>("minimum relative humidity [-]", 0.1);
  min_wind_speed_ = plist.get<double>("minimum wind speed [m s^-1]", 1.0);
  wind_speed_ref_ht_ = plist.get<double>("wind speed reference height [m]", 2.0);
}

void
SEBTwoComponentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                             const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  const Relations::ModelParams params(plist_);
  double snow_eps = 1.e-5;

  // collect met data
  const auto& qSW_in = *S->GetFieldData(met_sw_key_)->ViewComponent("cell",false);
  const auto& qLW_in = *S->GetFieldData(met_lw_key_)->ViewComponent("cell",false);
  const auto& air_temp = *S->GetFieldData(met_air_temp_key_)->ViewComponent("cell",false);
  const auto& rel_hum = *S->GetFieldData(met_rel_hum_key_)->ViewComponent("cell",false);
  const auto& wind_speed = *S->GetFieldData(met_wind_speed_key_)->ViewComponent("cell",false);
  const auto& Prain = *S->GetFieldData(met_prain_key_)->ViewComponent("cell",false);
  const auto& Psnow = *S->GetFieldData(met_psnow_key_)->ViewComponent("cell",false);

  // collect snow properties
  const auto& snow_depth = *S->GetFieldData(snow_depth_key_)->ViewComponent("cell",false);
  const auto& snow_dens = *S->GetFieldData(snow_dens_key_)->ViewComponent("cell",false);
  const auto& snow_death_rate = *S->GetFieldData(snow_death_rate_key_)->ViewComponent("cell",false);

  // collect skin properties
  const auto& mol_dens = *S->GetFieldData(mol_dens_key_)->ViewComponent("cell",false);
  const auto& mass_dens = *S->GetFieldData(mass_dens_key_)->ViewComponent("cell",false);
  const auto& ponded_depth = *S->GetFieldData(ponded_depth_key_)->ViewComponent("cell",false);
  const auto& unfrozen_fraction = *S->GetFieldData(unfrozen_fraction_key_)->ViewComponent("cell",false);
  const auto& sg_albedo = *S->GetFieldData(sg_albedo_key_)->ViewComponent("cell",false);
  const auto& emissivity = *S->GetFieldData(sg_emissivity_key_)->ViewComponent("cell",false);
  const auto& area_fracs = *S->GetFieldData(area_frac_key_)->ViewComponent("cell",false);
  const auto& surf_pres = *S->GetFieldData(surf_pres_key_)->ViewComponent("cell",false);
  const auto& surf_temp = *S->GetFieldData(surf_temp_key_)->ViewComponent("cell",false);

  // collect subsurface properties
  const auto& sat_gas = *S->GetFieldData(sat_gas_key_)->ViewComponent("cell",false);
  const auto& poro = *S->GetFieldData(poro_key_)->ViewComponent("cell",false);
  const auto& ss_pres = *S->GetFieldData(ss_pres_key_)->ViewComponent("cell",false);

  // collect output vecs
  auto& water_source = *results[0]->ViewComponent("cell",false);
  auto& energy_source = *results[1]->ViewComponent("cell",false);
  auto& ss_water_source = *results[2]->ViewComponent("cell",false);
  auto& ss_energy_source = *results[3]->ViewComponent("cell",false);
  auto& snow_source = *results[4]->ViewComponent("cell",false);
  auto& new_snow = *results[5]->ViewComponent("cell",false);
  water_source.PutScalar(0.);
  energy_source.PutScalar(0.);
  ss_water_source.PutScalar(0.);
  ss_energy_source.PutScalar(0.);
  snow_source.PutScalar(0.);
  new_snow.PutScalar(0.);

  const auto& mesh = *S->GetMesh(domain_);
  const auto& mesh_ss = *S->GetMesh(domain_ss_);

  Epetra_MultiVector *melt_rate(nullptr), *evap_rate(nullptr), *snow_temp(nullptr);
  Epetra_MultiVector *qE_sh(nullptr), *qE_lh(nullptr), *qE_sm(nullptr);
  Epetra_MultiVector *qE_lw_out(nullptr), *qE_cond(nullptr), *albedo(nullptr);
  if (diagnostics_) {
    albedo = S->GetFieldData(albedo_key_, albedo_key_)->ViewComponent("cell",false).get();
    albedo->PutScalar(0.);
    melt_rate = S->GetFieldData(melt_key_, melt_key_)->ViewComponent("cell",false).get();
    melt_rate->PutScalar(0.);
    evap_rate = S->GetFieldData(evap_key_, evap_key_)->ViewComponent("cell",false).get();
    evap_rate->PutScalar(0.);
    snow_temp = S->GetFieldData(snow_temp_key_, snow_temp_key_)->ViewComponent("cell",false).get();
    snow_temp->PutScalar(273.15);
    qE_sh = S->GetFieldData(qE_sh_key_, qE_sh_key_)->ViewComponent("cell",false).get();
    qE_sh->PutScalar(0.);
    qE_lh = S->GetFieldData(qE_lh_key_, qE_lh_key_)->ViewComponent("cell",false).get();
    qE_lh->PutScalar(0.);
    qE_sm = S->GetFieldData(qE_sm_key_, qE_sm_key_)->ViewComponent("cell",false).get();
    qE_sm->PutScalar(0.);
    qE_lw_out = S->GetFieldData(qE_lw_out_key_, qE_lw_out_key_)->ViewComponent("cell",false).get();
    qE_lw_out->PutScalar(0.);
    qE_cond = S->GetFieldData(qE_cond_key_, qE_cond_key_)->ViewComponent("cell",false).get();
    qE_cond->PutScalar(0.);
  }

  for (const auto& lc : land_cover_) {
    AmanziMesh::Entity_ID_List lc_ids;
    mesh.get_set_entities(lc.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

    for (auto c : lc_ids) {
      // get the top cell
      AmanziMesh::Entity_ID subsurf_f = mesh.entity_get_parent(AmanziMesh::CELL, c);
      AmanziMesh::Entity_ID_List cells;
      mesh_ss.face_get_cells(subsurf_f, AmanziMesh::Parallel_type::OWNED, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      // met data structure
      Relations::MetData met;
      met.Z_Us = wind_speed_ref_ht_;
      met.Us = std::max(wind_speed[0][c], min_wind_speed_);
      met.QswIn = qSW_in[0][c];
      met.QlwIn = qLW_in[0][c];
      met.air_temp = air_temp[0][c];
      met.relative_humidity = std::max(rel_hum[0][c], min_rel_hum_);
      met.Pr = Prain[0][c];

      std::cout << "surf_temp_key_ = " << surf_temp_key_ << std::endl;
      std::cout << "c = " << c << ", surf_temp[0][c] = " << surf_temp[0][c] << std::endl;

      // non-snow covered column
      if (area_fracs[0][c] > 0) {
        Relations::GroundProperties surf;
        surf.temp = surf_temp[0][c];
        surf.water_transition_depth = lc.second.water_transition_depth;
        if (ponded_depth[0][c] > lc.second.water_transition_depth) {
          surf.pressure = surf_pres[0][c];
          surf.porosity = 1.;
          surf.saturation_gas = 0.;
        } else {
          double factor = std::max(ponded_depth[0][c],0.) /
            lc.second.water_transition_depth;
          surf.pressure = factor*surf_pres[0][c] + (1-factor)*ss_pres[0][cells[0]];
          surf.porosity = factor + (1-factor)*poro[0][cells[0]];
          surf.saturation_gas = (1-factor)*sat_gas[0][cells[0]];
        }
        if (model_1p1_) surf.pressure = surf_pres[0][c];
        surf.ponded_depth = ponded_depth[0][c];
        surf.unfrozen_fraction = unfrozen_fraction[0][c];
        surf.roughness = lc.second.roughness_ground;
        if (model_1p1_) surf.density_w = 1000.;
        else surf.density_w = mass_dens[0][c];
        surf.dz = lc.second.dessicated_zone_thickness;
        surf.albedo = sg_albedo[0][c];
        surf.emissivity = emissivity[0][c];

        // must ensure that energy is put into melting snow precip, even if it
        // all melts so there is no snow column
        if (area_fracs[1][c] == 0.) {
          met.Ps = Psnow[0][c];
          surf.snow_death_rate = snow_death_rate[0][c]; // m H20 / s
        } else {
          met.Ps = 0.;
          surf.snow_death_rate = 0.;
        }

        // calculate the surface balance
        const Relations::EnergyBalance eb = Relations::UpdateEnergyBalanceWithoutSnow(surf, met, params);
        Relations::MassBalance mb = Relations::UpdateMassBalanceWithoutSnow(surf, params, eb);
        Relations::FluxBalance flux = Relations::UpdateFluxesWithoutSnow(surf, met, params, eb, mb, model_1p1_);

        // fQe, Me positive is condensation, water flux positive to surface
        water_source[0][c] += area_fracs[0][c] * flux.M_surf;
        energy_source[0][c] += area_fracs[0][c] * flux.E_surf * 1.e-6; // convert to MW/m^2

        double area_to_volume = mesh.cell_volume(c) / mesh_ss.cell_volume(cells[0]);
        double ss_water_source_l;
        if (model_1p1_) ss_water_source_l = flux.M_subsurf * area_to_volume * surf.density_w / 0.0180153; // convert from m/s to mol/m^3/s
        else ss_water_source_l = flux.M_subsurf * area_to_volume * mol_dens[0][c]; // convert from m/s to mol/m^3/s
        ss_water_source[0][cells[0]] += area_fracs[0][c] * ss_water_source_l;
        double ss_energy_source_l = flux.E_subsurf * area_to_volume * 1.e-6; // convert from W/m^2 to MW/m^3
        ss_energy_source[0][cells[0]] += area_fracs[0][c] * ss_energy_source_l;

        snow_source[0][c] += area_fracs[0][c] * flux.M_snow;
        new_snow[0][c] += area_fracs[0][c] * met.Ps;

        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "CELL " << c << " NO_SNOW"
                     << ": Ms = " << flux.M_surf << ", Es = " << flux.E_surf * 1.e-6
                     << ", Mss = " << ss_water_source_l << ", Ess = " << ss_energy_source_l
                     << ", Sn = " << flux.M_snow << std::endl;

        // diagnostics
        if (diagnostics_) {
          (*evap_rate)[0][c] -= area_fracs[0][c] * mb.Me;
          (*qE_sh)[0][c] += area_fracs[0][c] * eb.fQh;
          (*qE_lh)[0][c] += area_fracs[0][c] * eb.fQe;
          (*qE_lw_out)[0][c] += area_fracs[0][c] * eb.fQlwOut;
          (*qE_cond)[0][c] += area_fracs[0][c] * eb.fQc;
          (*albedo)[0][c] += area_fracs[0][c] * surf.albedo;

          if (area_fracs[1][c] == 0.) {
            (*qE_sm)[0][c] = eb.fQm;
            (*melt_rate)[0][c] = mb.Mm;
            (*snow_temp)[0][c] = 273.15;
          }
        }
      }

      // snow column
      if (area_fracs[1][c] > 0.) {
        Relations::GroundProperties surf;
        surf.temp = surf_temp[0][c];
        surf.pressure = surf_pres[0][c];
        surf.ponded_depth = ponded_depth[0][c];
        surf.porosity = 1.;
        surf.saturation_gas = 0.;
        surf.unfrozen_fraction = unfrozen_fraction[0][c];
        surf.roughness = lc.second.roughness_ground;
        if (model_1p1_) surf.density_w = 1000;
        else surf.density_w = mass_dens[0][c];
        surf.dz = lc.second.dessicated_zone_thickness;
        surf.albedo = sg_albedo[1][c];
        surf.emissivity = emissivity[1][c];
        surf.water_transition_depth = lc.second.water_transition_depth;

        met.Ps = Psnow[0][c] / area_fracs[1][c];

        Relations::SnowProperties snow;
        snow.height = snow_depth[0][c] / area_fracs[1][c]; // all snow on this patch
        AMANZI_ASSERT(snow.height >= lc.second.snow_transition_depth - 1.e-6);
        // area_fracs may have been set to 1 for snow depth < snow_ground_trans
        // due to min fractional area option in area_fractions evaluator.
        // Decreasing the tol by 1e-6 is about equivalent to a min fractional
        // area of 1e-5 (the default)
        snow.density = snow_dens[0][c];
        snow.albedo = surf.albedo;
        snow.emissivity = surf.emissivity;
        snow.roughness = lc.second.roughness_snow;

        const Relations::EnergyBalance eb = Relations::UpdateEnergyBalanceWithSnow(surf, met, params, snow);
        const Relations::MassBalance mb = Relations::UpdateMassBalanceWithSnow(surf, params, eb);
        Relations::FluxBalance flux = Relations::UpdateFluxesWithSnow(surf, met, params, snow, eb, mb);

        // fQe, Me positive is condensation, water flux positive to surface.  No
        // need for subsurf as there is snow present.
        water_source[0][c] += area_fracs[1][c] * flux.M_surf;
        energy_source[0][c] += area_fracs[1][c] * flux.E_surf * 1.e-6; // convert to MW/m^2 from W/m^2
        snow_source[0][c] += area_fracs[1][c] * flux.M_snow;
        new_snow[0][c] += std::max(met.Ps + mb.Me, 0.) * area_fracs[1][c];

        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "CELL " << c << " SNOW"
                     << ": Ms = " << flux.M_surf << ", Es = " << flux.E_surf * 1.e-6
                     << ", Mss = " << 0. << ", Ess = " << 0.
                     << ", Sn = " << flux.M_snow << std::endl;

        // diagnostics
        if (diagnostics_) {
          (*evap_rate)[0][c] -= area_fracs[1][c] * mb.Me;
          (*qE_sh)[0][c] += area_fracs[1][c] * eb.fQh;
          (*qE_lh)[0][c] += area_fracs[1][c] * eb.fQe;
          (*qE_lw_out)[0][c] += area_fracs[1][c] * eb.fQlwOut;
          (*qE_cond)[0][c] += area_fracs[1][c] * eb.fQc;

          (*qE_sm)[0][c] = area_fracs[1][c] * eb.fQm;
          (*melt_rate)[0][c] = area_fracs[1][c] * mb.Mm;
          (*snow_temp)[0][c] = snow.temp;
          (*albedo)[0][c] += area_fracs[1][c] * surf.albedo;
        }
      }
    }
  }

  // debugging
  if (diagnostics_ && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Surface Balance calculation:" << std::endl;
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("area fractions");
    vecs.push_back(S->GetFieldData(area_frac_key_).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.clear();
    vecs.clear();
    vnames.push_back("air_temp");
    vecs.push_back(S->GetFieldData(met_air_temp_key_).ptr());
    vnames.push_back("rel_hum");
    vecs.push_back(S->GetFieldData(met_rel_hum_key_).ptr());
    vnames.push_back("precip_rain");
    vecs.push_back(S->GetFieldData(met_prain_key_).ptr());
    vnames.push_back("precip_snow");
    vecs.push_back(S->GetFieldData(met_psnow_key_).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("p_ground");
    vecs.push_back(S->GetFieldData(surf_pres_key_).ptr());
    vnames.push_back("ponded_depth");
    vecs.push_back(S->GetFieldData(ponded_depth_key_).ptr());
    vnames.push_back("unfrozen_fraction");
    vecs.push_back(S->GetFieldData(unfrozen_fraction_key_).ptr());
    vnames.push_back("snow_depth");
    vecs.push_back(S->GetFieldData(snow_depth_key_).ptr());
    vnames.push_back("snow_death");
    vecs.push_back(S->GetFieldData(snow_death_rate_key_).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("T_ground");
    vecs.push_back(S->GetFieldData(surf_temp_key_).ptr());
    vnames.push_back("snow_temp");
    vecs.push_back(S->GetFieldData(snow_temp_key_).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("inc shortwave radiation");
    vecs.push_back(S->GetFieldData(met_sw_key_).ptr());
    vnames.push_back("inc longwave radiation");
    vecs.push_back(S->GetFieldData(met_lw_key_).ptr());
    vnames.push_back("inc latent heat");
    vecs.push_back(S->GetFieldData(qE_lh_key_).ptr());
    vnames.push_back("inc sensible heat");
    vecs.push_back(S->GetFieldData(qE_sh_key_).ptr());
    vnames.push_back("out longwave radiation");
    vecs.push_back(S->GetFieldData(qE_lw_out_key_).ptr());
    vnames.push_back("out conducted energy");
    vecs.push_back(S->GetFieldData(qE_cond_key_).ptr());
    vnames.push_back("out melting energy");
    vecs.push_back(S->GetFieldData(qE_sm_key_).ptr());

    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("water_source");
    vecs.push_back(S->GetFieldData(water_source_key_).ptr());
    vnames.push_back("evap flux");
    vecs.push_back(S->GetFieldData(evap_key_).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();
    vnames.push_back("mass src");
    vnames.push_back("energy src");
    vnames.push_back("snow src");
    vnames.push_back("new snow");
    vecs.push_back(results[0]);
    vecs.push_back(results[1]);
    vecs.push_back(results[4]);
    vecs.push_back(results[5]);
    db_->WriteVectors(vnames, vecs, true);

    vnames.clear();
    vecs.clear();
    vnames.push_back("sub mass src");
    vnames.push_back("sub energy src");
    vecs.push_back(results[2]);
    vecs.push_back(results[3]);
    db_ss_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();
  }
}

void
SEBTwoComponentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) {}


void
SEBTwoComponentEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  if (!compatible_) {
    if (db_ == Teuchos::null)
      db_ = Teuchos::rcp(new Debugger(S->GetMesh(domain_), my_keys_[0], plist_));
    if (db_ss_ == Teuchos::null) {
      Teuchos::ParameterList plist(plist_);
      plist.remove("debug cells", false);
      plist.remove("debug faces", false);
      if (plist.isParameter("subsurface debug cells")) {
        plist.set("debug cells", plist.get<Teuchos::Array<int>>("subsurface debug cells"));
      } else {
      }
      if (plist.isParameter("subsurface debug faces"))
        plist.set("debug faces", plist.get<Teuchos::Array<int>>("subsurface debug faces"));
      db_ss_ = Teuchos::rcp(new Debugger(S->GetMesh(domain_ss_), my_keys_[0], plist));
    }

    if (land_cover_.size() == 0)
      land_cover_ = getLandCover(S->ICList().sublist("land cover types"),
              {"roughness_snow", "roughness_ground",
               "water_transition_depth", "snow_transition_depth",
               "dessicated_zone_thickness"});

    CompositeVectorSpace domain_fac;
    domain_fac.SetMesh(S->GetMesh(domain_))
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    CompositeVectorSpace domain_fac_owned;
    domain_fac_owned.SetMesh(S->GetMesh(domain_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    CompositeVectorSpace domain_fac_owned_ss;
    domain_fac_owned_ss.SetMesh(S->GetMesh(domain_ss_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    CompositeVectorSpace domain_fac_owned_snow;
    domain_fac_owned_snow.SetMesh(S->GetMesh(domain_snow_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    CompositeVectorSpace domain_fac_2;
    domain_fac_2.SetMesh(S->GetMesh(domain_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 2);

    CompositeVectorSpace domain_fac_ss;
    domain_fac_ss.SetMesh(S->GetMesh(domain_ss_))
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    CompositeVectorSpace domain_fac_snow;
    domain_fac_snow.SetMesh(S->GetMesh(domain_snow_))
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    for (auto my_key : my_keys_) {
      auto my_fac = S->RequireField(my_key, my_key);
      if (Keys::getDomain(my_key) == domain_snow_) {
        my_fac->Update(domain_fac_owned_snow);
      } else if (Keys::getDomain(my_key) == domain_) {
        my_fac->Update(domain_fac_owned);
      } else if (Keys::getDomain(my_key) == domain_ss_) {
        my_fac->Update(domain_fac_owned_ss);
      } else {
        Errors::Message message("SEBTwoComponentEvaluator: Key requested with unrecognizable domain name.");
        Exceptions::amanzi_throw(message);
      }

      // Check plist for vis or checkpointing control.
      bool io_my_key = plist_.get<bool>("visualize", true);
      S->GetField(my_key, my_key)->set_io_vis(io_my_key);
      bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
      S->GetField(my_key, my_key)->set_io_checkpoint(checkpoint_my_key);
    }

    if (diagnostics_) {
      S->RequireField(albedo_key_, albedo_key_)->Update(domain_fac_owned);
      S->RequireField(melt_key_, melt_key_)->Update(domain_fac_owned);
      S->RequireField(evap_key_, evap_key_)->Update(domain_fac_owned);
      S->RequireField(snow_temp_key_, snow_temp_key_)->Update(domain_fac_owned);
      S->RequireField(qE_sh_key_, qE_sh_key_)->Update(domain_fac_owned);
      S->RequireField(qE_lh_key_, qE_lh_key_)->Update(domain_fac_owned);
      S->RequireField(qE_sm_key_, qE_sm_key_)->Update(domain_fac_owned);
      S->RequireField(qE_lw_out_key_, qE_lw_out_key_)->Update(domain_fac_owned);
      S->RequireField(qE_cond_key_, qE_cond_key_)->Update(domain_fac_owned);
    }

    for (auto dep_key : dependencies_) {
      auto fac = S->RequireField(dep_key);
      std::cout << "Keys::getDomain(dep_key) = " << Keys::getDomain(dep_key) << std::endl;
      std::cout << "domain_ss_ = " << domain_ss_ << std::endl;
      std::cout << "domain_snow_ = " << domain_snow_ << std::endl;
      if (Keys::getDomain(dep_key) == domain_ss_) {
        fac->Update(domain_fac_ss);
      } else if (dep_key == sg_albedo_key_ ||
                 dep_key == sg_emissivity_key_ ||
                 dep_key == area_frac_key_) {
        fac->Update(domain_fac_2);
      } else if (Keys::getDomain(dep_key) == domain_snow_) {
        fac->Update(domain_fac_snow);
      } else {
        std::cout << "before Update" << std::endl;
        fac->Update(domain_fac);
        std::cout << "after Update" << std::endl;
      }

      std::cout << "dep_key = " << dep_key << std::endl;
      S->RequireFieldEvaluator(dep_key)->EnsureCompatibility(S);
      std::cout << "compatible" << std::endl;
    }

    // additionally MANUALLY require the area frac, because it is not in the
    // list of dependencies :ISSUE:#8
    //S->RequireField(area_frac_key_)->Update(domain_fac_2);

    compatible_ = true;
  }
}

void
SEBTwoComponentEvaluator::UpdateFieldDerivative_(const Teuchos::Ptr<State>& S, Key wrt_key)
{
  Errors::Message message("SEBTwoComponentEvaluator: cannot differentiate with respect to anything.");
  Exceptions::amanzi_throw(message);
}

}  // namespace Relations
}  // namespace SurfaceBalance
}  // namespace Amanzi
