/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   Amanzi Chemistry

   License: see COPYRIGHT
   Author: mostly from old State

   Interface layer between Chemistry_PK and State, this is a harness for
   accessing the new state-dev from the Chemistry PK.  It also manages the
   mineralogy and names.

   ------------------------------------------------------------------------- */

#include "Chemistry_State.hh"

#include "beaker.hh"

namespace Amanzi {
namespace AmanziChemistry {

Chemistry_State::Chemistry_State(Teuchos::ParameterList& plist,
        const Teuchos::RCP<State>& S) :
    PK_State(std::string("state"), S),
    plist_(plist),
    number_of_aqueous_components_(-1),
    number_of_minerals_(-1),
    number_of_ion_exchange_sites_(-1),
    number_of_sorption_sites_(-1),
    using_sorption_(false),
    using_sorption_isotherms_(false) {

  SetupSoluteNames_();
  SetupMineralNames_();
  SetupSorptionSiteNames_();
  ParseMeshBlocks_();
  RequireData_();
}


Chemistry_State::Chemistry_State(const Teuchos::RCP<State>& S,
        int number_of_aqueous_components,
        int number_of_minerals,
        int number_of_ion_exchange_sites,
        int number_of_sorption_sites,
        bool using_sorption,
        bool using_sorption_isotherms) :
    PK_State(std::string("state"), S),
    number_of_aqueous_components_(number_of_aqueous_components),
    number_of_minerals_(number_of_minerals),
    number_of_ion_exchange_sites_(number_of_ion_exchange_sites),
    number_of_sorption_sites_(number_of_sorption_sites),
    using_sorption_(using_sorption),
    using_sorption_isotherms_(using_sorption_isotherms) {
  RequireData_();
}


void Chemistry_State::SetupSoluteNames_() {
  // get the number of component concentrations from the
  // parameter list
  if (plist_.isParameter("Number of component concentrations")) {
    number_of_aqueous_components_ =
        plist_.get<int>("Number of component concentrations");
  } else {
    // if the parameter list does not contain this key, then assume we
    // are being called from a state constructor w/o a valid parameter
    // list (vis/restart) and the number_of_components variable was
    // already set to a valid (non-zero?) value....
  }

  if (number_of_aqueous_components_ > 0) {
    // read the component names if they are spelled out
    Teuchos::Array<std::string> comp_names;
    if (plist_.isParameter("Component Solutes")) {
      comp_names = plist_.get<Teuchos::Array<std::string> >("Component Solutes");
    }
  }
}  // end SetupSoluteNames()

void Chemistry_State::SetupMineralNames_() {
  // do we need to worry about minerals?
  mineral_names_.clear();
  Teuchos::Array<std::string> data;
  if (plist_.isParameter("Minerals")) {
    data = plist_.get<Teuchos::Array<std::string> >("Minerals");
  }

  // the mineral_names_ list should be the order expected by the chemistry....
  mineral_names_.clear();
  mineral_name_id_map_.clear();
  for (int m = 0; m < data.size(); ++m) {
    mineral_name_id_map_[data.at(m)] = m;
    mineral_names_.push_back(data.at(m));
  }

  if (mineral_names_.size() > 0) {
    // we read some mineral names, so override any value that may have
    // been set in the constructor
    number_of_minerals_ = mineral_names_.size();
  }
}  // end SetupMineralNames()


void Chemistry_State::SetupSorptionSiteNames_() {
  // could almost generalize the SetupMineralNames and
  // SetupSorptionSiteNames into a single function w/ different
  // parameters, but using_sorption needs to be set...

  // do we need to worry about sorption sites?
  sorption_site_names_.clear();
  Teuchos::Array<std::string> data;
  if (plist_.isParameter("Sorption Sites")) {
    data = plist_.get<Teuchos::Array<std::string> >("Sorption Sites");
  }

  // the sorption_site_names_ list should be the order expected by the chemistry...
  sorption_site_names_.clear();
  sorption_site_name_id_map_.clear();
  for (int s = 0; s < data.size(); ++s) {
    sorption_site_name_id_map_[data.at(s)] = s;
    sorption_site_names_.push_back(data.at(s));
  }

  if (sorption_site_names_.size() > 0) {
    // we read some sorption site names, so override any value that
    // may have been set in the constructor and set the sorption flag
    // so we allocate the correct amount of memory
    number_of_sorption_sites_ = sorption_site_names_.size();
    using_sorption_ = true;
  } else if (number_of_sorption_sites() > 0 && sorption_site_names_.size() == 0) {
    // assume we are called from the constructor w/o a valid parameter
    // list and the sorption names will be set later....?
  }
}  // end SetupSorptionSiteNames()


void Chemistry_State::ParseMeshBlocks_() {
  const std::string block_key("Mesh block ");

  // loop through the state parameter list looking for mesh blocks
  for (Teuchos::ParameterList::ConstIterator item = plist_.begin();
       item != plist_.end(); ++item) {

    std::string item_name = plist_.name(item);
    size_t found_block = item_name.find(block_key);
    if (found_block != std::string::npos) {
      std::string block_name = item_name.substr(found_block + block_key.length(),
              item_name.length());
      Teuchos::ParameterList mesh_block_data = plist_.sublist(plist_.name(item));
      // block_name should match mesh_block_data.region...
      std::string region_name = mesh_block_data.get<std::string>("Region");

      //
      // check for chemistry related data in the block:
      //

      if (mesh_block_data.isSublist("Mineralogy")) {
        VerifyMineralogy_(region_name, mesh_block_data.sublist("Mineralogy"));
      }

      if (mesh_block_data.isSublist("Sorption Isotherms")) {
        VerifySorptionIsotherms_(region_name,
                mesh_block_data.sublist("Sorption Isotherms"));
      }

      if (mesh_block_data.isSublist("Surface Complexation Sites")) {
        VerifySorptionSites_(region_name,
                             mesh_block_data.sublist("Surface Complexation Sites"));
      }

      if (mesh_block_data.isParameter("Cation Exchange Capacity")) {
        // limit to one ion exchange site for now....
        using_sorption_ = true;
        number_of_ion_exchange_sites_ = 1;
      }
    }  // end if(mesh_block)
  }  // end for(plist_)
}


void Chemistry_State::VerifyMineralogy_(const std::string& region_name,
                             const Teuchos::ParameterList& minerals_list) {
  // loop through each mineral, verify that the mineral name is known
  for (Teuchos::ParameterList::ConstIterator mineral_iter = minerals_list.begin();
       mineral_iter != minerals_list.end(); ++mineral_iter) {
    std::string mineral_name = minerals_list.name(mineral_iter);
    if (!mineral_name_id_map_.count(mineral_name)) {
      std::stringstream message;
      message << "Error: Chemistry_State::VerifyMineralogy(): " << mineral_name
              << " was specified in the mineralogy for region "
              << region_name << " but was not listed in the minerals phase list.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));
    }

    // all minerals will have a volume fraction and specific surface
    // area, but sane defaults can be provided, so we don't bother
    // with them here.

  }  // end for(minerals)
}  // end VerifyMineralogy()

void Chemistry_State::VerifySorptionIsotherms_(const std::string& region_name,
        const Teuchos::ParameterList& isotherms_list) {
  // verify that every species listed is in the component names list.
  // verify that every listed species has a Kd value (no sane default)
  // langmuir and freundlich values are optional (sane defaults)
  using_sorption_ = true;
  using_sorption_isotherms_ = true;

  // loop through each species in the isotherm list
  for (Teuchos::ParameterList::ConstIterator species_iter = isotherms_list.begin();
       species_iter != isotherms_list.end(); ++species_iter) {
    std::string species_name = isotherms_list.name(species_iter);

    // verify that the name is a known species
    if (!comp_name_id_map_.count(species_name)) {
      std::stringstream message;
      message << "Error: Chemistry_State::VerifySorptionIsotherms(): region: "
              << region_name << " contains isotherm data for solute \'"
              << species_name
              << "\' but it is not specified in the component solutes list.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));
    }


    // check that this item is a sublist:
    Teuchos::ParameterList species_data;
    if (!isotherms_list.isSublist(species_name)) {
      std::stringstream message;
      message << "Error: Chemistry_State::VerifySorptionIsotherms(): region: "
              << region_name << " ; species : " << species_name
              << " ; must be a named \'ParameterList\' of isotherm data.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));
    } else {
      species_data = isotherms_list.sublist(species_name);
    }

    // verify that the required parameters are present
    if (!species_data.isParameter("Kd")) {
      std::stringstream message;
      message << "Error: Chemistry_State::VerifySorptionIsotherms(): region: "
              << region_name << " ; species name: " << species_name
              << " ; each isotherm must have a 'Kd' parameter.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));
    }
    // langmuir and freundlich parameters are optional, we'll assign
    // sane defaults.
  }  // end for(species)
}  // end VerifySorptionIsotherms()

void Chemistry_State::VerifySorptionSites_(const std::string& region_name,
        const Teuchos::ParameterList& sorption_site_list) {
  using_sorption_ = true;
  // loop through each sorption site, verify that the site name is known
  for (Teuchos::ParameterList::ConstIterator site_iter = sorption_site_list.begin();
       site_iter != sorption_site_list.end(); ++site_iter) {
    std::string site_name = sorption_site_list.name(site_iter);
    if (!sorption_site_name_id_map_.count(site_name)) {
      std::stringstream message;
      message << "Error: Chemistry_State::VerifySorptionSites(): " << site_name
              << " was specified in the 'Surface Complexation Sites' list for region "
              << region_name << " but was not listed in the sorption sites phase list.\n";
      Exceptions::amanzi_throw(Errors::Message(message.str()));
    }

    // all sorption sites will have a site density
    // but we can default to zero, so don't do any further checking

  }  // end for(sorption_sites)
}  // end VerifySorptionSites()



void Chemistry_State::RequireData_() {
  // Apparently Chemistry and Transport/Flow do not agree upon what water density is?
  water_density_ = Teuchos::rcp(new Epetra_Vector(S_->GetMesh()->cell_map(false)));
  water_density_initialized_ = false;

  // Require data from flow
  S_->RequireField("porosity", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireField("water_saturation", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireScalar("fluid_density", name_);
  S_->RequireFieldEvaluator("cell_volume");

  // Require my data
  if (number_of_aqueous_components_ > 0) {
    // TCC
    if (compnames_.size() != number_of_aqueous_components_) {
      compnames_.clear();
      for (int i=0; i!=number_of_aqueous_components_; ++i) {
        std::stringstream ss;
        ss << "Component " << i;
        compnames_.push_back(ss.str());
      }
    }

    // set the names for vis
    std::vector<std::vector<std::string> > conc_names_cv(1);
    for (std::vector<std::string>::const_iterator compname = compnames_.begin();
         compname != compnames_.end(); ++compname) {
      conc_names_cv[0].push_back(*compname + std::string(" conc"));
    }
    S_->RequireField("total_component_concentration", name_, conc_names_cv)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);

    // now create the map
    for (int i=0; i!=number_of_aqueous_components_; ++i) {
      comp_name_id_map_[compnames_[i]] = i;
    }


    // CreateStoragePrimarySpecies()
    S_->RequireField("free_ion_species", name_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);
    S_->RequireField("primary_activity_coeff", name_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);

    // CreateStorageTotalSorbed()
    if (using_sorption_) {
      S_->RequireField("total_sorbed", name_)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);
    }

    // CreateStorageSorptionIsotherms()
    if (using_sorption_isotherms_) {
      S_->RequireField("isotherm_kd", name_)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);
      S_->RequireField("isotherm_freundlich_n", name_)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);
      S_->RequireField("isotherm_langmuir_b", name_)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);
    }
  }

  // CreateStorageMinerals()
  if (number_of_minerals_ > 0) {
    if (mineral_names_.size() > 0) {
      // set the names for vis
      ASSERT(mineral_names_.size() == number_of_minerals_);
      std::vector<std::vector<std::string> > vf_names_cv(1);
      std::vector<std::vector<std::string> > ssa_names_cv(1);
      for (std::vector<std::string>::const_iterator mineral_name=mineral_names_.begin();
           mineral_name!=mineral_names_.end(); ++mineral_name) {
        vf_names_cv[0].push_back(*mineral_name + std::string(" vol frac"));
        ssa_names_cv[0].push_back(*mineral_name + std::string(" spec surf area"));
      }
      S_->RequireField("mineral_volume_fractions", name_, vf_names_cv)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_minerals_);
      S_->RequireField("mineral_specific_surface_area", name_, ssa_names_cv)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_minerals_);
    } else {
      S_->RequireField("mineral_volume_fractions", name_)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_minerals_);
      S_->RequireField("mineral_specific_surface_area", name_)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_minerals_);
    }
  }

  // CreateStorageIonExchange()
  if (number_of_ion_exchange_sites_ > 0) {
    S_->RequireField("ion_exchange_sites", name_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_of_ion_exchange_sites_);
    S_->RequireField("ion_exchange_ref_cation_conc", name_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_of_ion_exchange_sites_);
  }

  // CreateStorageSurfaceComplexation()
  if (number_of_sorption_sites_ > 0) {
    if (sorption_site_names_.size() > 0) {
      // set the names for vis
      ASSERT(sorption_site_names_.size() == number_of_sorption_sites_);
      std::vector<std::vector<std::string> > ss_names_cv(1);
      std::vector<std::vector<std::string> > scfsc_names_cv(1);
      for (std::vector<std::string>::const_iterator sorption_site_name=
               sorption_site_names_.begin();
           sorption_site_name!=sorption_site_names_.end(); ++sorption_site_name) {
        ss_names_cv[0].push_back(*sorption_site_name + std::string(" sorption site"));
        scfsc_names_cv[0].push_back(*sorption_site_name +
                std::string(" surface complex free site conc"));
      }
      S_->RequireField("sorption_sites", name_, ss_names_cv)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_sorption_sites_);
      S_->RequireField("surface_complex_free_site_conc", name_, scfsc_names_cv)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_sorption_sites_);
    } else {
      S_->RequireField("sorption_sites", name_)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_sorption_sites_);
      S_->RequireField("surface_complex_free_site_conc", name_)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, number_of_sorption_sites_);
    }
  }
}


void Chemistry_State::Initialize() {
  // initialize
  if (S_->HasField("free_ion_species")) {
    free_ion_species()->PutScalar(1.0e-9);
    S_->GetField("free_ion_species",name_)->set_initialized();
  }

  if (S_->HasField("primary_activity_coeff")) {
    primary_activity_coeff()->PutScalar(1.0);
    S_->GetField("primary_activity_coeff",name_)->set_initialized();
  }

  if (S_->HasField("total_sorbed")) {
    total_sorbed()->PutScalar(1.0);
    S_->GetField("total_sorbed",name_)->set_initialized();
  }

  if (S_->HasField("isotherm_kd")) {
    isotherm_kd()->PutScalar(0.0);
    S_->GetField("isotherm_kd",name_)->set_initialized();
  }

  if (S_->HasField("isotherm_freundlich_n")) {
    isotherm_freundlich_n()->PutScalar(1.0);
    S_->GetField("isotherm_freundlich_n",name_)->set_initialized();
  }

  if (S_->HasField("isotherm_langmuir_b")) {
    isotherm_langmuir_b()->PutScalar(1.0);
    S_->GetField("isotherm_langmuir_b",name_)->set_initialized();
  }

  if (S_->HasField("mineral_volume_fractions")) {
    mineral_volume_fractions()->PutScalar(0.0);
    S_->GetField("mineral_volume_fractions",name_)->set_initialized();
  }

  if (S_->HasField("mineral_specific_surface_area")) {
    mineral_specific_surface_area()->PutScalar(1.0);
    S_->GetField("mineral_specific_surface_area",name_)->set_initialized();
  }

  if (S_->HasField("ion_exchange_sites")) {
    ion_exchange_sites()->PutScalar(1.);
    S_->GetField("ion_exchange_sites",name_)->set_initialized();
  }

  if (S_->HasField("ion_exchange_ref_cation_conc")) {
    ion_exchange_ref_cation_conc()->PutScalar(1.);
    S_->GetField("ion_exchange_ref_cation_conc",name_)->set_initialized();
  }

  if (S_->HasField("sorption_sites")) {
    sorption_sites()->PutScalar(1.);
    S_->GetField("sorption_sites",name_)->set_initialized();
  }

  if (S_->HasField("surface_complex_free_site_conc")) {
    surface_complex_free_site_conc()->PutScalar(1.);
    S_->GetField("surface_complex_free_site_conc",name_)->set_initialized();
  }
}


// This can only be done AFTER the chemistry is initialized and fully set up?
void Chemistry_State::AllocateAdditionalChemistryStorage(
    const Beaker::BeakerComponents& components) {
  int n_secondary_comps = components.secondary_activity_coeff.size();
  if (n_secondary_comps > 0) {
    // CreateStorageSecondaryActivityCoeff()
    S_->RequireField("secondary_activity_coeff", name_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, n_secondary_comps);
    S_->GetField("secondary_activity_coeff",name_)->CreateData();
    S_->GetFieldData("secondary_activity_coeff",name_)->PutScalar(1.0);
  }
}

} // namespace AmanziChemistry
} // namespace Amanzi
