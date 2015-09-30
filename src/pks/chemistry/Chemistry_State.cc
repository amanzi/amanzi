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

namespace Amanzi {
namespace AmanziChemistry {

Chemistry_State::Chemistry_State(Teuchos::ParameterList& plist,
                                 const std::vector<std::string>& component_names,
                                 const Teuchos::RCP<State>& S) :
    S_(S),
    ghosted_(true),
    name_("state"),
    plist_(plist),
    number_of_aqueous_components_(component_names.size()),
    number_of_minerals_(0),
    number_of_ion_exchange_sites_(0),
    number_of_sorption_sites_(0),
    using_sorption_(false),
    using_sorption_isotherms_(false),
    compnames_(component_names),
    num_aux_data_(-1) {

  mesh_ = S_->GetMesh();
  // SetupSoluteNames_();
  SetupMineralNames_();
  SetupSorptionSiteNames_();

  // in the old version, this was only in the Block sublist... may need work?
  if (plist_.isParameter("Cation Exchange Capacity")) {
    using_sorption_ = true;
    number_of_ion_exchange_sites_ = 1;
  }

  ParseMeshBlocks_();
  // RequireData_();
  // RequireAuxData_();
}

void Chemistry_State::Setup() {
  RequireData_();  
  RequireAuxData_();
}


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

#if 0
void Chemistry_State::SetupSoluteNames_() {
  // get the number of component concentrations from the
  // parameter list
  if (plist_.isParameter("number of component concentrations")) {
    number_of_aqueous_components_ =
        plist_.get<int>("number of component concentrations");
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
}  // end SetupSoluteNames_()
#endif

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

  // check if there is an initial condition for ion_exchange_sites
  if (plist_.sublist("initial conditions").isSublist("ion_exchange_sites")) {
    // there is currently only at most one site...
    using_sorption_ = true;
    number_of_ion_exchange_sites_ = 1;
  }

  if (plist_.sublist("initial conditions").isSublist("isotherm_kd")) {
    using_sorption_ = true;
    using_sorption_isotherms_ = true;
  }

  if (plist_.sublist("initial conditions").isSublist("sorption_sites")) {
    using_sorption_ = true;
    Teuchos::Array<std::string> ss_names_ = plist_.get<Teuchos::Array<std::string> >("Sorption Sites");
    number_of_sorption_sites_ = ss_names_.size();
  }


  // const std::string block_key("Mesh block ");
  // // loop through the state parameter list looking for mesh blocks
  // for (Teuchos::ParameterList::ConstIterator item = plist_.begin();
  //      item != plist_.end(); ++item) {

  //   std::string item_name = plist_.name(item);
  //   size_t found_block = item_name.find(block_key);
  //   if (found_block != std::string::npos) {
  //     std::string block_name = item_name.substr(found_block + block_key.length(),
  //             item_name.length());
  //     Teuchos::ParameterList mesh_block_data = plist_.sublist(plist_.name(item));
  //     // block_name should match mesh_block_data.region...
  //     std::string region_name = mesh_block_data.get<std::string>("Region");

  //     //
  //     // check for chemistry related data in the block:
  //     //

  //     if (mesh_block_data.isSublist("Mineralogy")) {
  //       VerifyMineralogy_(region_name, mesh_block_data.sublist("Mineralogy"));
  //     }

  //     if (mesh_block_data.isSublist("Sorption Isotherms")) {
  //       VerifySorptionIsotherms_(region_name,
  //               mesh_block_data.sublist("Sorption Isotherms"));
  //     }

  //     if (mesh_block_data.isSublist("Surface Complexation Sites")) {
  //       VerifySorptionSites_(region_name,
  //                            mesh_block_data.sublist("Surface Complexation Sites"));
  //     }

  //     if (mesh_block_data.isParameter("Cation Exchange Capacity")) {
  //       // limit to one ion exchange site for now....
  //       using_sorption_ = true;
  //       number_of_ion_exchange_sites_ = 1;
  //     }
  //   }  // end if(mesh_block)
  // }  // end for(plist_)
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
  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("saturation_liquid")) {
    S_->RequireField("saturation_liquid", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  
  S_->RequireScalar("fluid_density", name_);
  if (!S_->HasFieldEvaluator("cell_volume")){
    S_->RequireFieldEvaluator("cell_volume");
  }

  // Require my data
  if (number_of_aqueous_components_ > 0) {

    // Make dummy names if our component names haven't already been set.
    if (compnames_.empty())
    {
      for (int i = 0; i < number_of_aqueous_components_; ++i)
      {
        std::stringstream ss;
        ss << "Component " << i << std::ends;
        compnames_.push_back(ss.str());
      }
    }

    // TCC
    if (!S_->HasField("total_component_concentration")) {    
      // set the names for vis
      std::vector<std::vector<std::string> > conc_names_cv(1);
      for (std::vector<std::string>::const_iterator compname = compnames_.begin();
           compname != compnames_.end(); ++compname) {
        conc_names_cv[0].push_back(*compname + std::string(" conc"));
      }
      S_->RequireField("total_component_concentration", name_, conc_names_cv)
        ->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);
    }

    // now create the map
    for (int i=0; i!=number_of_aqueous_components_; ++i) {
      comp_name_id_map_[compnames_[i]] = i;
    }


    // CreateStoragePrimarySpecies()
    {
      std::vector<std::vector<std::string> > species_names_cv(1);
      for (std::vector<std::string>::const_iterator compname = compnames_.begin();
           compname != compnames_.end(); ++compname) {
        species_names_cv[0].push_back(*compname);
      }
      S_->RequireField("free_ion_species", name_, species_names_cv)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);
      S_->RequireField("primary_activity_coeff", name_, species_names_cv)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_of_aqueous_components_);
    }

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


void Chemistry_State::SetAuxDataNames(const std::vector<std::string>& aux_data_names) {
  for (size_t i = 0; i < aux_data_names.size(); ++i) {
    std::vector<std::vector<std::string> > subname(1);
    subname[0].push_back("0");
    if (!S_->HasField(aux_data_names[i])) {
      Teuchos::RCP<CompositeVectorSpace> fac = S_->RequireField(aux_data_names[i], name_, subname);
      fac->SetMesh(mesh_);
      fac->SetGhosted(false);
      fac->SetComponent("cell", AmanziMesh::CELL, 1);
      Teuchos::RCP<CompositeVector> sac = Teuchos::rcp(new CompositeVector(*fac));

      // Zero the field.
      S_->GetField(aux_data_names[i], name_)->SetData(sac);
      S_->GetField(aux_data_names[i], name_)->CreateData();
      S_->GetFieldData(aux_data_names[i], name_)->PutScalar(0.0);
      S_->GetField(aux_data_names[i], name_)->set_initialized();
    }
  }
}

void Chemistry_State::RequireAuxData_() {
  if (plist_.isParameter("auxiliary data"))  {
    Teuchos::Array<std::string> names = plist_.get<Teuchos::Array<std::string> >("auxiliary data");  
    
    for (Teuchos::Array<std::string>::const_iterator name = names.begin(); name != names.end(); ++name) {

      // Insert the field into the state.
      std::vector<std::vector<std::string> > subname(1);
      subname[0].push_back("0");
      S_->RequireField(*name, name_, subname)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, 1);
    }
  }
}

void Chemistry_State::InitializeField_(Teuchos::ParameterList& ic_plist,
    std::string fieldname, bool sane_default, double default_val) {
  // Initialize mineral volume fractions
  // -- first initialize to a default: this should be a valid default if the
  // parameter is optional, and non-valid if it is not.

  if (S_->HasField(fieldname)) {
    if (!S_->GetField(fieldname)->initialized())
      S_->GetFieldData(fieldname, name_)->PutScalar(default_val);
  }

  // -- initialize from the ParameterList
  if (ic_plist.isSublist(fieldname)) {
    if (!S_->GetField(fieldname)->initialized())
      S_->GetField(fieldname, name_)->Initialize(ic_plist.sublist(fieldname));
  } else if (sane_default) {
    // -- sane default provided, functional initialization not necessary
    S_->GetField(fieldname, name_)->set_initialized();
  }
}

void Chemistry_State::Initialize() {
  // Most things are initialized through State, but State can only manage that
  // if they are always initialized.  If sane defaults are available, or they
  // can be derived from other initialized quantities, they are initialized
  // here, where we can manage that logic.

  // initialize list
  Teuchos::ParameterList ic_plist = plist_.sublist("initial conditions");

  // Aqueous species:
  if (number_of_aqueous_components_ > 0) {
    if (!S_->GetField("total_component_concentration",name_)->initialized()) {
      InitializeField_(ic_plist, "total_component_concentration", false, -1.0);
    }
    InitializeField_(ic_plist, "free_ion_species", false, 0.0);
    InitializeField_(ic_plist, "primary_activity_coeff", false, 1.0);

    // Sorption sites: all will have a site density, but we can default to zero
    if (using_sorption_) {
      InitializeField_(ic_plist, "total_sorbed", false, 0.0);
    }

    // Sorption isotherms: Kd required, Langmuir and Freundlich optional
    if (using_sorption_isotherms_) {
      InitializeField_(ic_plist, "isotherm_kd", false, -1.0);
      InitializeField_(ic_plist, "isotherm_freundlich_n", false, 1.0);
      InitializeField_(ic_plist, "isotherm_langmuir_b", false, 1.0);
    }
  }

  // Minerals: vol frac and surface areas
  if (number_of_minerals_ > 0) {
    InitializeField_(ic_plist, "mineral_volume_fractions", false, 0.0);
    InitializeField_(ic_plist, "mineral_specific_surface_area", false, 1.0);
  }

  // Ion exchange sites: default to 1
  if (number_of_ion_exchange_sites_ > 0) {
    InitializeField_(ic_plist, "ion_exchange_sites", false, 1.0);
    InitializeField_(ic_plist, "ion_exchange_ref_cation_conc", false, 1.0);
  }

  if (number_of_sorption_sites_ > 0) {
    InitializeField_(ic_plist, "sorption_sites", false, 1.0);
    InitializeField_(ic_plist, "surface_complex_free_site_conc", false, 1.0);
  }

  // initialize auxiliary fields
  if (plist_.isParameter("auxiliary data"))  {
    Teuchos::Array<std::string> names = plist_.get<Teuchos::Array<std::string> >("auxiliary data");  
    
    for (Teuchos::Array<std::string>::const_iterator name = names.begin(); name != names.end(); ++name) {
      S_->GetFieldData(*name, name_)->PutScalar(0.0);
      S_->GetField(*name, name_)->set_initialized();
    }  
  }
}

// This can only be done AFTER the chemistry is initialized and fully set up?
void Chemistry_State::AllocateAdditionalChemistryStorage(
    const Beaker::BeakerComponents& components) {
  int n_secondary_comps = components.secondary_activity_coeff.size();
  if (n_secondary_comps > 0) {
    // CreateStorageSecondaryActivityCoeff()
    Teuchos::RCP<CompositeVectorSpace> fac =
        S_->RequireField("secondary_activity_coeff", name_);
    fac->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, n_secondary_comps);
    Teuchos::RCP<CompositeVector> sac = Teuchos::rcp(new CompositeVector(*fac));
    S_->GetField("secondary_activity_coeff",name_)->SetData(sac);
    S_->GetField("secondary_activity_coeff",name_)->CreateData();
    S_->GetFieldData("secondary_activity_coeff",name_)->PutScalar(1.0);
    S_->GetField("secondary_activity_coeff",name_)->set_initialized();
  }
}

// This can only be done AFTER the chemistry is initialized and fully set up?
// NOTE: This is the version of the above method that interacts with Alquimia.
void Chemistry_State::AllocateAdditionalChemistryStorage(int num_aqueous_components) {
  if (num_aqueous_components > 0) {
    // CreateStorageSecondaryActivityCoeff()
    Teuchos::RCP<CompositeVectorSpace> fac =
        S_->RequireField("secondary_activity_coeff", name_);
    fac->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, num_aqueous_components);
    Teuchos::RCP<CompositeVector> sac = Teuchos::rcp(new CompositeVector(*fac));    
    S_->GetField("secondary_activity_coeff",name_)->SetData(sac);
    S_->GetField("secondary_activity_coeff",name_)->CreateData();
    S_->GetFieldData("secondary_activity_coeff",name_)->PutScalar(1.0);
    S_->GetField("secondary_activity_coeff",name_)->set_initialized();
  }
}

#ifdef ALQUIMIA_ENABLED

void Chemistry_State::CopyToAlquimia(const int cell_id,
                                     AlquimiaMaterialProperties& mat_props,
                                     AlquimiaState& state,
                                     AlquimiaAuxiliaryData& aux_data)
{
  CopyToAlquimia(cell_id, total_component_concentration(), mat_props, state, aux_data);
}

void Chemistry_State::CopyToAlquimia(const int cell_id,
                                     Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                                     AlquimiaMaterialProperties& mat_props,
                                     AlquimiaState& state,
                                     AlquimiaAuxiliaryData& aux_data)
{
  state.water_density = (*this->water_density())[cell_id];
  state.porosity = (*this->porosity())[cell_id];

  for (unsigned int c = 0; c < number_of_aqueous_components(); c++) 
  {
    double* cell_components = (*aqueous_components)[c];
    double total_mobile = cell_components[cell_id];
    state.total_mobile.data[c] = total_mobile;
    if (using_sorption()) 
    {
      double* cell_total_sorbed = (*this->total_sorbed())[c];
      double immobile = cell_total_sorbed[cell_id];
      state.total_immobile.data[c] = immobile;
    }  // end if(using_sorption)
  }

  // minerals
  assert(state.mineral_volume_fraction.size == number_of_minerals());
  assert(state.mineral_specific_surface_area.size == number_of_minerals());
  for (unsigned int m = 0; m < number_of_minerals(); m++) 
  {
    double* cell_minerals = (*this->mineral_volume_fractions())[m];
    state.mineral_volume_fraction.data[m] = cell_minerals[cell_id];
    if (this->mineral_specific_surface_area() != Teuchos::null) 
    {
      double* cells_ssa = (*this->mineral_specific_surface_area())[m];
      state.mineral_specific_surface_area.data[m] = cells_ssa[cell_id];
    }
  }

  // ion exchange
  assert(state.cation_exchange_capacity.size == number_of_ion_exchange_sites());
  if (number_of_ion_exchange_sites() > 0) 
  {
    for (unsigned int i = 0; i < number_of_ion_exchange_sites(); i++) 
    {
      double* cell_ion_exchange_sites = (*this->ion_exchange_sites())[i];
      state.cation_exchange_capacity.data[i] = cell_ion_exchange_sites[cell_id];
    }
  }
  
  // surface complexation
  if (number_of_sorption_sites() > 0) 
  {
    assert(number_of_sorption_sites() == state.surface_site_density.size);
    for (int s = 0; s < number_of_sorption_sites(); ++s) 
    {
      // FIXME: Need site density names, too?
      double* cell_sorption_sites = (*this->sorption_sites())[s];
      state.surface_site_density.data[s] = cell_sorption_sites[cell_id];
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }

  // Auxiliary data -- block copy.
  if (S_->HasField("alquimia_aux_data"))
    aux_data_ = S_->GetField("alquimia_aux_data", name_)->GetFieldData()->ViewComponent("cell");
  if (num_aux_data_ != -1) 
  {
    int num_aux_ints = aux_data.aux_ints.size;
    int num_aux_doubles = aux_data.aux_doubles.size;
    for (int i = 0; i < num_aux_ints; i++) 
    {
      double* cell_aux_ints = (*aux_data_)[i];
      aux_data.aux_ints.data[i] = (int)cell_aux_ints[cell_id];
    }
    for (int i = 0; i < num_aux_doubles; i++) 
    {
      double* cell_aux_doubles = (*aux_data_)[i + num_aux_ints];
      aux_data.aux_doubles.data[i] = cell_aux_doubles[cell_id];
    }
  }

  mat_props.volume = (*this->volume())[cell_id];
  mat_props.saturation = (*this->water_saturation())[cell_id];

  // sorption isotherms
  if (using_sorption_isotherms()) 
  {
    for (unsigned int i = 0; i < number_of_aqueous_components(); ++i) 
    {
      double* cell_data = (*this->isotherm_kd())[i];
      mat_props.isotherm_kd.data[i] = cell_data[cell_id];
      
      cell_data = (*this->isotherm_freundlich_n())[i];
      mat_props.freundlich_n.data[i] = cell_data[cell_id];
      
      cell_data = (*this->isotherm_langmuir_b())[i];
      mat_props.langmuir_b.data[i] = cell_data[cell_id];
    }
  }
}

void Chemistry_State::CopyFromAlquimia(const int cell_id,
                                       const AlquimiaMaterialProperties& mat_props,
                                       const AlquimiaState& state,
                                       const AlquimiaAuxiliaryData& aux_data,
                                       const AlquimiaAuxiliaryOutputData& aux_output,
                                       Teuchos::RCP<const Epetra_MultiVector> aqueous_components)
{
  // If the chemistry has modified the porosity and/or density, it needs to 
  // be updated here.
  //(this->water_density())[cell_id] = state.water_density;
  //(this->porosity())[cell_id] = state.porosity;
  for (unsigned int c = 0; c < number_of_aqueous_components(); c++) 
  {
    double mobile = state.total_mobile.data[c];
    double* cell_components = (*aqueous_components)[c];
    cell_components[cell_id] = mobile;

    if (using_sorption())
    {
      double immobile = state.total_immobile.data[c];
      double* cell_total_sorbed = (*this->total_sorbed())[c];
      cell_total_sorbed[cell_id] = immobile;
    }
  }

  // Free ion species.
  for (unsigned int c = 0; c < number_of_aqueous_components(); c++) {
    double* cell_free_ion = (*this->free_ion_species())[c];
    cell_free_ion[cell_id] = aux_output.primary_free_ion_concentration.data[c];
  }

  // Mineral properties.
  for (unsigned int m = 0; m < number_of_minerals(); m++) 
  {
    double* cell_minerals = (*this->mineral_volume_fractions())[m];
    cell_minerals[cell_id] = state.mineral_volume_fraction.data[m];
    if (this->mineral_specific_surface_area() != Teuchos::null) 
    {
      cell_minerals = (*this->mineral_specific_surface_area())[m];
      cell_minerals[cell_id] = state.mineral_specific_surface_area.data[m];
    }
  }

  // ion exchange
  for (unsigned int i = 0; i < number_of_ion_exchange_sites(); i++) 
  {
    double* cell_ion_exchange_sites = (*this->ion_exchange_sites())[i];
    cell_ion_exchange_sites[cell_id] = state.cation_exchange_capacity.data[i];
  }

  // surface complexation
  if (number_of_sorption_sites() > 0)
  {
    for (unsigned int i = 0; i < number_of_sorption_sites(); i++) 
    {
      double* cell_sorption_sites = (*this->sorption_sites())[i];
      cell_sorption_sites[cell_id] = state.surface_site_density.data[i];
    }
  }

  // Auxiliary data -- block copy.
  int num_aux_ints = aux_data.aux_ints.size;
  int num_aux_doubles = aux_data.aux_doubles.size;
  if (num_aux_data_ == -1) 
  {
    // Set things up and register a vector in the State.
    assert(num_aux_ints >= 0);
    assert(num_aux_doubles >= 0);
    num_aux_data_ = num_aux_ints + num_aux_doubles;
    if (!S_->HasField("alquimia_aux_data"))
    {
      Teuchos::RCP<CompositeVectorSpace> fac = S_->RequireField("alquimia_aux_data", name_);
      fac->SetMesh(mesh_);
      fac->SetGhosted(false);
      fac->SetComponent("cell", AmanziMesh::CELL, num_aux_data_);
      Teuchos::RCP<CompositeVector> sac = Teuchos::rcp(new CompositeVector(*fac));

      // Zero the field.
      Teuchos::RCP<Field> F = S_->GetField("alquimia_aux_data", name_);
      F->SetData(sac);
      F->CreateData();
      F->GetFieldData()->PutScalar(0.0);
      F->set_initialized();
    }
    aux_data_ = S_->GetField("alquimia_aux_data", name_)->GetFieldData()->ViewComponent("cell");
  }
  else
  {
    assert(num_aux_data_ == num_aux_ints + num_aux_doubles);
  }
  for (int i = 0; i < num_aux_ints; i++) 
  {
    double* cell_aux_ints = (*aux_data_)[i];
    cell_aux_ints[cell_id] = (double)aux_data.aux_ints.data[i];
  }
  for (int i = 0; i < num_aux_doubles; i++) 
  {
    double* cell_aux_doubles = (*aux_data_)[i + num_aux_ints];
    cell_aux_doubles[cell_id] = aux_data.aux_doubles.data[i];
  }

  if (using_sorption_isotherms()) 
  {
    for (unsigned int i = 0; i < number_of_aqueous_components(); ++i) 
    {
      double* cell_data = (*this->isotherm_kd())[i];
      cell_data[cell_id] = mat_props.isotherm_kd.data[i];

      cell_data = (*this->isotherm_freundlich_n())[i];
      cell_data[cell_id] = mat_props.freundlich_n.data[i];

      cell_data = (*this->isotherm_langmuir_b())[i];
      cell_data[cell_id] = mat_props.langmuir_b.data[i];
    }
  }
}

void Chemistry_State::InitFromAlquimia(const int cell_id,
                                       const AlquimiaMaterialProperties& mat_props,
                                       const AlquimiaState& state,
                                       const AlquimiaAuxiliaryData& aux_data,
                                       const AlquimiaAuxiliaryOutputData& aux_output)
{
  // If the chemistry has modified the porosity and/or density, it needs to 
  // be updated here.
  //(this->water_density())[cell_id] = state.water_density;
  //(this->porosity())[cell_id] = state.porosity;

 
  for (unsigned int c = 0; c < number_of_aqueous_components(); c++) {

    if (!S_->GetField("total_component_concentration")->initialized()){
      double mobile = state.total_mobile.data[c];
      double* cell_components = (*total_component_concentration())[c];
      cell_components[cell_id] = mobile;
    }

    if (using_sorption()){
      if (S_->HasField("total_sorbed")){
        if (!S_->GetField("total_sorbed")->initialized()) {
          double immobile = state.total_immobile.data[c];
          double* cell_total_sorbed = (*this->total_sorbed())[c];
          cell_total_sorbed[cell_id] = immobile;
        }
      }
    }
  }

  
  
  // Free ion species.
  for (unsigned int c = 0; c < number_of_aqueous_components(); c++) {
    if (!S_->GetField("free_ion_species")->initialized()) {
      double* cell_free_ion = (*this->free_ion_species())[c];
      cell_free_ion[cell_id] = aux_output.primary_free_ion_concentration.data[c];
    }
  }

  if (S_->HasField("mineral_volume_fractions")){
    if (!S_->GetField("mineral_volume_fractions")->initialized()) {
      // Mineral properties.
      for (unsigned int m = 0; m < number_of_minerals(); m++){
        double* cell_minerals = (*this->mineral_volume_fractions())[m];
        cell_minerals[cell_id] = state.mineral_volume_fraction.data[m];
        if (this->mineral_specific_surface_area() != Teuchos::null) {
          if (!S_->GetField("mineral_specific_surface_area")->initialized()) {
            cell_minerals = (*this->mineral_specific_surface_area())[m];
            cell_minerals[cell_id] = state.mineral_specific_surface_area.data[m];
          }
        }
      }
    }
  }



  if (S_->HasField("ion_exchange_sites")){
    if (!S_->GetField("ion_exchange_sites")->initialized()) {
      // ion exchange
      for (unsigned int i = 0; i < number_of_ion_exchange_sites(); i++) {
        double* cell_ion_exchange_sites = (*this->ion_exchange_sites())[i];
        cell_ion_exchange_sites[cell_id] = state.cation_exchange_capacity.data[i];
      }
    }
  }

  // surface complexation
  if (number_of_sorption_sites() > 0){
    if (!S_->GetField("sorption_sites")->initialized()) {    
      for (unsigned int i = 0; i < number_of_sorption_sites(); i++){
        double* cell_sorption_sites = (*this->sorption_sites())[i];
        cell_sorption_sites[cell_id] = state.surface_site_density.data[i];
      }
    }
  }


  // Auxiliary data -- block copy.
  int num_aux_ints = aux_data.aux_ints.size;
  int num_aux_doubles = aux_data.aux_doubles.size;
  if (num_aux_data_ == -1){
    // Set things up and register a vector in the State.
    assert(num_aux_ints >= 0);
    assert(num_aux_doubles >= 0);
    num_aux_data_ = num_aux_ints + num_aux_doubles;
    if (!S_->HasField("alquimia_aux_data")){
     
      Teuchos::RCP<CompositeVectorSpace> fac = S_->RequireField("alquimia_aux_data", name_);
      fac->SetMesh(mesh_);
      fac->SetGhosted(false);
      fac->SetComponent("cell", AmanziMesh::CELL, num_aux_data_);
      Teuchos::RCP<CompositeVector> sac = Teuchos::rcp(new CompositeVector(*fac));

      // Zero the field.
      Teuchos::RCP<Field> F = S_->GetField("alquimia_aux_data", name_);
      F->SetData(sac);
      F->CreateData();
      F->GetFieldData()->PutScalar(0.0);
      F->set_initialized(false);
        
    }

    aux_data_ = S_->GetField("alquimia_aux_data", name_)->GetFieldData()->ViewComponent("cell");
  }
  else {
    assert(num_aux_data_ == num_aux_ints + num_aux_doubles);
  }


  if (!S_->GetField("alquimia_aux_data", name_)->initialized()){
    for (int i = 0; i < num_aux_ints; i++){
      double* cell_aux_ints = (*aux_data_)[i];
      cell_aux_ints[cell_id] = (double)aux_data.aux_ints.data[i];
    }

    for (int i = 0; i < num_aux_doubles; i++){
      double* cell_aux_doubles = (*aux_data_)[i + num_aux_ints];
      cell_aux_doubles[cell_id] = aux_data.aux_doubles.data[i];
    }
  }

  if (using_sorption_isotherms()) {
    for (unsigned int i = 0; i < number_of_aqueous_components(); ++i){
       double* cell_data;
       if (!S_->GetField("isotherm_kd")->initialized() ){
         cell_data = (*this->isotherm_kd())[i];
         cell_data[cell_id] = mat_props.isotherm_kd.data[i];
       }

       if (!S_->GetField("isotherm_freundlich_n")->initialized() ){
         cell_data = (*this->isotherm_freundlich_n())[i];
         cell_data[cell_id] = mat_props.freundlich_n.data[i];
       }

       if (!S_->GetField("isotherm_langmuir_b")->initialized() ){
         cell_data = (*this->isotherm_langmuir_b())[i];
         cell_data[cell_id] = mat_props.langmuir_b.data[i];
       }

    }
  }

}



#endif

void Chemistry_State::InitFromBeakerStructure(const int cell_id,
                                              Beaker::BeakerComponents beaker_components){

  if (!S_->GetField("total_component_concentration")->initialized()){   
    for (unsigned int c = 0; c < number_of_aqueous_components(); c++) {
      double* cell_components = (*total_component_concentration())[c];
      cell_components[cell_id] = beaker_components.total.at(c);
    }
  }

  if (!S_->GetField("free_ion_species")->initialized()) {
    for (unsigned int c = 0; c < number_of_aqueous_components(); c++) {
      double* cell_free_ion = (*free_ion_species())[c];
      cell_free_ion[cell_id] = beaker_components.free_ion.at(c);
    }
  }

  //
  // activity coefficients
  //
  int n_primary_comps = beaker_components.primary_activity_coeff.size();
  if (n_primary_comps > 0) {
    if (!S_->GetField("primary_activity_coeff")->initialized()) {
      for (unsigned int i = 0; i < beaker_components.primary_activity_coeff.size(); ++i) {
        double* cells = (*primary_activity_coeff())[i];
        cells[cell_id] = beaker_components.primary_activity_coeff.at(i);
      }
    }
  }
  int n_secondary_comps = beaker_components.secondary_activity_coeff.size();
  if (n_secondary_comps > 0) {
    if (!S_->GetField("secondary_activity_coeff")->initialized()) {
      for (unsigned int i = 0; i < beaker_components.secondary_activity_coeff.size(); ++i) {
        double* cells = (*secondary_activity_coeff())[i];
        cells[cell_id] =  beaker_components.secondary_activity_coeff.at(i);
      }
    }
  }
  //
  // minerals
  //
  if (!S_->GetField("mineral_volume_fractions")->initialized()) {
    for (unsigned int m = 0; m < number_of_minerals_; m++) {
      double* cell_minerals = (*mineral_volume_fractions())[m];
      cell_minerals[cell_id] = beaker_components.mineral_volume_fraction.at(m);
      if (mineral_specific_surface_area() != Teuchos::null) {
        cell_minerals = (*mineral_specific_surface_area())[m];
        cell_minerals[cell_id] = beaker_components.mineral_specific_surface_area.at(m);
      }
    }
  }

  //
  // sorption
  //
  if (using_sorption()) {
    if (!S_->GetField("total_sorbed")->initialized()) {
      for (unsigned int c = 0; c < number_of_aqueous_components(); c++) {
        double* cell_total_sorbed = (*total_sorbed())[c];
        cell_total_sorbed[cell_id] = beaker_components.total_sorbed.at(c);
      }
    }
  }

  //
  // surface complexation
  //
  

  for (unsigned int i = 0; i < number_of_sorption_sites(); i++) {
    if (!S_->GetField("sorption_sites")->initialized()) {
      double* cell_sorption_sites = (*sorption_sites())[i];
      cell_sorption_sites[cell_id] = beaker_components.surface_site_density.at(i);
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }


 
  for (unsigned int i = 0; i < beaker_components.surface_complex_free_site_conc.size(); ++i) {
    if (!S_->GetField("surface_complex_free_site_conc")->initialized()) { 
     double* cells = (*surface_complex_free_site_conc())[i];
      cells[cell_id] = beaker_components.surface_complex_free_site_conc.at(i);
    }
  }

  //
  // ion exchange
  //

  for (unsigned int i = 0; i < number_of_ion_exchange_sites(); i++) {
    if (!S_->GetField("ion_exchange_sites")->initialized()) {
      double* cell_ion_exchange_sites = (*ion_exchange_sites())[i];
      cell_ion_exchange_sites[cell_id] = beaker_components.ion_exchange_sites.at(i);
      // TODO(bandre): need to save ion exchange ref cation conc here!
    }
  }


  for (unsigned int i = 0; i < beaker_components.ion_exchange_ref_cation_conc.size(); ++i) {
    if (!S_->GetField("ion_exchange_ref_cation_conc")->initialized()) {
      double* cells = (*ion_exchange_ref_cation_conc())[i];
      cells[cell_id] = beaker_components.ion_exchange_ref_cation_conc.at(i);
    }
  }

  //
  // sorption isotherms
  //
  if (using_sorption_isotherms()) {
    for (unsigned int i = 0; i < number_of_aqueous_components(); ++i) {
      if (!S_->GetField("isotherm_kd")->initialized()) {
        double* cell_data = (*isotherm_kd())[i];
        cell_data[cell_id] = beaker_components.isotherm_kd.at(i);
      }

      if (!S_->GetField("isotherm_freundlich_n")->initialized()) {
        double* cell_data = (*isotherm_freundlich_n())[i];
        cell_data[cell_id] = beaker_components.isotherm_freundlich_n.at(i);
      }

      if (!S_->GetField("isotherm_langmuir_b")->initialized()) {
        double *cell_data = (*isotherm_langmuir_b())[i];
        cell_data[cell_id] = beaker_components.isotherm_langmuir_b.at(i);
      }
    }
  }

}

void Chemistry_State::SetAllFieldsInitialized(){
  if (using_sorption_isotherms()) {
    S_->GetField("isotherm_langmuir_b", name_)->set_initialized();
    S_->GetField("isotherm_kd", name_)->set_initialized();
    S_->GetField("isotherm_freundlich_n", name_)->set_initialized();
  }
    
  if (S_->HasField("alquimia_aux_data")){
    S_->GetField("alquimia_aux_data", name_)->set_initialized();
  }
  if (number_of_sorption_sites() > 0){
    S_->GetField("sorption_sites", name_)->set_initialized();
  }
  if (S_->HasField("ion_exchange_sites")){
    S_->GetField("ion_exchange_sites", name_)->set_initialized();
  }
  if (S_->HasField("mineral_volume_fractions")){
    S_->GetField("mineral_volume_fractions", name_)->set_initialized();
    if (this->mineral_specific_surface_area() != Teuchos::null) {
      S_->GetField("mineral_specific_surface_area", name_)->set_initialized();
    }
  }

  if (S_->HasField("mineral_volume_fractions")){
    S_->GetField("mineral_volume_fractions", name_)->set_initialized();
  }

  if (number_of_aqueous_components() > 0){
    S_->GetField("free_ion_species", name_)->set_initialized();
    S_->GetField("total_component_concentration", name_)->set_initialized();
    if (S_->HasField("primary_activity_coeff")){
       S_->GetField("primary_activity_coeff", name_)->set_initialized();
    }
    if (using_sorption()){
      S_->GetField("total_sorbed", name_)->set_initialized();
    }
  }

}


} // namespace AmanziChemistry
} // namespace Amanzi
