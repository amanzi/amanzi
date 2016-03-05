/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors:

  Trilinos based process kernel for chemistry. Geochemistry
  calculations live in the chemistry library. The PK stores the 
  instance of the chemistry object and drives the chemistry 
  calculations on a cell by cell basis. It handles the movement of
  data back and forth between the State and the chemistry library 
  data structures.
*/
 
#include <string>
#include <algorithm>

// TPLs
#include "boost/mpi.hpp"
#include "boost/algorithm/string.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Mesh.hh"
#include "beaker.hh"
#include "chemistry_exception.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "simple_thermo_database.hh"
#include "VerboseObject.hh"

// Chemistry
#include "Amanzi_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

// This should go away
extern VerboseObject* chem_out;

/* ******************************************************************
* Constructor
******************************************************************* */
Amanzi_PK::Amanzi_PK(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& glist,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln) :
    soln_(soln),
    max_time_step_(9.9e9),
    chem_(NULL),
    current_time_(0.0),
    saved_time_(0.0)
{
  S_ = S;
  mesh_ = S_->GetMesh();

  // extract pk name
  std::string pk_name = pk_tree.name();
  boost::iterator_range<std::string::iterator> res = boost::algorithm::find_last(pk_name, "->"); 
  if (res.end() - pk_name.end() != 0) boost::algorithm::erase_head(pk_name, res.end() - pk_name.begin());

  // create pointer to the chemistry parameter list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  cp_list_ = Teuchos::sublist(pk_list, "Chemistry", true);

  // collect high-level information about the problem
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(glist, "State", true);

  InitializeMinerals(cp_list_);
  InitializeSorptionSites(cp_list_, state_list);

  // grab the component names
  comp_names_.clear();
  Teuchos::RCP<Teuchos::ParameterList> cd_list = Teuchos::sublist(glist, "Cycle Driver", true);
  if (cd_list->isParameter("component names")) {
    comp_names_ = cd_list->get<Teuchos::Array<std::string> >("component names").toVector();
  } else{
    Errors::Message msg("Amanzi_PK: Cycle Driver has no input parameter component names.");
    Exceptions::amanzi_throw(msg);
  }
  number_aqueous_components_ = comp_names_.size();
  number_free_ion_ = number_aqueous_components_;
  number_total_sorbed_ = number_aqueous_components_;

  // verbosity object
  vo_ = Teuchos::rcp(new VerboseObject("Chem::Amanzi", *cp_list_)); 
  chem_out = &*vo_;
}


/* *******************************************************************
* Destructor
******************************************************************* */
Amanzi_PK::~Amanzi_PK() {
  delete chem_;
}


/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void Amanzi_PK::Setup(const Teuchos::Ptr<State>& S)
{
  // use common registration steps
  Chemistry_PK::Setup(S);

  // no additional steps to be done.
}


/* ******************************************************************
* Can this be done during Setup phase?
******************************************************************* */
void Amanzi_PK::AllocateAdditionalChemistryStorage_(
    const Beaker::BeakerComponents& components)
{
  int n_secondary_comps = components.secondary_activity_coeff.size();
  if (n_secondary_comps > 0) {
    Teuchos::RCP<CompositeVectorSpace> fac = S_->RequireField("secondary_activity_coeff", passwd_);
    fac->SetMesh(mesh_)->SetGhosted(false)
       ->SetComponent("cell", AmanziMesh::CELL, n_secondary_comps);

    Teuchos::RCP<CompositeVector> sac = Teuchos::rcp(new CompositeVector(*fac));
    S_->GetField("secondary_activity_coeff", passwd_)->SetData(sac);
    S_->GetField("secondary_activity_coeff", passwd_)->CreateData();
    S_->GetFieldData("secondary_activity_coeff", passwd_)->PutScalar(1.0);
    S_->GetField("secondary_activity_coeff", passwd_)->set_initialized();
  }
}


/* *******************************************************************
* Initialization
******************************************************************* */
void Amanzi_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // initialization using base class
  Chemistry_PK::Initialize(S);

  Teuchos::RCP<Epetra_MultiVector> tcc = 
      S->GetFieldData("total_component_concentration", passwd_)->ViewComponent("cell", true);

  XMLParameters();

  // TODO: some sort of check of the state object to see if mineral_ssa,
  // CEC, site density, etc is present.

  // initial conditions for minerals etc should be handled by the
  // state/chemistry_state object before we reach this point. We just
  // resize our local memory for migrating data here.

  SizeBeakerStructures_();

  // copy the first cell data into the beaker storage for
  // initialization purposes
  CopyCellStateToBeakerStructures(0, tcc);

  // finish setting up & testing the chemistry object
  int ierr(0);
  std::string internal_msg;
  try {
    vo_->Write(Teuchos::VERB_HIGH, "Initializing chemistry in cell 0...\n");
    chem_->Setup(beaker_components_, beaker_parameters_);
    chem_->Display();

    // solve for initial free-ion concentrations
    vo_->Write(Teuchos::VERB_HIGH, "Initial speciation calculations in cell 0...\n");
    chem_->Speciate(&beaker_components_, beaker_parameters_);

    vo_->Write(Teuchos::VERB_HIGH, "\nTest solution of initial conditions in cell 0:\n");
    chem_->DisplayResults();
  } catch (ChemistryException& geochem_err) {
    ierr = 1;
    internal_msg = geochem_err.message_;
  }

  ErrorAnalysis(ierr, internal_msg);

  // TODO(bandre): at this point we should know about any additional
  // storage that chemistry needs...
  AllocateAdditionalChemistryStorage_(beaker_components_);

  SetupAuxiliaryOutput();

  // solve for initial free-ion concentrations
  int num_cells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ierr = 0;
  for (int c = 0; c < num_cells; ++c) {
    CopyCellStateToBeakerStructures(c, tcc);

    try {
      chem_->Speciate(&beaker_components_, beaker_parameters_);
      CopyBeakerStructuresToCellState(c, tcc);

    } catch (ChemistryException& geochem_err) {
      ierr = 1;
      internal_msg = geochem_err.message_;
    }
  }

  // figure out if any of the processes threw an error, if so all processes will re-throw
  ErrorAnalysis(ierr, internal_msg);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initalization of PK was successful, T=" 
        << S->time() << vo_->reset() << std::endl << std::endl;
  }
}


/* *******************************************************************
* Initialization helper functions
******************************************************************* */
void Amanzi_PK::XMLParameters()
{
  // thermo file name and format, then create the database!
  if (cp_list_->isSublist("Thermodynamic Database")) {
    Teuchos::ParameterList& tdb_list_ = cp_list_->sublist("Thermodynamic Database");

    // currently we only support the simple format.
    if (tdb_list_.isParameter("format")) {
      std::string database_format = tdb_list_.get<std::string>("format");
      if (database_format == "simple") {
        chem_ = new SimpleThermoDatabase();
      } else {
        // invalid database format
        std::ostringstream msg;
        msg << ChemistryException::kChemistryError;
        msg << "Amanzi_PK::XMLParameters(): \n";
        msg << "  In sublist 'Thermodynamic Database', the parameter 'format' must be 'simple'.\n";
        Exceptions::amanzi_throw(ChemistryInvalidInput(msg.str()));  
      }
    } else {
      // invalid database format
      std::ostringstream msg;
      msg << ChemistryException::kChemistryError;
      msg << "Amanzi_PK::XMLParameters(): \n";
      msg << "  In sublist 'Thermodynamic Database', the parameter 'format' must be specified.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(msg.str()));
    }

    beaker_parameters_ = chem_->GetDefaultParameters();

    // get file name
    if (tdb_list_.isParameter("file")) {
      beaker_parameters_.thermo_database_file = tdb_list_.get<std::string>("file");
    } else {
      std::ostringstream msg;
      msg << ChemistryException::kChemistryError;
      msg << "Amanzi_PK::XMLParameters(): \n";
      msg << "  Input parameter 'file' in 'Thermodynamic Database' sublist must be specified.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(msg.str()));         
    }
  } else {
    std::ostringstream msg;
    msg << ChemistryException::kChemistryError;
    msg << "Amanzi_PK::XMLParameters(): \n";
    msg << "  'Thermodynamic Database' sublist must be specified.\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(msg.str()));    
  }

  // activity model
  beaker_parameters_.activity_model_name = cp_list_->get<std::string>("activity model", "unit");
  // -- Pitzer virial coefficients database
  if (beaker_parameters_.activity_model_name == "pitzer-hwm") {
    if (cp_list_->isParameter("Pitzer database file")) {
      beaker_parameters_.pitzer_database = cp_list_->get<std::string>("Pitzer database file");
    } else {
      std::ostringstream msg;
      msg << ChemistryException::kChemistryError;
      msg << "Amanzi_PK::XMLParameters():\n";
      msg << "  Input parameter 'Pitzer database file' must be specified if 'activity model' is 'pitzer-hwm'.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(msg.str()));
    }
  }

  // solver parameters
  beaker_parameters_.tolerance = cp_list_->get<double>("tolerance", 1.0e-12);
  beaker_parameters_.max_iterations = cp_list_->get<int>("maximum Newton iterations", 200);

  // auxiliary data
  aux_names_.clear();
  if (cp_list_->isParameter("auxiliary data")) {
    Teuchos::Array<std::string> names = cp_list_->get<Teuchos::Array<std::string> >("auxiliary data");
    for (Teuchos::Array<std::string>::const_iterator name = names.begin();
         name != names.end(); ++name) {
      if (*name == "pH") {
        aux_names_.push_back(*name);
      } else {
        std::stringstream message;
        message << "XMLParameters(): unknown value in 'auxiliary data' list: " 
                << *name << std::endl;
        vo_->WriteWarning(Teuchos::VERB_LOW, message);
      }
    }
  }

  // misc other chemistry flags
  max_time_step_ = cp_list_->get<double>("max time step (s)", 9.9e+9);
}


/* *******************************************************************
* Requires that Beaker::Setup() has already been called!
******************************************************************* */
void Amanzi_PK::SetupAuxiliaryOutput()
{
  // TODO: this indexing scheme will not be appropriate when
  // additional types of aux data are requested, e.g. mineral SI.....
  unsigned int nvars = aux_names_.size();
  std::string name;
  aux_index_.clear();
  for (unsigned int i = 0; i < nvars; i++) {
    if (aux_names_.at(i) == "pH") {
      name = "H+";
    } else {
      name = aux_names_.at(i);
    }
    int index = chem_->GetPrimaryIndex(name);
    if (index == -1) {
        // check to make sure it is not -1, an invalid name/index
      std::stringstream message;
      message << "ChemistryPK::SetupAuxiliaryOutput() : "
              << "Output was requested for '" << aux_names_.at(i) 
              << "' (" << name 
              << ") but no chemistry varibles of this name were found.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(message.str()));
    } else {
      aux_index_.push_back(index);
    }
  }

  // create the Epetra_MultiVector that will hold the data
  if (nvars > 0) {
    aux_data_ = Teuchos::rcp(new Epetra_MultiVector(mesh_->cell_map(false), nvars));
  } else {
    aux_data_ = Teuchos::null;
  }
}


/* *******************************************************************
* Initialize the beaker component data structure
******************************************************************* */
void Amanzi_PK::SizeBeakerStructures_()
{
  // NOTE: The beaker already has data for site density, sorption
  // isotherms, ssa. If we want to use that single global value, then
  // we leave these arrays empty as a flag to the beaker to use its
  // own value. If we want to over ride the global chemistry value
  // with cell by cell data, then we resize the containers here.

  beaker_components_.total.resize(number_aqueous_components_, 0.0);
  beaker_components_.free_ion.resize(number_aqueous_components_, 1.0e-9);
  beaker_components_.mineral_volume_fraction.resize(number_minerals_, 0.0);

  if (using_sorption_) {
    beaker_components_.total_sorbed.resize(number_total_sorbed_, 0.0);
  } else {
    beaker_components_.total_sorbed.clear();
  }

  if (number_minerals_ > 0) {
    beaker_components_.mineral_specific_surface_area.resize(number_minerals_, 0.0);
  } else {
    beaker_components_.mineral_specific_surface_area.clear();
  }

  if (number_ion_exchange_sites_ > 0) {
    beaker_components_.ion_exchange_sites.resize(number_ion_exchange_sites_, 0.0);
  } else {
    beaker_components_.ion_exchange_sites.clear();
  }

  if (number_sorption_sites_ > 0) {
    beaker_components_.surface_site_density.resize(number_sorption_sites_, 0.0);
  } else {
    beaker_components_.surface_site_density.clear();
  }

  if (using_sorption_isotherms_) {
    beaker_components_.isotherm_kd.resize(number_aqueous_components_, 0.0);
    beaker_components_.isotherm_freundlich_n.resize(number_aqueous_components_, 0.0);
    beaker_components_.isotherm_langmuir_b.resize(number_aqueous_components_, 0.0);
  } else {
    beaker_components_.isotherm_kd.clear();
    beaker_components_.isotherm_freundlich_n.clear();
    beaker_components_.isotherm_langmuir_b.clear();
  }
}


/* *******************************************************************
* We must use the aqueous totals value calculated from transport
* (aqueous_components), not the value stored in state!
******************************************************************* */
void Amanzi_PK::CopyCellStateToBeakerStructures(
    int cell_id, Teuchos::RCP<Epetra_MultiVector> aqueous_components)
{
  for (unsigned int i = 0; i < number_aqueous_components_; i++) {
    beaker_components_.total.at(i) = (*aqueous_components)[i][cell_id];
  }

  const Epetra_MultiVector& free_ion = *S_->GetFieldData("free_ion_species")->ViewComponent("cell", true);
  for (int i = 0; i < number_aqueous_components_; ++i) {
    beaker_components_.free_ion.at(i) = free_ion[i][cell_id];
  }

  // activity coefficients
  if (beaker_components_.primary_activity_coeff.size() > 0) {
    const Epetra_MultiVector& activity = *S_->GetFieldData("primary_activity_coeff")->ViewComponent("cell", true);
    for (int i = 0; i < beaker_components_.primary_activity_coeff.size(); ++i) {
      beaker_components_.primary_activity_coeff.at(i) = activity[i][cell_id];
    }
  }

  if (beaker_components_.secondary_activity_coeff.size() > 0) {
    const Epetra_MultiVector& activity = *S_->GetFieldData("secondary_activity_coeff")->ViewComponent("cell", true);
    for (int i = 0; i < beaker_components_.secondary_activity_coeff.size(); ++i) {
      beaker_components_.secondary_activity_coeff.at(i) = activity[i][cell_id];
    }
  }

  // minerals
  if (number_minerals_ > 0) {
    const Epetra_MultiVector& mineral_vf = *S_->GetFieldData("mineral_volume_fractions")->ViewComponent("cell", true);
    const Epetra_MultiVector& mineral_ssa = *S_->GetFieldData("mineral_specific_surface_area")->ViewComponent("cell", true);

    for (int i = 0; i < number_minerals_; ++i) {
      beaker_components_.mineral_volume_fraction[i] = mineral_vf[i][cell_id];
      beaker_components_.mineral_specific_surface_area.at(i) = mineral_ssa[i][cell_id];
    }
  }

  // general sorption storage
  if (using_sorption_) {
    const Epetra_MultiVector& sorbed = *S_->GetFieldData("total_sorbed")->ViewComponent("cell", true);
    for (int i = 0; i < number_aqueous_components_; ++i) {
      beaker_components_.total_sorbed.at(i) = sorbed[i][cell_id];
    }
  }

  // ion exchange
  // TODO: only allow one ion exchange site at the moment!
  if (number_ion_exchange_sites_ > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData("ion_exchange_sites")->ViewComponent("cell", true);
    for (unsigned int i = 0; i < number_ion_exchange_sites_; i++) {
      beaker_components_.ion_exchange_sites[i] = ion_exchange[i][cell_id];
      // TODO(bandre): need to save ion exchange ref cation conc here!
    }
  }
  
  if (beaker_components_.ion_exchange_ref_cation_conc.size() > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData("ion_exchange_ref_cation_conc")->ViewComponent("cell", true);
    for (unsigned int i = 0; i < beaker_components_.ion_exchange_ref_cation_conc.size(); ++i) {
      beaker_components_.ion_exchange_ref_cation_conc.at(i) = ion_exchange[i][cell_id];
    }
  }

  // surface complexation
  if (number_sorption_sites_ > 0) {
    const Epetra_MultiVector& sorption_sites = *S_->GetFieldData("sorption_sites")->ViewComponent("cell", true);

    for (int i = 0; i < number_sorption_sites_; ++i) {
      beaker_components_.surface_site_density[i] = sorption_sites[i][cell_id];
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }

  if (beaker_components_.surface_complex_free_site_conc.size() > 0) {
    const Epetra_MultiVector& surface_complex =
        *S_->GetFieldData("surface_complex_free_site_conc")->ViewComponent("cell", true);
    for (int i = 0; i < beaker_components_.surface_complex_free_site_conc.size(); ++i) {
      beaker_components_.surface_complex_free_site_conc.at(i) = surface_complex[i][cell_id];
    }
  }

  // sorption isotherms
  if (using_sorption_isotherms_) {
    const Epetra_MultiVector& isotherm_kd = *S_->GetFieldData("isotherm_kd")->ViewComponent("cell", true);
    const Epetra_MultiVector& isotherm_freundlich_n = *S_->GetFieldData("isotherm_freundlich_n")->ViewComponent("cell", true);
    const Epetra_MultiVector& isotherm_langmuir_b = *S_->GetFieldData("isotherm_langmuir_b")->ViewComponent("cell", true);

    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      beaker_components_.isotherm_kd.at(i) = isotherm_kd[i][cell_id];
      beaker_components_.isotherm_freundlich_n.at(i) = isotherm_freundlich_n[i][cell_id];
      beaker_components_.isotherm_langmuir_b.at(i) = isotherm_langmuir_b[i][cell_id];
    }
  }

  // copy data from state arrays into the beaker parameters
  const Epetra_MultiVector& porosity = *S_->GetFieldData("porosity")->ViewComponent("cell", true);
  const Epetra_MultiVector& water_saturation = *S_->GetFieldData("saturation_liquid")->ViewComponent("cell", true);
  double water_density = *S_->GetScalarData("fluid_density");

  beaker_parameters_.water_density = water_density;
  beaker_parameters_.porosity = porosity[0][cell_id];
  beaker_parameters_.saturation = water_saturation[0][cell_id];
  beaker_parameters_.volume = mesh_->cell_volume(cell_id);
}


/* *******************************************************************
* Copy data from the beaker back into the state arrays.
******************************************************************* */
void Amanzi_PK::CopyBeakerStructuresToCellState(
    int cell_id, Teuchos::RCP<Epetra_MultiVector> aqueous_components)
{
  for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
    (*aqueous_components)[i][cell_id] = beaker_components_.total.at(i);
  }

  const Epetra_MultiVector& free_ion = *S_->GetFieldData("free_ion_species")->ViewComponent("cell", true);
  for (int i = 0; i < number_aqueous_components_; ++i) {
    free_ion[i][cell_id] = beaker_components_.free_ion.at(i);
  }

  // activity coefficients
  if (beaker_components_.primary_activity_coeff.size() > 0) {
    const Epetra_MultiVector& activity = *S_->GetFieldData("primary_activity_coeff")->ViewComponent("cell", true);
    for (int i = 0; i < beaker_components_.primary_activity_coeff.size(); ++i) {
      activity[i][cell_id] = beaker_components_.primary_activity_coeff.at(i);
    }
  }

  if (beaker_components_.secondary_activity_coeff.size() > 0) {
    const Epetra_MultiVector& activity = *S_->GetFieldData("secondary_activity_coeff")->ViewComponent("cell", true);
    for (int i = 0; i < beaker_components_.secondary_activity_coeff.size(); ++i) {
      activity[i][cell_id] =  beaker_components_.secondary_activity_coeff.at(i);
    }
  }

  // minerals
  if (number_minerals_ > 0) {
    const Epetra_MultiVector& mineral_vf = *S_->GetFieldData("mineral_volume_fractions")->ViewComponent("cell", true);
    const Epetra_MultiVector& mineral_ssa = *S_->GetFieldData("mineral_specific_surface_area")->ViewComponent("cell", true);

    for (int i = 0; i < number_minerals_; ++i) {
      mineral_vf[i][cell_id] = beaker_components_.mineral_volume_fraction.at(i);
      mineral_ssa[i][cell_id] = beaker_components_.mineral_specific_surface_area.at(i);
    }
  }

  // sorption
  if (using_sorption_) {
    const Epetra_MultiVector& sorbed = *S_->GetFieldData("total_sorbed")->ViewComponent("cell", true);
    for (int i = 0; i < number_aqueous_components_; ++i) {
      sorbed[i][cell_id] = beaker_components_.total_sorbed.at(i);
    }
  }

  // surface complexation
  if (number_sorption_sites_ > 0) {
    const Epetra_MultiVector& sorption_sites = *S_->GetFieldData("sorption_sites")->ViewComponent("cell", true);
    for (int i = 0; i < number_sorption_sites_; i++) {
      sorption_sites[i][cell_id] = beaker_components_.surface_site_density.at(i);
      // TODO: need to save surface complexation free site conc here!
    }
  }

  if (beaker_components_.surface_complex_free_site_conc.size() > 0) {
    const Epetra_MultiVector& surface_complex =
        *S_->GetFieldData("surface_complex_free_site_conc")->ViewComponent("cell", true);
    for (int i = 0; i < beaker_components_.surface_complex_free_site_conc.size(); ++i) {
      surface_complex[i][cell_id] = beaker_components_.surface_complex_free_site_conc.at(i);
    }
  }

  // ion exchange
  if (number_ion_exchange_sites_ > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData("ion_exchange_sites")->ViewComponent("cell", true);
    for (int i = 0; i < number_ion_exchange_sites_; ++i) {
      ion_exchange[i][cell_id] = beaker_components_.ion_exchange_sites.at(i);
      // TODO(bandre): need to save ion exchange ref cation conc here!
    }
  }

  if (beaker_components_.ion_exchange_ref_cation_conc.size() > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData("ion_exchange_ref_cation_conc")->ViewComponent("cell", true);
    for (int i = 0; i < beaker_components_.ion_exchange_ref_cation_conc.size(); ++i) {
      ion_exchange[i][cell_id] = beaker_components_.ion_exchange_ref_cation_conc.at(i);
    }
  }

  // sorption isotherms
  if (using_sorption_isotherms_) {
    const Epetra_MultiVector& isotherm_kd = *S_->GetFieldData("isotherm_kd")->ViewComponent("cell", true);
    const Epetra_MultiVector& isotherm_freundlich_n = *S_->GetFieldData("isotherm_freundlich_n")->ViewComponent("cell", true);
    const Epetra_MultiVector& isotherm_langmuir_b = *S_->GetFieldData("isotherm_langmuir_b")->ViewComponent("cell", true);

    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      isotherm_kd[i][cell_id] = beaker_components_.isotherm_kd.at(i);
      isotherm_freundlich_n[i][cell_id] = beaker_components_.isotherm_freundlich_n.at(i);
      isotherm_langmuir_b[i][cell_id] = beaker_components_.isotherm_langmuir_b.at(i);
    }
  }

  // TODO(bandre): if chemistry can modify the porosity or density,
  // then they should be updated here!
}


/* ******************************************************************
* This function advances concentrations in the auxialiry vector 
* aqueous_components_ (defined in the base class). This vector must 
* be set up using routine set_aqueous_components(). Tipically, it 
* contains values advected by the transport PK.
******************************************************************* */
bool Amanzi_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool failed(false);
  std::string internal_msg;

  double dt = t_new - t_old;
  current_time_ = saved_time_ + dt;

  int num_itrs, max_itrs(0), min_itrs(10000000), avg_itrs(0);
  int cmax(-1), cmin(-1), ierr(0);

  int num_cells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < num_cells; ++c) {
    CopyCellStateToBeakerStructures(c, aqueous_components_);
    try {
      // create a backup copy of the components
      chem_->CopyComponents(beaker_components_, &beaker_components_copy_);

      // chemistry computations for this cell
      num_itrs = chem_->ReactionStep(&beaker_components_, beaker_parameters_, dt);

      if (max_itrs < num_itrs) {
        max_itrs = num_itrs;
        cmax = c;
      }
      if (min_itrs > num_itrs) {
        min_itrs = num_itrs;
        cmin = c;
      }
      avg_itrs += num_itrs;
    } catch (ChemistryException& geochem_err) {
      ierr = 1;
      internal_msg = geochem_err.kChemistryError;
      break;
    }

    if (ierr == 0) CopyBeakerStructuresToCellState(c, aqueous_components_);
    // TODO: was porosity etc changed? copy someplace
  }

  ErrorAnalysis(ierr, internal_msg);
  
  std::stringstream ss;
  ss << "Newton iterations: " << min_itrs << "/" << max_itrs << "/" 
     << avg_itrs / num_cells << ", maximum in gid=" << mesh_->GID(cmax, AmanziMesh::CELL) << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, ss.str());

  // dumping the values of the final cell. not very helpful by itself,
  // but can be move up into the loops....
  // chem_->DisplayTotalColumnHeaders(true);
  // chem_->DisplayTotalColumns(current_time_, beaker_components_, true);

  return failed;
}


/* ******************************************************************
* The MPC will call this function to signal to the process kernel 
* that it has accepted the state update, thus, the PK should update
* possible auxilary state variables here
******************************************************************* */
void Amanzi_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  saved_time_ = t_new;

  // debug output
  // chem_->Speciate(&beaker_components_, beaker_parameters_);
  // chem_->DisplayResults();
  // chem_->DisplayTotalColumnHeaders(bool display_free_columns);
  // chem_->DisplayTotalColumns(saved_time_, beaker_components_, true);
}


/* ******************************************************************
*
******************************************************************* */
Teuchos::RCP<Epetra_MultiVector> Amanzi_PK::extra_chemistry_output_data() {
  if (aux_data_ != Teuchos::null) {
    const Epetra_MultiVector& free_ion = *S_->GetFieldData("free_ion_species")->ViewComponent("cell", true);
    const Epetra_MultiVector& activity = *S_->GetFieldData("primary_activity_coeff")->ViewComponent("cell", true);
    int num_cells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

    for (int cell = 0; cell < num_cells; cell++) {
      // populate aux_data_ by copying from the appropriate internal storage
      for (unsigned int i = 0; i < aux_names_.size(); i++) {
        if (aux_names_.at(i) == "pH") {
          double* cell_aux_data = (*aux_data_)[i];
          double* cell_free_ion = free_ion[aux_index_.at(i)];
          double* activity_coeff = activity[aux_index_.at(i)];
          double activity = cell_free_ion[cell] * activity_coeff[cell];
          cell_aux_data[cell] = -std::log10(activity);
        } else {
          // don't support anything else at this time....
        }
      }
    }
  }
  return aux_data_;
}


/* ******************************************************************
*
******************************************************************* */
void Amanzi_PK::set_chemistry_output_names(std::vector<std::string>* names) {
  names->clear();
  for (std::vector<std::string>::const_iterator name = aux_names_.begin();
       name != aux_names_.end(); name++) {
    names->push_back(*name);
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
