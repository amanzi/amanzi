#include <sstream>
#include <string>
#include <algorithm>

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputParserIS.hh"
#include "InputParserIS_Defs.hh"
#include "XMLParameterListWriter.hh"

namespace Amanzi {
namespace AmanziInput {

/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputParserIS::Translate(Teuchos::ParameterList* input_plist, int numproc)
{
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::rcp(input_plist, false);
  numproc_ = numproc;

  // first make sure the version is correct
  CheckAmanziInputVersion_(plist);
  InitGlobalInfo_(plist);

  // unstructured header
  Teuchos::ParameterList new_list, tmp_list, pks_list, cd_list;
  
  new_list.set<bool>("Native Unstructured Input", true);
  new_list.set<std::string>("grid_option", "Unstructured");
  new_list.set<std::string>("input file name",
                            plist->get<std::string>("input file name", "unit_test.xml"));

  // checkpoint list is optional
  tmp_list = CreateCheckpointDataList_(plist);
  if (tmp_list.begin() != tmp_list.end()) {
    new_list.sublist("Checkpoint Data") = tmp_list;
  }

  // walkabout list is optional
  tmp_list = CreateWalkaboutDataList_(plist);
  if (tmp_list.begin() != tmp_list.end()) {
    new_list.sublist("Walkabout Data") = tmp_list;
  }

  // visualization list is optional
  tmp_list = CreateVisualizationDataList_(plist);
  if (tmp_list.begin() != tmp_list.end()) {
    new_list.sublist("Visualization Data") = tmp_list;
  }

  if (plist->sublist("Output").isSublist("Observation Data")) {
    Teuchos::ParameterList& od_list = plist->sublist("Output").sublist("Observation Data");
    if (od_list.begin() != od_list.end()) {
      new_list.sublist("Observation Data") = CreateObservationDataList_(plist);
    }
  }

  // mesh and geometry
  new_list.sublist("Regions") = CopyRegionsList_(plist);
  new_list.sublist("Mesh") = CreateMeshList_(plist);
  new_list.sublist("Domain") = CopyDomainList_(plist);

  // cycle driver
  cd_list = CreateCycleDriverList_(plist);
  new_list.sublist("Cycle Driver") = cd_list;
  CreatePKslist_(cd_list, pks_list);

  FillPKslist_(plist, pks_list);

  new_list.sublist("State") = CreateStateList_(plist);
  new_list.sublist("Preconditioners") = CreatePreconditionersList_(plist);
  new_list.sublist("Solvers") = CreateSolversList_(plist);

  new_list.sublist("PKs") = pks_list;

  // analysis list
  new_list.sublist("Analysis") = CreateAnalysisList_();

  // hack (additional transport diagnostics)
  new_list.sublist("PKs").sublist("Transport")
          .set<Teuchos::Array<std::string> >("runtime diagnostics: regions", transport_diagnostics_);

  // output unused parameters (experimental)
  if (vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    // PrintUnused_(plist->sublist("Execution Control").sublist("Numerical Control Parameters"), vo_);
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("red") << "\n\nEND-of-LINE Statement.\n\n"
               << "Series 1.2.x will be terminated in the next release of Amanzi.\n"
               << "Please upgrade to spec 2.2 or higher to use new features." << vo_->reset() << "\n\n\n";
  }

  return new_list;
}


/* ******************************************************************
* Verify that we use XML file with the correct amanzi version.
****************************************************************** */
void InputParserIS::CheckAmanziInputVersion_(Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  std::string version = plist->get<std::string>("Amanzi Input Format Version", "FAIL");
  if (version == "FAIL") {
    Exceptions::amanzi_throw(Errors::Message("The input file does not specify an \"Amanzi Input Format Version\""));
  }

  int major, minor, micro;

  std::stringstream ss;
  ss << version;
  std::string ver;

  try {
    getline(ss,ver,'.');
    major = boost::lexical_cast<int>(ver);

    getline(ss,ver,'.');
    minor = boost::lexical_cast<int>(ver);

    getline(ss,ver);
    micro = boost::lexical_cast<int>(ver);
  }
  catch (...) {
    Exceptions::amanzi_throw(Errors::Message("The version string in the input file '"+version+"' has the wrong format, please use X.Y.Z, where X, Y, and Z are integers."));
  }

  if ((major != AMANZI_OLD_INPUT_VERSION_MAJOR) ||
      (minor != AMANZI_OLD_INPUT_VERSION_MINOR) ||
      ((micro != AMANZI_OLD_INPUT_VERSION_MICRO) && micro != (AMANZI_OLD_INPUT_VERSION_MICRO + 1)) ) {
    std::stringstream ss_ver_reqd;
    ss_ver_reqd << AMANZI_OLD_INPUT_VERSION_MAJOR << "." << AMANZI_OLD_INPUT_VERSION_MINOR << "." << AMANZI_OLD_INPUT_VERSION_MICRO;
    std::stringstream ss_ver_inp;
    ss_ver_inp << major << "." << minor << "." << micro;

    Exceptions::amanzi_throw(Errors::Message("The input format version "+ss_ver_inp.str()+" does not match the required version "+ss_ver_reqd.str()));
  }
}


/* ******************************************************************
* Initizialize some global information.
****************************************************************** */
void InputParserIS::InitGlobalInfo_(Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  Errors::Message msg;

  // spatial dimension
  if (plist->isSublist("Domain")) {
    spatial_dimension_ = plist->sublist("Domain").get<int>("Spatial Dimension",0);
  } else {
    spatial_dimension_ = 0;
  }

  // create an all region
  if (!plist->sublist("Regions").isSublist("All")) {
    Teuchos::ParameterList& allreg = plist->sublist("Regions").sublist("All")
        .sublist("Region: Box");

    Teuchos::Array<double> low, high;
    low.push_back(-1e99);
    high.push_back(1e99);

    if (spatial_dimension_ >= 2) {
      low.push_back(-1e99);
      high.push_back(1e99);
    }

    if (spatial_dimension_ == 3) {
      low.push_back(-1e99);
      high.push_back(1e99);
    }

    allreg.set<Teuchos::Array<double> >("Low Coordinate", low);
    allreg.set<Teuchos::Array<double> >("High Coordinate", high);
  }

  // define verbosity
  if (plist->isSublist("Execution Control")) {
    std::string verbosity = plist->sublist("Execution Control").get<std::string>("Verbosity", VERBOSITY_DEFAULT);

    if (verbosity == "None" || verbosity == "none") {
      verbosity_level = "none";
    } else if (verbosity == "Low" || verbosity == "low") {
      verbosity_level = "low";
    } else if (verbosity == "Medium" ||verbosity == "medium") {
      verbosity_level = "medium";
    } else if (verbosity == "High" || verbosity == "high") {
      verbosity_level = "high";
    } else if (verbosity == "Extreme" || verbosity == "extreme") {
      verbosity_level = "extreme";
    } else {
      Exceptions::amanzi_throw(Errors::Message("Verbosity must be one of None, Low, Medium, High, or Extreme."));
    }
  }

  Teuchos::ParameterList vo_list;
  vo_list.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);
  vo_ = new VerboseObject("InputParser1.2.3", vo_list); 

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    *vo_->os() << "verbosity level=" <<  vo_->getVerbLevel()
               << " (low=" << Teuchos::VERB_LOW 
               << ", medium, high, extreme=" << Teuchos::VERB_EXTREME << ")" << std::endl;
  }

  // check if Transport is Off
  std::string transport_model = plist->sublist("Execution Control").get<std::string>("Transport Model");
  std::string chemistry_model = plist->sublist("Execution Control").get<std::string>("Chemistry Model");

  phases_.resize(3);
  phases_[0].name = "Aqueous";
  phases_[1].name = "Gaseous";

  // don't know the history of these variables, clear them just to be safe.
  comp_names_.clear();
  mineral_names_.clear();
  sorption_site_names_.clear();

  Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");
  Teuchos::ParameterList::ConstIterator item;

  for (item = phase_list.begin(); item != phase_list.end(); ++item) {
    std::string phase_name = phase_list.name(item);
    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      *vo_->os() << "Found phase: " << phase_name << std::endl;
    }

    if (transport_model != "Off" || chemistry_model != "Off") {
      int i(0);
      if (phase_name == "Gaseous") i = 1;
      if (phase_name == "Solid") i = 2;
      phases_[i].name = phase_name;

      Teuchos::ParameterList& phase_sublist = phase_list.sublist(phase_name);

      if (phase_sublist.isSublist("Phase Components")) {
        Teuchos::ParameterList& phase_components = phase_sublist.sublist("Phase Components");
        // for now there should only be one sublist here, we allow it to be named something
        // the user chooses, e.g. Water
        Teuchos::ParameterList::ConstIterator pcit = phase_components.begin();
        ++pcit;
        if (pcit != phase_components.end()) {
          Exceptions::amanzi_throw(Errors::Message("Currently Amanzi only supports one phase component, e.g. Water"));
        }
        pcit = phase_components.begin();
        if (!pcit->second.isList()) {
          msg << "The Phase Components list must only have one sublist, but you have specified"
              << " a parameter instead.";
          Exceptions::amanzi_throw(msg);
        }
        phases_[i].solute_name = pcit->first;
        Teuchos::ParameterList& components = phase_components.sublist(phases_[i].solute_name);
        if (components.isParameter("Component Solutes")) {
          phases_[i].solute_comp_names = components.get<Teuchos::Array<std::string> >("Component Solutes").toVector();
        }
      }
    }

    if (chemistry_model != "Off") {
      if (phase_list.name(item) == "Solid") {
        Teuchos::ParameterList solid_list = phase_list.sublist("Solid");
        // this is the order that the chemistry expects
        if (solid_list.isParameter("Minerals")) {
          mineral_names_ = solid_list.get<Teuchos::Array<std::string> >("Minerals");
        }
        if (solid_list.isParameter("Sorption Sites")) {
          sorption_site_names_ = solid_list.get<Teuchos::Array<std::string> >("Sorption Sites");
        }
      }  // end Solid phase
    }

    if ((phase_list.name(item) != "Aqueous" ) && 
        (phase_list.name(item) != "Solid") && 
        (phase_list.name(item) != "Gaseous")) {
      std::stringstream message;
      msg << "Error: InputParserIS: The only phases supported on unstructured"
          << " meshes at this time are 'Aqueous', 'Gaseous', and and 'Solid'.\n";
      Exceptions::amanzi_throw(msg);
    }
  }

  // create lists of components
  comp_names_ = phases_[0].solute_comp_names;
  comp_names_all_ = phases_[0].solute_comp_names;
  comp_names_all_.insert(comp_names_all_.end(),
      phases_[1].solute_comp_names.begin(), phases_[1].solute_comp_names.end());

  // dispersion (this is going to be used to translate to the transport list as well as the state list)
  // check if we need to write a dispersivity sublist
  need_dispersion_ = false;
  if (plist->isSublist("Material Properties")) {
    for (Teuchos::ParameterList::ConstIterator it = plist->sublist("Material Properties").begin();
         it != plist->sublist("Material Properties").end(); ++it) {
      if ((it->second).isList()) {
        Teuchos::ParameterList& mat_sublist = plist->sublist("Material Properties").sublist(it->first);
        for (Teuchos::ParameterList::ConstIterator jt = mat_sublist.begin(); jt != mat_sublist.end(); ++jt) {
          if ((jt->second).isList()) {
            const std::string pname = jt->first;
            if (pname.find("Dispersion Tensor") == 0 ||
                pname.find("Tortuosity") == 0) {
              need_dispersion_ = true;
            }
          }
        }
      }
    }
  }
}


/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::Array<std::string> InputParserIS::TranslateForms_(Teuchos::Array<std::string>& forms)
{
  Teuchos::Array<std::string> target_forms;

  for (Teuchos::Array<std::string>::const_iterator i = forms.begin();
       i != forms.end(); ++i) {

    if (*i == "Constant") {
      target_forms.push_back("constant");
    } else if (*i == "Linear") {
      target_forms.push_back("linear");
    } else {
      Exceptions::amanzi_throw(Errors::Message("Cannot translate the tabular function form "+*i));
    }
  }
  return target_forms;
}


/* ******************************************************************
* Create verbosity list.
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateVerbosityList_(const std::string& vlevel)
{
  Teuchos::ParameterList vlist;

  if (vlevel == "low") {
    vlist.set<std::string>("Verbosity Level","low");
  } else if (vlevel == "medium") {
    vlist.set<std::string>("Verbosity Level","medium");
  } else if (vlevel == "high") {
    vlist.set<std::string>("Verbosity Level","high");
  } else if (vlevel == "extreme") {
    vlist.set<std::string>("Verbosity Level","extreme");
  } else if (vlevel == "none") {
    vlist.set<std::string>("Verbosity Level","none");
  }

  return vlist;
}


/* ******************************************************************
* Analysis list can used by special tools.
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateAnalysisList_()
{
  Teuchos::ParameterList alist;
  alist.set<Teuchos::Array<std::string> >("used boundary condition regions", vv_bc_regions);
  alist.set<Teuchos::Array<std::string> >("used source regions", vv_src_regions);
  alist.set<Teuchos::Array<std::string> >("used observation regions", vv_obs_regions);

  alist.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);

  return alist;
}


/* ******************************************************************
* Analysis of unused parameters
****************************************************************** */
void InputParserIS::PrintUnused_(const Teuchos::ParameterList& p, VerboseObject* vo) const
{
  Teuchos::OSTab tab = vo_->getOSTab();
  // print parameters first
  for (Teuchos::ParameterList::ConstIterator i = p.begin(); i != p.end(); ++i) {
    if (!(p.entry(i).isUsed())) {
      *vo->os() << vo_->color("yellow") << "Unused: \"" << p.name(i) << "\"" << vo_->reset() << std::endl;
    }
  }
  
  // print sublists second
  for (Teuchos::ParameterList::ConstIterator i = p.begin(); i != p.end(); ++i) {
    const Teuchos::ParameterEntry& entry_i = p.entry(i);
    if (!entry_i.isList()) continue;

    const std::string& docString = entry_i.docString();
    const std::string& name_i = p.name(i);
    PrintUnused_(Teuchos::getValue<Teuchos::ParameterList>(entry_i), vo);
  }
}

}  // namespace AmanziInput
}  // namespace Amanzi
