/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Erin Barker (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <fstream>

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "InputConverterU.hh"
#include "XMLParameterListWriter.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Main driver for the new translator.
****************************************************************** */
Teuchos::ParameterList InputConverterU::Translate(int rank, int num_proc)
{
  rank_ = rank;
  num_proc_ = num_proc;
  Teuchos::ParameterList out_list;
  
  // grab verbosity early
  verb_list_ = TranslateVerbosity_();
  Teuchos::ParameterList tmp_list(verb_list_);
  vo_ = new VerboseObject("InputConverter", tmp_list);
  Teuchos::OSTab tab = vo_->getOSTab();

  // checks that input XML is structurally sound
  VerifyXMLStructure_();

  // checks that input XML has valid version
  ParseVersion_();

  // parsing of miscalleneous lists
  ParseSolutes_();
  ParseConstants_();
  ModifyDefaultPhysicalConstants_();
  ParseModelDescription_();

  out_list.set<bool>("Native Unstructured Input", "true");

  out_list.sublist("units") = TranslateUnits_();  
  out_list.sublist("mesh") = TranslateMesh_();
  out_list.sublist("domain").set<int>("spatial dimension", dim_);
  out_list.sublist("regions") = TranslateRegions_();

  out_list.sublist("state") = TranslateState_();

  const Teuchos::ParameterList& tmp = TranslateOutput_();
  for (auto it = tmp.begin(); it != tmp.end(); ++it) {
    out_list.sublist(it->first) = tmp.sublist(it->first);
  }

  out_list.sublist("cycle driver") = TranslateCycleDriver_();
  Teuchos::ParameterList& cd_list = out_list.sublist("cycle driver");
  out_list.sublist("PKs") = TranslatePKs_(cd_list);

  out_list.sublist("solvers") = TranslateSolvers_();
  out_list.sublist("preconditioners") = TranslatePreconditioners_();

  // analysis list
  out_list.sublist("analysis") = CreateAnalysis_();
  FilterEmptySublists_(out_list);

  // post-processing (some may go away)
  FinalizeMPC_PKs_(out_list);
  MergeInitialConditionsLists_(out_list);

  // miscalleneous cross-list information
  // -- initialization file name
  if (init_filename_.size() > 0) {
    out_list.sublist("state").set<std::string>("initialization filename", init_filename_);
  }

  // -- additional saturated flow fields
  if (pk_model_["flow"] == "darcy") {
    Teuchos::Array<std::string> regions(1, "All");
    out_list.sublist("state").sublist("initial conditions").sublist("saturation_liquid")
        .sublist("function").sublist("ALL")
        .set<Teuchos::Array<std::string> >("regions", regions)
        .set<std::string>("component", "cell")
        .sublist("function").sublist("function-constant")
        .set<double>("value", 1.0);

    if (fracture_regions_.size() > 0) {
      out_list.sublist("state").sublist("initial conditions").sublist("fracture-saturation_liquid")
        .sublist("function").sublist("ALL")
        .set<Teuchos::Array<std::string> >("regions", regions)
        .set<std::string>("component", "cell")
        .sublist("function").sublist("function-constant")
        .set<double>("value", 1.0);
    }
  }

  // -- additional transport diagnostics
  if (transport_diagnostics_.size() > 0) {
    out_list.sublist("PKs").sublist("transport")
        .set<Teuchos::Array<std::string> >("runtime diagnostics: regions", transport_diagnostics_);
  }

  // -- walkabout enforces non-contiguous maps
  if (io_walkabout_) {
    out_list.sublist("mesh").sublist("unstructured").sublist("expert")
        .set<bool>("contiguous global ids", false);
  }

  // -- single region for fracture network
  if (fracture_regions_.size() > 0) {
    // std::string network = CreateUniqueName_(fracture_regions_);
    std::string network("FRACTURE_NETWORK_INTERNAL");
    out_list.sublist("regions").sublist(network).sublist("region: logical")
        .set<std::string>("operation", "union")
        .set<Teuchos::Array<std::string> >("regions", fracture_regions_);

    if (out_list.isSublist("visualization data")) {
      Teuchos::ParameterList& aux = out_list.sublist("visualization data");
      out_list.sublist("visualization data fracture") = aux;

      std::string name = aux.get<std::string>("file name base");
      aux.set<std::string>("file name base", name + "_matrix");
      out_list.sublist("visualization data fracture")
          .set<std::string>("file name base", name + "_fracture");
    }

    if (io_mesh_info_) {
      out_list.sublist("mesh info fracture") = out_list.sublist("mesh info");
      std::string filename = out_list.sublist("mesh info").get<std::string>("filename");
      out_list.sublist("mesh info fracture").set<std::string>("filename", filename + "_fracture");
    }
  }

  // -- final I/O
  PrintStatistics_();

  // save the translated file
  if (rank_ == 0) SaveXMLFile(out_list, xmlfilename_);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "CPU time stamp: " << vo_->clock() << std::endl;

  return out_list;
}
  

/* ******************************************************************
* Check that XML has required objects that frequnetly used.
****************************************************************** */
void InputConverterU::VerifyXMLStructure_()
{
  MemoryManager mm;

  std::vector<std::string> names;
  names.push_back("execution_controls");
  names.push_back("materials");
  names.push_back("process_kernels");
  names.push_back("phases");
  names.push_back("mesh");

  for (std::vector<std::string>::iterator it = names.begin(); it != names.end(); ++it) {
    DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode(it->c_str()));
    IsEmpty(node_list, *it); 
  }
}


/* ******************************************************************
* Extract information of solute components.
****************************************************************** */
void InputConverterU::ParseSolutes_()
{
  bool flag;
  char* tagname;
  char* text_content;

  MemoryManager mm;

  DOMNode* node;
  DOMNode* knode = doc_->getElementsByTagName(mm.transcode("phases"))->item(0);

  // liquid phase (try solutes, then primaries)
  std::string species("solute");
  node = GetUniqueElementByTagsString_(knode, "liquid_phase, dissolved_components, solutes", flag);
  if (!flag) {
    node = GetUniqueElementByTagsString_(knode, "liquid_phase, dissolved_components, primaries", flag);
    species = "primary";
  }

  DOMNodeList* children = node->getChildNodes();
  int nchildren = children->getLength();

  for (int i = 0; i < nchildren; ++i) {
    DOMNode* inode = children->item(i);
    tagname = mm.transcode(inode->getNodeName());
    std::string name = TrimString_(mm.transcode(inode->getTextContent()));

    if (species == tagname) {
      phases_["water"].push_back(name);

      DOMElement* element = static_cast<DOMElement*>(inode);
      // Polyethylene glycol has a molar mass of 1,000,000 g/mol
      double mol_mass = GetAttributeValueD_(element, "molar_mass", TYPE_NUMERICAL, 0.0, 1000.0, "kg/mol", false, 1.0);
      solute_molar_mass_[name] = mol_mass;
    }
  }
  
  comp_names_all_ = phases_["water"];

  // gas phase
  node = GetUniqueElementByTagsString_(knode, "gas_phase, dissolved_components, solutes", flag);
  if (flag) {
    children = node->getChildNodes();
    nchildren = children->getLength();

    for (int i = 0; i < nchildren; ++i) {
      DOMNode* inode = children->item(i);
      tagname = mm.transcode(inode->getNodeName());
      text_content = mm.transcode(inode->getTextContent());

      if (strcmp(tagname, "solute") == 0) {
        phases_["air"].push_back(TrimString_(text_content));
      }
    }

    comp_names_all_.insert(comp_names_all_.end(), phases_["air"].begin(), phases_["air"].end());
  }

  // output
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int nsolutes = phases_["water"].size();
    *vo_->os() << "Phase 'water' has " << nsolutes << " solutes\n";
    for (int i = 0; i < nsolutes; ++i) {
      *vo_->os() << " solute: " << phases_["water"][i] << std::endl;
    }
  }
}


/* ******************************************************************
* Adjust default values for physical constants
****************************************************************** */
void InputConverterU::ModifyDefaultPhysicalConstants_()
{
  // gravity magnitude
  const_gravity_ = GRAVITY_MAGNITUDE;
  auto it = constants_.find("gravity");
  if (it != constants_.end()) const_gravity_ = std::strtod(it->second.c_str(), NULL);

  // reference atmospheric pressure
  const_atm_pressure_ = ATMOSPHERIC_PRESSURE;
  it = constants_.find("atmospheric_pressure");
  if (it != constants_.end()) const_atm_pressure_ = std::strtod(it->second.c_str(), NULL);
}


/* ******************************************************************
* Extract generic verbosity object for all sublists.
****************************************************************** */
void InputConverterU::ParseModelDescription_()
{
  MemoryManager mm;
  DOMNodeList* node_list;
  DOMNode* node;

  bool flag;
  node_list = doc_->getElementsByTagName(mm.transcode("model_description"));
  node = GetUniqueElementByTagsString_(node_list->item(0), "coordinate_system", flag);

  if (flag) {
    coords_ = CharToStrings_(mm.transcode(node->getTextContent()));
  } else { 
    coords_.push_back("x"); 
    coords_.push_back("y"); 
    coords_.push_back("z"); 
  }

  node = GetUniqueElementByTagsString_(node_list->item(0), "author", flag);
  if (flag && vo_->getVerbLevel() >= Teuchos::VERB_HIGH) 
    *vo_->os() << "AUTHOR: " << mm.transcode(node->getTextContent()) << std::endl;
}


/* ******************************************************************
* Extract generic verbosity object for all sublists.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateVerbosity_()
{
  Teuchos::ParameterList vlist;

  MemoryManager mm;

  DOMNodeList* node_list;
  DOMNode* node_attr;
  DOMNamedNodeMap* attr_map;
  char* text_content;
 
  // get execution contorls node
  node_list = doc_->getElementsByTagName(mm.transcode("execution_controls"));
  
  for (int i = 0; i < node_list->getLength(); i++) {
    DOMNode* inode = node_list->item(i);

    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      DOMNodeList* children = inode->getChildNodes();

      for (int j = 0; j < children->getLength(); j++) {
        DOMNode* jnode = children->item(j);

        if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
          char* tagname = mm.transcode(jnode->getNodeName());
          if (std::string(tagname) == "verbosity") {
            attr_map = jnode->getAttributes();
            node_attr = attr_map->getNamedItem(mm.transcode("level"));
            if (node_attr) {
              text_content = mm.transcode(node_attr->getNodeValue());
              vlist.sublist("verbose object").set<std::string>("verbosity level", TrimString_(text_content));
              break;
            } else {
              ThrowErrorIllformed_("verbosity", "value", "level");
            }
          }
        }
      }
    }
  }
  return vlist;
}


/* ******************************************************************
* Units
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateUnits_()
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;  
  DOMNode* node;

  bool flag;
  std::string length("m"), time("s"), mass("kg"), concentration("molar"), temperature("K");

  node = GetUniqueElementByTagsString_("model_description, units, length_unit", flag);
  // if (flag) length = GetTextContentS_(node, "cm, in, ft, yd, m, km");
  if (flag) length = GetTextContentS_(node, "m");

  node = GetUniqueElementByTagsString_("model_description, units, time_unit", flag);
  // if (flag) time = GetTextContentS_(node, "s, min, h, d, y");
  if (flag) time = GetTextContentS_(node, "s");

  node = GetUniqueElementByTagsString_("model_description, units, mass_unit", flag);
  // if (flag) mass = GetTextContentS_(node, "g, lb, kg, ton");
  if (flag) mass = GetTextContentS_(node, "kg");

  node = GetUniqueElementByTagsString_("model_description, units, conc_unit", flag);
  if (flag) concentration = GetTextContentS_(node, "molar, SI, ppm, ppb");

  node = GetUniqueElementByTagsString_("model_description, units, temperature_unit", flag);
  // if (flag) temperature = GetTextContentS_(node, "K, C, F");
  if (flag) temperature = GetTextContentS_(node, "K");

  out_list.set<std::string>("length", length);
  out_list.set<std::string>("time", time);
  out_list.set<std::string>("mass", mass);
  out_list.set<std::string>("concentration", concentration);
  out_list.set<std::string>("amount", "mol");
  out_list.set<std::string>("temperature", temperature);
  
  // update default system units
  // units_.system().concentration = concentration;
  units_.Init(out_list);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating units: " << length << " " << time << " " 
               << mass << " " << concentration << " " << temperature << std::endl;

  return out_list;
}


/* ******************************************************************
* Analysis list can used by special tools.
****************************************************************** */
Teuchos::ParameterList InputConverterU::CreateAnalysis_()
{
  Teuchos::ParameterList out_list;
  auto& domain_list = out_list.sublist("domain");

  domain_list.set<Teuchos::Array<std::string> >("used boundary condition regions", vv_bc_regions_);
  domain_list.set<Teuchos::Array<std::string> >("used source regions", vv_src_regions_);
  domain_list.set<Teuchos::Array<std::string> >("used observation regions", vv_obs_regions_);

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");

  return out_list;
}


/* ******************************************************************
* Filters out empty sublists starting with node "parent".
****************************************************************** */
void InputConverterU::MergeInitialConditionsLists_(Teuchos::ParameterList& plist)
{
  if (plist.sublist("PKs").isSublist("chemistry")) {
    Teuchos::ParameterList& ics = plist.sublist("state")
                                       .sublist("initial conditions");
    Teuchos::ParameterList& icc = plist.sublist("PKs").sublist("chemistry")
                                       .sublist("initial conditions");

    for (auto it = icc.begin(); it != icc.end(); ++it) {
      if (icc.isSublist(it->first)) {
        Teuchos::ParameterList& slist = icc.sublist(it->first);
        if (slist.isSublist("function")) {
          ics.sublist(it->first) = slist;
          slist.set<std::string>("function", "list was moved to state");
        }
      }
    }
  }  
}


/* ******************************************************************
* Filters out empty sublists starting with node "parent".
****************************************************************** */
void InputConverterU::FilterEmptySublists_(Teuchos::ParameterList& plist)
{
  for (auto it = plist.begin(); it != plist.end(); ++it) {
    if (plist.isSublist(it->first)) {
      Teuchos::ParameterList& slist = plist.sublist(it->first);
      if (slist.numParams() == 0) {
        plist.remove(it->first);
      } else {
        FilterEmptySublists_(slist);
      }
    }
  }
}


/* ******************************************************************
* Search tools
****************************************************************** */
bool InputConverterU::FindNameInVector_(
   const std::string& name, const std::vector<std::string>& list)
{
  for (int i = 0; i < list.size(); ++i) {
    if (name == list[i]) return true;
  }
  return false;
}


/* ******************************************************************
* Create automatic string name.
****************************************************************** */
std::string InputConverterU::CreateNameFromVector_(const std::vector<std::string>& list)
{
  std::string str;
  for (auto it = list.begin(); it != list.end(); ++it) {
    str = str + *it;
  }
  return str;
}


/* ******************************************************************
* Output of XML
****************************************************************** */
void InputConverterU::SaveXMLFile(
    Teuchos::ParameterList& out_list, std::string& xmlfilename)
{
  bool flag;
  int precision(0);
  std::string filename("");

  DOMNode* node = GetUniqueElementByTagsString_("misc, echo_translated_input", flag);
  if (flag) {
    DOMElement* element = static_cast<DOMElement*>(node);
    filename = GetAttributeValueS_(element, "file_name", TYPE_NONE, false, "skip");
    precision = GetAttributeValueL_(element, "output_precision", TYPE_NONE, 0, 10, false, 0);
  }

  if (filename == "") {
    filename = xmlfilename;
    std::string new_extension("_native_v8.xml");
    size_t pos = filename.find(".xml");
    filename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)14);
  }

  if (filename != "skip") {
    if (vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "Writing the translated XML to " << filename.c_str() << std::endl;;
    }

    Teuchos::Amanzi_XMLParameterListWriter XMLWriter;
    if (precision > 0) XMLWriter.set_precision(precision);
    Teuchos::XMLObject XMLobj = XMLWriter.toXML(out_list);

    std::ofstream xmlfile;
    xmlfile.open(filename.c_str());
    xmlfile << XMLobj;
  }
}


/* ******************************************************************
* Create a name by concatenating names in a list
****************************************************************** */
std::string InputConverterU::CreateUniqueName_(const Teuchos::Array<std::string>& list)
{
  std::string name;
  for (auto it = list.begin(); it != list.end(); ++it) {
    name.append(*it);
  }
  return name;
}


/* ******************************************************************
* Print final comments
****************************************************************** */
void InputConverterU::PrintStatistics_()
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Final comments:\n found units: ";
    int n(0);
    for (std::set<std::string>::iterator it = found_units_.begin(); it != found_units_.end(); ++it) {
      *vo_->os() << *it << " ";
      if (++n > 10) {
        n = 0;
        *vo_->os() << "\n continue:    ";
      }
    }
    *vo_->os() << std::endl;
  }
}

}  // namespace AmanziInput
}  // namespace Amanzi
