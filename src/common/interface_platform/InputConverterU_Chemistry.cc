/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <sstream>
#include <string>

//TPLs
#include "Teuchos_ParameterList.hpp"
#include "xercesc/dom/DOM.hpp"

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create flow list.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateChemistry_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating chemistry" << std::endl;

  XString mm;
  char* text;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;

  bool flag;
  node = getUniqueElementByTagNames_("process_kernels", "chemistry", flag);
  std::string engine = GetAttributeValueS_(static_cast<DOMElement*>(node), "engine");

  // process engine
  bool native(false);
  if (strcmp(engine.c_str(), "amanzi") == 0) {
    out_list.set<std::string>("chemistry model", "Amanzi");
    native = true;
  } else if (strcmp(engine.c_str(), "pflotran") == 0) {
    out_list.set<std::string>("chemistry model", "Alquimia");
    out_list.set<std::string>("Engine", "PFloTran");
  }
  
  // general parameters
  out_list.set<int>("Number of component concentrations", comp_names_all_.size());

  // minerals
  node = getUniqueElementByTagsString_("phases, solid_phase, minerals", flag);
  children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode("mineral"));

  std::vector<std::string> minerals;
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    text = mm.transcode(inode->getTextContent());
    minerals.push_back(text);
  }

  // region specific initial conditions
  Teuchos::ParameterList& ic_list = out_list.sublist("initial conditions");

  node_list = doc_->getElementsByTagName(mm.transcode("materials"));
  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (minerals.size() > 0) {
      out_list.set<Teuchos::Array<std::string> >("Minerals", minerals);

      node = getUniqueElementByTagNames_(inode, "assigned_regions", flag);
      text = mm.transcode(node->getTextContent());
      std::vector<std::string> regions = CharToStrings_(text);

      Teuchos::ParameterList &volfrac = ic_list.sublist("mineral_volume_fractions");
      Teuchos::ParameterList &surfarea = ic_list.sublist("mineral_specific_surface_area");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); it++) {
        Teuchos::ParameterList& aux1_list = volfrac.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux1_list.set<int>("Number of DoFs", minerals.size())
            .set("Function type", "composite function");

        Teuchos::ParameterList& aux2_list = surfarea.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux2_list.set<int>("Number of DoFs", minerals.size())
            .set("Function type", "composite function");

        for (int j = 0; j < minerals.size(); ++j) {
          std::stringstream ss;
          ss << "DoF " << j + 1 << " Function";

          node = getUniqueElementByTagNames_(inode, "minerals", "mineral", flag);  // FIX ME
          if (flag) {
            element = static_cast<DOMElement*>(node);
            double mvf = GetAttributeValueD_(element, "volume_fraction", false, 0.0);
            double msa = GetAttributeValueD_(element, "specific_surface_area", false, 0.0);

            aux1_list.sublist(ss.str()).sublist("function-constant").set<double>("value", mvf);
            aux2_list.sublist(ss.str()).sublist("function-constant").set<double>("value", msa);
          }
        }
      }
    }
  }

/*
      if (matprop_list.sublist(matprop_list.name(i)).isParameter("Cation Exchange Capacity")) {
        double cec = matprop_list.sublist(matprop_list.name(i)).get<double>("Cation Exchange Capacity");

        Teuchos::ParameterList &ion_exchange_sites_ic = ic_list.sublist("ion_exchange_sites");
        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {
          Teuchos::ParameterList& aux1_list =
              ion_exchange_sites_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          // this needs to be more general... for now we initialize with one DoF
          aux1_list.set<int>("Number of DoFs", 1)
              .set("Function type", "composite function");

          aux1_list.sublist("DoF 1 Function").sublist("function-constant")
              .set<double>("value", cec);
        }

        Teuchos::ParameterList &ion_exchange_ref_cation_conc_ic = ic_list.sublist("ion_exchange_ref_cation_conc");
        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {
          Teuchos::ParameterList& aux1_list =
              ion_exchange_ref_cation_conc_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          // this needs to be more general... for now we initialize with one DoF
          aux1_list.set<int>("Number of DoFs", 1)
              .set("Function type", "composite function");

          aux1_list.sublist("DoF 1 Function").sublist("function-constant")
              .set<double>("value", 1.0);  // this should be read from the input file??... TODO
        }
      }

      if (matprop_list.sublist(matprop_list.name(i)).isSublist("Sorption Isotherms")) {
        Teuchos::ParameterList& isotherm_kd_ic = ic_list.sublist("isotherm_kd");
        Teuchos::ParameterList& isotherm_langmuir_b_ic = ic_list.sublist("isotherm_langmuir_b");
        Teuchos::ParameterList& isotherm_freundlich_n_ic = ic_list.sublist("isotherm_freundlich_n");

        Teuchos::ParameterList& sorption_isotherms_list = matprop_list.sublist(matprop_list.name(i)).sublist("Sorption Isotherms");

        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {
          // Kd
          {
            Teuchos::ParameterList& aux1_list =
                isotherm_kd_ic.sublist("function").sublist(*ir)
                .set<std::string>("region",*ir)
                .set<std::string>("component","cell")
                .sublist("function");

            aux1_list.set<int>("Number of DoFs", comp_names_.size())
                .set("Function type", "composite function");

            for ( int ic = 0; ic != comp_names_.size(); ++ic) {

              std::stringstream ss;
              ss << "DoF " << ic + 1 << " Function";

              double kd(0.0);
              if (sorption_isotherms_list.isSublist(comp_names_[ic])) {
                kd = sorption_isotherms_list.sublist(comp_names_[ic]).get<double>("Kd",0.0);
              }

              aux1_list.sublist(ss.str()).sublist("function-constant")
                  .set<double>("value", kd);
            }
          }

          // Langmuir
          {
            Teuchos::ParameterList& aux2_list =
                isotherm_langmuir_b_ic.sublist("function").sublist(*ir)
                .set<std::string>("region",*ir)
                .set<std::string>("component","cell")
                .sublist("function");

            aux2_list.set<int>("Number of DoFs", comp_names_.size())
                .set("Function type", "composite function");

            for ( int ic = 0; ic != comp_names_.size(); ++ic) {

              std::stringstream ss;
              ss << "DoF " << ic + 1 << " Function";

              double langmuir_b(1.0);
              if (sorption_isotherms_list.isSublist(comp_names_[ic])) {
                langmuir_b = sorption_isotherms_list.sublist(comp_names_[ic]).get<double>("Langmuir b",1.0);
              }

              aux2_list.sublist(ss.str()).sublist("function-constant")
                  .set<double>("value", langmuir_b);
            }
          }

          // Freundlich
          {
            Teuchos::ParameterList& aux3_list =
                isotherm_freundlich_n_ic.sublist("function").sublist(*ir)
                .set<std::string>("region",*ir)
                .set<std::string>("component","cell")
                .sublist("function");

            aux3_list.set<int>("Number of DoFs", comp_names_.size())
                .set("Function type", "composite function");

            for ( int ic = 0; ic != comp_names_.size(); ++ic) {

              std::stringstream ss;
              ss << "DoF " << ic + 1 << " Function";

              double freundlich_n(1.0);
              if (sorption_isotherms_list.isSublist(comp_names_[ic])) {
                freundlich_n = sorption_isotherms_list.sublist(comp_names_[ic]).get<double>("Freundlich n",1.0);
              }

              aux3_list.sublist(ss.str()).sublist("function-constant")
                  .set<double>("value", freundlich_n);
            }
          }
        }
      }

      if (sorption_site_names_.size() > 0) {
        Teuchos::ParameterList& sorption_sites_ic = ic_list.sublist("sorption_sites");

        Teuchos::ParameterList& sorption_sites_list = matprop_list.sublist(matprop_list.name(i)).sublist("Surface Complexation Sites");

        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {

          Teuchos::ParameterList& aux1_list =
              sorption_sites_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          aux1_list.set<int>("Number of DoFs", sorption_site_names_.size())
              .set("Function type", "composite function");

          for ( int ic = 0; ic != sorption_site_names_.size(); ++ic) {

            std::stringstream ss;
            ss << "DoF " << ic + 1 << " Function";

            double value(0.0);
            if (sorption_sites_list.isSublist(sorption_site_names_[ic])) {
              value = sorption_sites_list.sublist(sorption_site_names_[ic]).get<double>("Site Density",0.0);
            }

            aux1_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", value);
          }
        }
      }
    }

    Teuchos::ParameterList& ic_list = plist->sublist("Initial Conditions");

    for (Teuchos::ParameterList::ConstIterator ic = ic_list.begin(); ic != ic_list.end(); ++ic) {
      if (ic_list.isSublist(ic->first)) {
        Teuchos::ParameterList& ics = ic_list.sublist(ic->first);

        Teuchos::Array<std::string> ass_regions = ics.get<Teuchos::Array<std::string> >("Assigned Regions");

        Teuchos::ParameterList &free_ion_species_ic = ic_list.sublist("free_ion_species");

        for (Teuchos::Array<std::string>::const_iterator ir=ass_regions.begin(); ir!=ass_regions.end(); ir++) {
          Teuchos::ParameterList& aux1_list =
              free_ion_species_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          aux1_list.set<int>("Number of DoFs", comp_names_.size())
              .set("Function type", "composite function");

          for (int j = 0; j<comp_names_.size(); ++j) {
            std::stringstream ss;
            ss << "DoF " << j + 1 << " Function";

            double value(1.0e-9);
            value = ics.sublist("Solute IC").sublist(phases_[0].name)
                .sublist(phases_[0].solute_name).sublist(comp_names_[j])
                .sublist("IC: Uniform Concentration").get<double>("Free Ion Guess",1.0e-9);

            aux1_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", value);
          }
        }
      }
    }
  }
*/

  out_list.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");
  return out_list;
}


/* ******************************************************************
* Empty
****************************************************************** */
void InputConverterU::CreateBDGFile(Teuchos::ParameterList sorption_list, Teuchos::ParameterList def_list)
{

  std::ofstream bgd_file;
  std::stringstream species_string;
  std::stringstream isotherms_string;

  // build streams
  for (Teuchos::ParameterList::ConstIterator i = sorption_list.begin(); i != sorption_list.end(); i++) {
    Teuchos::ParameterList& tmpList = sorption_list.sublist(sorption_list.name(i)) ;
    species_string << sorption_list.name(i) << " ;   0.00 ;   0.00 ;   1.00 \n";
    if (tmpList.isParameter("Langmuir b")) {
      isotherms_string << sorption_list.name(i) << " ; langmuir ; " << tmpList.get<double>("Kd")<< " " <<tmpList.get<double>("Langmuir b") << std::endl;
    } else if (tmpList.isParameter("Freundlich n")) {
      isotherms_string << sorption_list.name(i) << " ; freundlich ; " << tmpList.get<double>("Kd")<< " " <<tmpList.get<double>("Freundlich n") << std::endl;
    } else {
      isotherms_string << sorption_list.name(i) << " ; linear ; " << tmpList.get<double>("Kd")<< std::endl;
    }
  }
  
  // build bgd filename
  std::string bgdfilename;
  if (def_list.isParameter("xmlfilename") ) {
    bgdfilename = def_list.get<std::string>("xmlfilename");
    std::string new_extension(".bgd");
    size_t pos = bgdfilename.find(".xml");
    bgdfilename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)4);
  } else {
    bgdfilename = "isotherms.bgd";
  }

  // open output bgd file
  bgd_file.open(bgdfilename.c_str());

  // <Primary Species
  bgd_file << "<Primary Species\n";
  bgd_file << species_string.str();

  //<Isotherms
  bgd_file << "<Isotherms\n" ;
  bgd_file << "# Note, these values will be overwritten by the xml file\n" ;
  bgd_file << isotherms_string.str();

  // close output bdg file
  bgd_file.close();
}

}  // namespace AmanziInput
}  // namespace Amanzi


