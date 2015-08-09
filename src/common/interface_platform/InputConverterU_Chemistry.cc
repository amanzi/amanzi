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

  if (minerals.size() > 0) {
    out_list.set<Teuchos::Array<std::string> >("Minerals", minerals);
  }

  // region specific initial conditions
  Teuchos::ParameterList& ic_list = out_list.sublist("initial conditions");

/*
    Teuchos::ParameterList& matprop_list = plist->sublist("Material Properties");
    DOMNode* inode = children->item(i);

    for (Teuchos::ParameterList::ConstIterator i = matprop_list.begin(); i != matprop_list.end(); i++) {
      // get the regions
      Teuchos::Array<std::string> regions = matprop_list.sublist(matprop_list.name(i)).get<Teuchos::Array<std::string> >("Assigned Regions");

      if (minerals.size() > 0) {
        double mvf(0.0), msa(0.0);

        // mineral volume fractions
        Teuchos::ParameterList &mineral_volfrac_ic = ic_list.sublist("mineral_volume_fractions");
        // mineral specific surface area
        Teuchos::ParameterList &mineral_surfarea_ic = ic_list.sublist("mineral_specific_surface_area");

        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {

          // mineral volume fractions and specific surface area

          Teuchos::ParameterList& aux1_list =
              mineral_volfrac_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          Teuchos::ParameterList& aux2_list =
              mineral_surfarea_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          aux1_list.set<int>("Number of DoFs", minerals.size())
              .set("Function type", "composite function");

          aux2_list.set<int>("Number of DoFs", minerals.size())
              .set("Function type", "composite function");

          for (int j = 0; j<minerals.size(); ++j) {
            std::stringstream ss;
            ss << "DoF " << j + 1 << " Function";

            mvf = matprop_list.sublist(matprop_list.name(i)).sublist("Mineralogy").sublist(minerals[j])
                .get<double>("Volume Fraction");

            msa = matprop_list.sublist(matprop_list.name(i)).sublist("Mineralogy").sublist(minerals[j])
                .get<double>("Specific Surface Area");

            aux1_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", mvf);

            aux2_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", msa);
          }
        }
      }


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

}  // namespace AmanziInput
}  // namespace Amanzi


