/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Input Converter

*/

#include <algorithm>
#include <sstream>
#include <string>

//TPLs
#include <xercesc/dom/DOM.hpp>

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"
#include "StringExt.hh"
#include "Units.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create list of porosity models.
****************************************************************** */
void
InputConverterU::TranslatePOM_(const std::string& domain,
                               Teuchos::ParameterList& out_ic,
                               Teuchos::ParameterList& out_ev)
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating porosity models" << std::endl;

  MemoryManager mm;
  Errors::Message msg;

  DOMNodeList* children;
  DOMNode *node, *knode;
  DOMElement* element;

  bool flag, poroelasticity(false), thermoelasticity(false);
  double compres, ref_pressure;

  std::string model;
  Key porosity_key = Keys::getKey(domain, "porosity");

  use_porosity_model_ = false;

  node = (domain == "fracture") ?
           GetUniqueElementByTagsString_("fracture_network, materials", flag) :
           GetUniqueElementByTagsString_("materials", flag);
  element = static_cast<DOMElement*>(node);
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);

    // get assigned regions
    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));
    std::string reg_str = CreateNameFromVector_(regions);

    // get compressibility
    node = GetUniqueElementByTagsString_(inode, "mechanical_properties, porosity", flag);
    if (flag) {
      TranslateFieldEvaluator_(node, porosity_key, "-", reg_str, regions, out_ic, out_ev);
    } else {
      msg << "Porosity element must be specified under mechanical_properties";
      Exceptions::amanzi_throw(msg);
    }

    model =
      GetAttributeValueS_(node, "model", "constant, compressible, poroelastic, thermoporoelastic");

    std::string type = GetAttributeValueS_(node, "type", TYPE_NONE, false, "");
    if (type == "h5file") {
      use_porosity_model_ = false;
      break;
    }

    double phi = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, 1.0);
    ref_pressure = GetAttributeValueD_(
      node, "reference_pressure", TYPE_NUMERICAL, 0.0, DBL_MAX, "Pa", false, const_atm_pressure_);

    double dilation_rock(0.0), dilation_liquid(0.0), biot(1.0), bulk(0.0);

    // optional poroelasticity, Biot-Willis coefficient
    if (model == "thermoporoelastic" || model == "poroelastic") {
      poroelasticity = true;
      use_porosity_model_ = true;
      knode = GetUniqueElementByTagsString_(inode, "mechanical_properties, biot_coefficient", flag);
      if (flag) biot = GetAttributeValueD_(knode, "value", TYPE_NUMERICAL, 0.0, 1.0, "");
    }

    // optional poroelasticity
    if (model == "thermoporoelastic" || model == "poroelastic" || model == "compressible") {
      use_porosity_model_ = true;
      compres = GetAttributeValueD_(node, "compressibility", TYPE_NUMERICAL, 0.0, 1.0, "Pa^-1");
    }

    // optional thermoporoelasticity
    if (model == "thermoporoelastic") {
      poroelasticity = true;
      thermoelasticity = true;
      use_porosity_model_ = true;

      knode =
        GetUniqueElementByTagsString_(inode, "mechanical_properties, rock_thermal_dilation", flag);
      if (flag)
        dilation_rock = GetAttributeValueD_(knode, "value", TYPE_NUMERICAL, 0.0, 1.0, "K^-1");

      knode = GetUniqueElementByTagsString_(
        inode, "mechanical_properties, liquid_thermal_dilation", flag);
      if (flag)
        dilation_liquid = GetAttributeValueD_(knode, "value", TYPE_NUMERICAL, 0.0, 1.0, "K^-1");

      if (flag)
        bulk = GetAttributeValueD_(node, "solid_bulk_modulus", TYPE_NUMERICAL, 0.0, DVAL_MAX, "Pa");
      compres = (biot - phi) / bulk;
    }

    std::stringstream ss;
    ss << "POM " << i;

    Teuchos::ParameterList& pom_list = out_list.sublist(ss.str());
    pom_list.set<Teuchos::Array<std::string>>("regions", regions);

    // we can have either uniform of compressible rock
    pom_list.set<std::string>("porosity model", model);
    if (model == "constant") {
      pom_list.set<double>("value", phi);
    } else {
      pom_list.set<double>("undeformed soil porosity", phi)
        .set<double>("reference pressure", ref_pressure)
        .set<double>("pore compressibility", compres)
        .set<double>("biot coefficient", biot)
        .set<double>("rock thermal dilation", dilation_rock)
        .set<double>("liquid thermal dilation", dilation_liquid);
      use_porosity_model_ = true;
    }
  }

  if (use_porosity_model_) {
    out_ev.sublist(porosity_key).sublist("parameters") = out_list;

    Key pressure_key = Keys::getKey(domain, "pressure");
    Key temperature_key = Keys::getKey(domain, "temperature");
    Key vol_strain_key = Keys::getKey(domain, "volumetric_strain");

    out_ev.sublist(porosity_key)
      .set<std::string>("evaluator type", "porosity")
      .set<std::string>("porosity key", porosity_key)
      .set<std::string>("pressure key", pressure_key)
      .set<std::string>("temperature key", temperature_key)
      .set<std::string>("volumetric strain key", vol_strain_key)
      .set<bool>("thermoelasticity", thermoelasticity)
      .set<bool>("poroelasticity", poroelasticity)
      .set<std::string>("tag", "");
  }
}

} // namespace AmanziInput
} // namespace Amanzi
