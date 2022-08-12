/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Sergi Molins <smolins@lbl.gov>
*/

#include <algorithm>
#include <sstream>
#include <string>
#include <sys/stat.h>

//TPLs
#include "Teuchos_ParameterList.hpp"
#include "xercesc/dom/DOM.hpp"

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputConverter.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create thermodynamic database
****************************************************************** */
Teuchos::ParameterList InputConverter::TranslateThermodynamicDatabase_() 
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;
  DOMNode *node, *inode, *knode;
  DOMNodeList *children;
  DOMElement* element;

  Errors::Message msg("error in thermodynamic database format");

  bool flag;
  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, primary_species", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("primary"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);
      std::string name = GetAttributeValueS_(inode, "name");
      
      knode = GetUniqueElementByTagsString_(inode, "species_data", flag);
      if (!flag) Exceptions::amanzi_throw(msg);

      double ion_size = GetAttributeValueD_(knode, "ion_size", TYPE_NUMERICAL, 0.0, DVAL_MAX);
      int z = GetAttributeValueL_(knode, "z", TYPE_NUMERICAL, -100, 100);
      double mol_weight = GetAttributeValueD_(knode, "weight", TYPE_NUMERICAL, 0.0, 10000.0);

      out_list.sublist("primary species").sublist(name)
          .set<double>("ion size parameter", ion_size)
          .set<int>("charge", z)
          .set<double>("gram molecular weight", mol_weight);
    }
  }

  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, aqueous_equilibrium_complexes", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("complex"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);
      std::string name = GetAttributeValueS_(inode, "name");

      knode = GetUniqueElementByTagsString_(inode, "species_data", flag);
      if (!flag) Exceptions::amanzi_throw(msg);

      double ion_size = GetAttributeValueD_(knode, "ion_size", TYPE_NUMERICAL, 0.0, DVAL_MAX);
      int z = GetAttributeValueL_(knode, "z", TYPE_NUMERICAL, -100, 100);
      double mol_weight = GetAttributeValueD_(knode, "weight", TYPE_NUMERICAL, 0.0, 10000.0);

      knode = GetUniqueElementByTagsString_(inode, "reaction", flag);
      if (!flag) Exceptions::amanzi_throw(msg);
      std::string reaction = TrimString_(mm.transcode(knode->getTextContent()));

      knode = GetUniqueElementByTagsString_(inode, "equilibrium_constant", flag);
      if (!flag) Exceptions::amanzi_throw(msg);
      double lnKeq = strtod(mm.transcode(knode->getTextContent()), NULL);

      out_list.sublist("aqueous equilibrium complexes").sublist(name)
          .set<double>("ion size parameter", ion_size)
          .set<int>("charge", z)
          .set<double>("gram molecular weight", mol_weight)
          .set<std::string>("reaction", reaction)
          .set<double>("equilibrium constant", lnKeq);
    }
  }

  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, mineral_kinetics", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("mineral"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);
      std::string name = GetAttributeValueS_(inode, "name");

      knode = GetUniqueElementByTagsString_(inode, "species_data", flag);
      if (!flag) Exceptions::amanzi_throw(msg);
      double mol_weight = GetAttributeValueD_(knode, "weight", TYPE_NUMERICAL, 0.0, 1000.0);

      knode = GetUniqueElementByTagsString_(inode, "reaction", flag);
      if (!flag) Exceptions::amanzi_throw(msg);
      std::string reaction = TrimString_(mm.transcode(knode->getTextContent()));

      knode = GetUniqueElementByTagsString_(inode, "equilibrium_constant", flag);
      if (!flag) Exceptions::amanzi_throw(msg);
      double lnKeq = strtod(mm.transcode(knode->getTextContent()), NULL);

      knode = GetUniqueElementByTagsString_(inode, "kinetics_data", flag);
      if (!flag) Exceptions::amanzi_throw(msg);

      std::string model = GetAttributeValueS_(knode, "model", "TST");
      double rate = GetAttributeValueD_(knode, "rate", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);
      std::string modifiers = GetAttributeValueS_(knode, "modifiers");
      double mol_volume = GetAttributeValueD_(knode, "molar_volume", TYPE_NUMERICAL, 0.0, DVAL_MAX);
      double ssa = GetAttributeValueD_(knode, "specific_surface_area", TYPE_NUMERICAL, 0.0, DVAL_MAX);

      out_list.sublist("mineral kinetics").sublist(name)
          .set<std::string>("rate model", model)
          .set<double>("rate constant", rate)
          .set<std::string>("modifiers", modifiers)
          .set<double>("gram molecular weight", mol_weight)
          .set<std::string>("reaction", reaction)
          .set<double>("equilibrium constant", lnKeq)
          .set<double>("molar volume", mol_volume)
          .set<double>("specific surface area", ssa);
    }
  }

  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, general_kinetics", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("general"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);

      knode = GetUniqueElementByTagsString_(inode, "reaction", flag);
      if (!flag) Exceptions::amanzi_throw(msg);

      std::string reactants = GetAttributeValueS_(knode, "reactants");
      std::string products = GetAttributeValueS_(knode, "products");

      knode = GetUniqueElementByTagsString_(inode, "kinetics_data", flag);
      if (!flag) Exceptions::amanzi_throw(msg);

      double frate = GetAttributeValueD_(knode, "forward_rate", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);
      double brate = GetAttributeValueD_(knode, "backward_rate", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);

      out_list.sublist("general kinetics").sublist("general_" + std::to_string(i)) 
          .set<std::string>("reactants", reactants)
          .set<std::string>("products", products)
          .set<double>("forward rate", frate)
          .set<double>("backward rate", brate);
    }
  }

  // ion exchange sites and ion complexes
  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, ion_exchange_sites", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("site"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);

      std::string name = GetAttributeValueS_(inode, "name");
      std::string location = GetAttributeValueS_(inode, "location");
      int z = GetAttributeValueL_(inode, "z", TYPE_NUMERICAL, -100, 100);

      out_list.sublist("ion exchange sites").sublist(name) 
          .set<std::string>("location", location)
          .set<int>("charge", z);
    }
  }

  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, ion_exchange_complexes", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("complex"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);
      std::string name = GetAttributeValueS_(inode, "name");

      knode = GetUniqueElementByTagsString_(inode, "reaction", flag);
      if (!flag) Exceptions::amanzi_throw(msg);
      std::string reaction = TrimString_(mm.transcode(knode->getTextContent()));

      knode = GetUniqueElementByTagsString_(inode, "equilibrium_constant", flag);
      if (!flag) Exceptions::amanzi_throw(msg);
      double lnKeq = strtod(mm.transcode(knode->getTextContent()), NULL);

      out_list.sublist("ion exchange complexes").sublist(name) 
          .set<std::string>("reaction", reaction)
          .set<double>("equilibrium constant", lnKeq);
    }
  }

  // surface sites and surface complexes
  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, surface_complex_sites", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("site"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);

      std::string name = GetAttributeValueS_(inode, "name");
      double density = GetAttributeValueD_(inode, "density", TYPE_NUMERICAL, 0.0, DVAL_MAX);

      out_list.sublist("surface complex sites").sublist(name) 
          .set<double>("density", density);
    }
  }

  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, surface_complexes", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("complex"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);

      std::string name = GetAttributeValueS_(inode, "name");
      int z = GetAttributeValueL_(inode, "z", TYPE_NUMERICAL, -100, 100);

      knode = GetUniqueElementByTagsString_(inode, "reaction", flag);
      if (!flag) Exceptions::amanzi_throw(msg);
      std::string reaction = TrimString_(mm.transcode(knode->getTextContent()));

      knode = GetUniqueElementByTagsString_(inode, "equilibrium_constant", flag);
      if (!flag) Exceptions::amanzi_throw(msg);
      double lnKeq = strtod(mm.transcode(knode->getTextContent()), NULL);

      out_list.sublist("surface complexes").sublist(name) 
          .set<int>("charge", z)
          .set<std::string>("reaction", reaction)
          .set<double>("equilibrium constant", lnKeq);
    }
  }

  // isotherms
  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, isotherms", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("isotherm"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);

      std::string name = GetAttributeValueS_(inode, "species");
      std::string model = GetAttributeValueS_(inode, "model", "linear, langmuir, freundlich");
      std::vector<double> params = GetAttributeVectorD_(inode, "parameters");

      out_list.sublist("isotherms").sublist(name) 
          .set<std::string>("model", model)
          .set<Teuchos::Array<double> >("parameters", params);
    }
  }

  // radioactive decay
  node = GetUniqueElementByTagsString_("geochemistry, thermodynamic_database, radioactive_decay", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    children = element->getElementsByTagName(mm.transcode("reaction"));

    for (int i = 0; i < children->getLength(); ++i) {
      inode = children->item(i);

      std::string reactant = GetAttributeValueS_(inode, "reactant");
      std::string product = GetAttributeValueS_(inode, "product");
      double half_life = GetAttributeValueD_(inode, "half_life", TYPE_NUMERICAL, 0.0, DVAL_MAX);

      out_list.sublist("radioactive decay").sublist("complex_" + std::to_string(i))
          .set<std::string>("reactant", reactant)
          .set<std::string>("product", product)
          .set<double>("half life", half_life);
    }
  }

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi


