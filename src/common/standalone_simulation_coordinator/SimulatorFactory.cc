/*
  Simulator 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatskiy (original version)
*/

#include "xercesc/util/PlatformUtils.hpp"

#include "errors.hh"
#include "exceptions.hh"
#include "InputConverter.hh"
#include "SimulatorFactory.hh"

#ifdef ENABLE_Unstructured
#include "AmanziUnstructuredGridSimulationDriver.hh"
#endif
#ifdef ENABLE_Structured
#include "amanzi_structured_grid_simulation_driver.H"
#endif

using namespace std;
XERCES_CPP_NAMESPACE_USE

namespace Amanzi {
namespace SimulatorFactory {

std::unique_ptr<Simulator>
Create(const std::string& input_filename, const std::string& output_prefix)
{
  XercesDOMParser* parser = Amanzi::AmanziInput::CreateXMLParser();
  DOMDocument* doc = Amanzi::AmanziInput::OpenXMLInput(parser, input_filename);

  // Determine whether this is a structured or unstructured input and dispatch 
  // accordingly.
  char str[100];
  XMLCh xstr[100];
  DOMElement* element = doc->getDocumentElement();
  XMLString::transcode(element->getTagName(), str, 99);
  std::string version;
  if (strcmp(str, "amanzi_input") == 0) {
    version = "v2";
  } else {
    if (strcmp(str, "ParameterList") == 0)
      version = "v1";
    else
      Exceptions::amanzi_throw(Errors::Message("Invalid input file."));
  }

  // Figure out the type of input (structured, unstructured).
  std::string type;
  if (version == "v2") {
    XMLString::transcode("type", xstr, 99);
    if (not element->hasAttribute(xstr)) 
      Exceptions::amanzi_throw(Errors::Message("Invalid input file (no 'type' attribute in amanzi_input)."));
    XMLString::transcode(element->getAttribute(xstr), str, 99);
    type = str;
  } else {
    // We check the parameter lists for the unstructured tag. If we don't 
    // find this, we assume it to be structured.
    type = "structured";
    XMLString::transcode("name", xstr, 99);

    DOMNodeList* children = doc->getDocumentElement()->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* child = children->item(i);
      if (child->getNodeType() != DOMNode::ELEMENT_NODE) continue;
      element = static_cast<DOMElement*>(child);
      XMLString::transcode(element->getAttribute(xstr), str, 99);
      if (strcmp(str, "Native Unstructured Input") == 0) {
        type = "unstructured";
        break;
      }
    }
  }

  // Create the appropriate simulator.
  std::unique_ptr<Simulator> simulator = nullptr;
  if (type == "structured") {
#ifdef ENABLE_Structured
    if (version == "v2")
      simulator = std::make_unique<AmanziStructuredGridSimulationDriver>(input_filename, doc, output_prefix);
    else
      amanzi_throw(Errors::Message("Input spec v1 is no longer supported by Amanzi-S."));
#else
    amanzi_throw(Errors::Message("Structured not supported in current build"));
#endif
  }
  else if (type == "unstructured") {
#ifdef ENABLE_Unstructured
    if (version == "v2")
      simulator = std::make_unique<AmanziUnstructuredGridSimulationDriver>(input_filename, doc, output_prefix);
    else
      simulator = std::make_unique<AmanziUnstructuredGridSimulationDriver>(input_filename);
#else
    amanzi_throw(Errors::Message("Unstructured not supported in current build"));
#endif
  }

  // Clean up.
  delete parser;
  XMLPlatformUtils::Terminate();

  return simulator;
}

} // namespace SimulatorFactory
} // namespace Amanzi

