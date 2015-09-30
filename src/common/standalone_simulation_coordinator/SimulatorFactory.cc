#include "errors.hh"
#include "exceptions.hh"

#ifdef ENABLE_Unstructured
#include "AmanziUnstructuredGridSimulationDriver.hh"
#endif
#ifdef ENABLE_Structured
#include "amanzi_structured_grid_simulation_driver.H"
#endif

#include "InputConverter.hh"
#include "SimulatorFactory.hh"

using namespace std;
using namespace xercesc;

namespace Amanzi
{

namespace SimulatorFactory
{

Simulator* Create(const std::string& input_filename)
{
  XercesDOMParser* parser = Amanzi::AmanziInput::CreateXMLParser();
  DOMDocument* doc = Amanzi::AmanziInput::OpenXMLInput(parser, input_filename);

  // Determine whether this is a structured or unstructured input and dispatch 
  // accordingly.
  char str[100];
  XMLCh xstr[100];
  DOMElement* element = doc->getDocumentElement();
  XMLString::transcode(element->getTagName(), str, 99);
  string version;
  if (strcmp(str, "amanzi_input") == 0)
    version = "v2";
  else
  {
    if (strcmp(str, "ParameterList") == 0)
      version = "v1";
    else
      Exceptions::amanzi_throw(Errors::Message("Invalid input file."));
  }

  // Figure out the type of input (structured, unstructured).
  string type;
  if (version == "v2")
  {
    XMLString::transcode("type", xstr, 99);
    if (not element->hasAttribute(xstr)) 
      Exceptions::amanzi_throw(Errors::Message("Invalid input file (no 'type' attribute in amanzi_input)."));
    XMLString::transcode(element->getAttribute(xstr), str, 99);
    type = str;
  }
  else
  {
    // We check the parameter lists for the Structured tag. If we don't 
    // find this, we assume it to be unstructured.
    type = "unstructured";
    XMLString::transcode("ParameterList", xstr, 99);
    DOMNodeList* nodes = doc->getElementsByTagName(xstr);
    for (int i = 0; i < nodes->getLength(); ++i)
    {
      DOMElement* element = static_cast<DOMElement*>(nodes->item(i));
      if (element != NULL)
      {
        XMLString::transcode("name", xstr, 99);
        XMLString::transcode(element->getAttribute(xstr), str, 99);
        if (strcmp(str, "Structured") == 0)
        {
          type = "structured";
          break;
        }
      }
    }
  }

  // Create the appropriate simulator.
  Simulator* simulator = NULL;
  if (type == "structured")
  {
#ifdef ENABLE_Structured
    // Uncomment the following lines when the new v2 -> PP translator works.
//    if (version == "v2")
//      simulator = new AmanziStructuredGridSimulationDriver(doc);
//    else 
      simulator = new AmanziStructuredGridSimulationDriver(input_filename);
#else
    amanzi_throw(Errors::Message("Structured not supported in current build"));
#endif
  }
  else if (type == "unstructured")
  {
#ifdef ENABLE_Unstructured
    if (version == "v2")
      simulator = new AmanziUnstructuredGridSimulationDriver(input_filename, doc);
    else 
      simulator = new AmanziUnstructuredGridSimulationDriver(input_filename);
#else
    amanzi_throw(Errors::Message("Unstructured not supported in current build"));
#endif
  }

  // Clean up.
  delete parser;
  XMLPlatformUtils::Terminate();

  return simulator;
}

} // end namespace SimulatorFactory
} // end namespace amanzi

