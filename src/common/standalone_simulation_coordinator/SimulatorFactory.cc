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
  DOMDocument* doc = Amanzi::AmanziInput::OpenXMLInput(input_filename);

  // Determine whether this is a structured or unstructured input and dispatch 
  // accordingly.
  char str[100];
  XMLCh xstr[100];
  DOMElement* element = doc->getDocumentElement();
  XMLString::transcode(element->getTagName(), str, 99);
  if (strcmp(str, "amanzi_input") != 0)
    Exceptions::amanzi_throw(Errors::Message("Invalid input file."));

  // Figure out the version of the input spec used.
  XMLString::transcode("version", xstr, 99);
  if (not element->hasAttribute(xstr)) 
    Exceptions::amanzi_throw(Errors::Message("Invalid input file (no 'version' attribute in amanzi_input)."));
  XMLString::transcode(element->getAttribute(xstr), str, 99);
  string version = str;

  // Figure out the type of input (structured, unstructured).
  XMLString::transcode("type", xstr, 99);
  if (not element->hasAttribute(xstr)) 
    Exceptions::amanzi_throw(Errors::Message("Invalid input file (no 'type' attribute in amanzi_input)."));
  XMLString::transcode(element->getAttribute(xstr), str, 99);
  string type = str;
  Simulator* simulator = NULL;
  if (type == "structured")
  {
#ifdef ENABLE_Structured
    if (version[0] == '2')
      simulator = new AmanziStructuredGridSimulationDriver(doc);
    else // Assume version 1
    {
      assert(version[0] == '1');
      simulator = new AmanziStructuredGridSimulationDriver(input_filename);
    }
#else
    amanzi_throw(Errors::Message("Structured not supported in current build"));
#endif
  }
  else if (strcmp(str, "unstructured") == 0)
  {
#ifdef ENABLE_Unstructured
    if (version[0] == '2')
      simulator = new AmanziUnstructuredGridSimulationDriver(doc);
    else // Assume version 1
    {
      assert(version[0] == '1');
      simulator = new AmanziUnstructuredGridSimulationDriver(input_filename);
    }
#else
    amanzi_throw(Errors::Message("Unstructured not supported in current build"));
#endif
  }

  // Clean up.
  delete doc;
  XMLPlatformUtils::Terminate();

  return simulator;
}

} // end namespace SimulatorFactory
} // end namespace amanzi

