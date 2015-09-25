#include "errors.hh"
#include "exceptions.hh"

// TPLs
#include "xercesc/parsers/XercesDOMParser.hpp"

#ifdef ENABLE_Unstructured
#include "AmanziUnstructuredGridSimulationDriver.hh"
#endif
#ifdef ENABLE_Structured
#include "amanzi_structured_grid_simulation_driver.H"
#endif

#include "SimulatorFactory.hh"

using namespace xercesc;

namespace Amanzi
{

namespace SimulatorFactory
{

Simulator* Create(const std::string& input_filename)
{
  XMLPlatformUtils::Initialize();

  // Set up an XML DOM parser.
  XercesDOMParser* parser = new XercesDOMParser();
  parser->setExitOnFirstFatalError(true);
  parser->setValidationConstraintFatal(true);
  parser->setValidationScheme(XercesDOMParser::Val_Never);
  parser->setDoNamespaces(true);
  parser->setCreateCommentNodes(false);

  AmanziErrorHandler* errorHandler = new AmanziErrorHandler();
  parser->setErrorHandler(errorHandler);
  parser->useCachedGrammarInParse(true);
 
  try {
    parser->parse(input_filename.c_str());
  }
  catch (const OutOfMemoryException& e) {
    std::cerr << "OutOfMemoryException" << std::endl;
    Exceptions::amanzi_throw(Errors::Message("Ran out of memory while parsing the input file. Aborting."));
  }
  catch (...) {
    Exceptions::amanzi_throw(Errors::Message("Errors occured while parsing the input file. Aborting."));
  }

  // Open the document.
  DOMDocument* doc = parser->getDocument();

  // Determine whether this is a structured or unstructured input and dispatch 
  // accordingly.
  char str[100];
  XMLCh xstr[100];
  DOMElement* element = doc_->getDocumentElement();
  XMLString::transcode(element->getTagName(), str, 99);
  if (strcmp(str, "amanzi_input") != 0)
    Exceptions::amanzi_throw(Errors::Message("Invalid input file."));
  XMLString::transcode("type", xstr, 99);
  if (not element->hasAttribute(xstr)) 
    Exceptions::amanzi_throw(Errors::Message("Invalid input file (no 'type' attribute in amanzi_input)."));
  XMLString::transcode(element->getAttribute(xstr), str, 99);
  Simulator* simulator = NULL;
  if (strcmp(str, "structured") == 0)
#ifdef ENABLE_Structured
    simulator = new AmanziStructuredGridSimulationDriver(doc);
#else
    amanzi_throw(Errors::Message("Structured not supported in current build"));
#endif
  else if (strcmp(str, "unstructured") == 0)
    simulator = new AmanziUnstructuredGridSimulationDriver(doc);
#else
    amanzi_throw(Errors::Message("Unstructured not supported in current build"));
#endif

  // Clean up.
  delete doc;
  delete errorHandler;
  delete parser;
  XMLPlatformUtils::Terminate();

  return simulator;
}

} // end namespace SimulatorFactory
} // end namespace amanzi

