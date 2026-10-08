/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

A simple utility that reads Xml and writes Yaml using Teuchos/Trilinos
capabilities.

*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_YamlParser_decl.hpp"


int
main(int argc, char* argv[])
{
  // input and output
  std::string output_yaml, input_xml;
  if ((argc >= 3) && (argv[argc - 1][0] != '-')) {
    output_yaml = std::string(argv[argc - 1]);
    argc--;
  }
  if ((argc >= 2) && (argv[argc - 1][0] != '-')) {
    input_xml = std::string(argv[argc - 1]);
    argc--;
  }

  // commandline parser
  Teuchos::CommandLineProcessor clp;
  clp.setDocString(
    "Convert Xml file to Yaml\n\nStandard usage: xml_to_yaml input.xml output.yaml\n");

  clp.throwExceptions(false);
  clp.recogniseAllOptions(true);

  auto parseReturn = clp.parse(argc, argv);
  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    return 0;
  }
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return 1;
  }

  // make sure we got input -- output can be std::cout
  if (input_xml.empty()) {
    std::cerr << "ERROR: no input file provided" << std::endl;
    clp.printHelpMessage("xml_to_yaml", std::cerr);
  }

  // convert
  // -- read
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(input_xml);

  // -- write
  if (output_yaml.empty()) {
    Teuchos::YAMLParameterList::writeYamlStream(std::cout, *plist);
  } else {
    Teuchos::YAMLParameterList::writeYamlFile(output_yaml, *plist);
  }
  return 0;
}
