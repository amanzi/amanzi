#include "InputTranslator.hh"
#include "InputParserIS_Defs.hh"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#define  BOOST_FILESYTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

#include <xercesc/dom/DOM.hpp>
//#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
//#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/parsers/DOMLSParserImpl.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include "ErrorHandler.hpp"
//#include <xercesc/dom/DOMError.hpp>
//#include <xercesc/dom/DOMErrorHandler.hpp>
//#include <xercesc/sax/ErrorHandler.hpp>

XERCES_CPP_NAMESPACE_USE

namespace Amanzi {
namespace AmanziNewInput {

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
//Teuchos::ParameterList translate(const std::string& xmlfilename, const std::string& xmlSchemafile) {
Teuchos::ParameterList translate(const std::string& xmlfilename, std::string& spec) {

  Teuchos::ParameterList new_list;
  Teuchos::ParameterList def_list;
  
  
  // check if this is a new or old file
  // if old, create PL and send back
  // otherwise continue

  // do validation and parsing 
  // (by passing in xmlfile this moves the all of this here, 
  // we keep all the xerces and error handling out of main
  // which creates tons of clutter)
  // TODO: error handler
  XMLPlatformUtils::Initialize();
  //const char* schemaFile(xmlSchemafile.c_str());
  const char* xmlFile(xmlfilename.c_str());

  XercesDOMParser *parser = new XercesDOMParser();
  parser->setExitOnFirstFatalError(true);
  parser->setValidationConstraintFatal(true);
  parser->setValidationScheme(XercesDOMParser::Val_Never);
  parser->setDoNamespaces(true);
  //parser->setDoSchema(true);
  AmanziErrorHandler* errorHandler = new AmanziErrorHandler();
  parser->setErrorHandler(errorHandler); //EIB - commented out until Xerces update to handle XSD 1.1
  //parser->setExternalNoNamespaceSchemaLocation(schemaFile);
  //parser->loadGrammar(XMLString::transcode(schemaFile), Grammar::SchemaGrammarType, true);
  parser->useCachedGrammarInParse(true);
 
  bool errorsOccured = false;

  // EIB - this validation portion is basically useless until Xerces C++ implements XSD 1.1
  //       has only been implemented in Java version, no public plans to implement in C++
  //       adding in as much error checking in "get_..." sections as I can, when I can
  try {
      parser->parse(xmlfilename.c_str());
  }
  catch (const OutOfMemoryException& e) {
      std::cerr << "OutOfMemoryException" << std::endl;
      errorsOccured = true;
      Exceptions::amanzi_throw(Errors::Message("Ran out of memory while parsing the input file. Aborting."));
  }
  catch (...) {
      errorsOccured = true;
      Exceptions::amanzi_throw(Errors::Message("Errors occured while parsing the input file. Aborting."));
  }

  // go through each section, if it exist in the file, translate it 
  // to the old format
  DOMDocument *doc = parser->getDocument();

  // check that XML has version 2,x or version 1.x spec
  DOMElement* root = doc->getDocumentElement();
  char* tmp = XMLString::transcode(root->getTagName());
  spec = "";
  if (strcmp(tmp, "amanzi_input") == 0) {
    spec = "v2";
    // check that input is for Amanzi-U
    char* type = XMLString::transcode(root->getAttribute(XMLString::transcode("type")));
    if (strcmp(type, "unstructured") == 0) {
      Teuchos::ParameterList& echo_list = new_list.sublist("Echo Translated Input");
      echo_list.set<std::string>("Format", "unstructured_native");
      delete errorHandler;
      XMLPlatformUtils::Terminate();
      return new_list;
    }
  } else {
    if (strcmp(tmp, "ParameterList") == 0) spec = "v1";
    delete errorHandler;
    XMLPlatformUtils::Terminate();
    return new_list;
  }
  XMLString::release(&tmp);

  // grab the version number attribute
  new_list.set<std::string>("Amanzi Input Format Version", get_amanzi_version(doc,def_list));

  // store the filename in the def_list for later use
  def_list.set<std::string>("xmlfilename",xmlfilename);

  // grab the simulation type structured vs. unstructured
  get_sim_type(doc, &def_list);

  // grab the mesh type
  //new_list.sublist(framework) = ...;

  // grab verbosity early
  Teuchos::ParameterList verb;
  verb = get_verbosity(doc);
  def_list.sublist("simulation") = verb;
  voI_ = Teuchos::rcp(new VerboseObject("InputTranslator", verb));
  Teuchos::OSTab tab = voI_->getOSTab();
  
  def_list.sublist("constants") = get_constants(doc, def_list);

  new_list.sublist("General Description") = get_model_description(doc, def_list);
  new_list.sublist("Echo Translated Input") = get_echo_translated(doc, def_list);

  new_list.sublist("Mesh") = get_Mesh(doc, def_list);
  new_list.sublist("Domain").set<int>("Spatial Dimension",dimension_);
  new_list.sublist("Execution Control") = get_execution_controls(doc, &def_list);
  new_list.sublist("Phase Definitions") = get_phases(doc, def_list);
  new_list.sublist("Regions") = get_regions(doc, &def_list);
  new_list.sublist("Material Properties") = get_materials(doc, def_list);
  new_list.sublist("Initial Conditions") = get_initial_conditions(doc, def_list);
  new_list.sublist("Boundary Conditions") = get_boundary_conditions(doc, def_list);
  new_list.sublist("Sources") = get_sources(doc, def_list);
  new_list.sublist("Output") = get_output(doc, def_list);
  
  // hack to go back and add chemistry list (not geochemistry, for kd problem)
  if ( def_list.isSublist("Chemistry") ) {
    new_list.sublist("Chemistry") = def_list.sublist("Chemistry");
  } else { 
    new_list.sublist("Chemistry") = make_chemistry(def_list);
  }
  
  if (def_list.isParameter("petsc_options_file")) {
    new_list.set<std::string>("Petsc Options File",def_list.get<std::string>("petsc_options_file"));
  }

  delete errorHandler;
  XMLPlatformUtils::Terminate();
  //def_list.print(std::cout,true,false);

  // return the completely translated input file as a parameter list
  return new_list;
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_verbosity(DOMDocument* xmlDoc) {
    
    DOMNodeList* nodeList;
    DOMNode* nodeAttr;
    DOMNamedNodeMap* attrMap;
    char* textContent;
    
    // get execution contorls node
    nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("execution_controls"));
    Teuchos::ParameterList simPL;
  
    for (int i=0; i<nodeList->getLength(); i++) {
        DOMNode* ecNode = nodeList->item(i);
        if (DOMNode::ELEMENT_NODE == ecNode->getNodeType()) {
            //loop over children
            DOMNodeList* children = ecNode->getChildNodes();
            for (int j=0; j<children->getLength(); j++) {
                DOMNode* currentNode = children->item(j) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                    char* tagname = XMLString::transcode(currentNode->getNodeName());
                    if (std::string(tagname) == "verbosity") {
                        attrMap = currentNode->getAttributes();
                        nodeAttr = attrMap->getNamedItem(XMLString::transcode("level"));
			if (nodeAttr) {
                          textContent = XMLString::transcode(nodeAttr->getNodeValue());
                          simPL.sublist("VerboseObject").set<std::string>("Verbosity Level",trim_string(textContent));
                          simPL.set<std::string>("verbosity",trim_string(textContent));
			} else {
                          throw_error_illformed("verbosity", "value", "level");
			}
                        XMLString::release(&textContent);
                    }
                }
            }
        }
    }
    return simPL;
    

}
    
/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_constants(DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  Teuchos::ParameterList list;

  DOMNamedNodeMap *attrMap;
  DOMNode *namedNode;
  char* name;
  char* type;
  char* value;
  char* char_array;
  double time;
  Errors::Message msg;

  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *voI_->os() << "Getting Constants" << std::endl;
  }
  
  // read in new stuff
  DOMNodeList* nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("definitions"));

  if (nodeList->getLength() > 0) {
    DOMNode* nodeD = nodeList->item(0);
    DOMNodeList* childern = nodeD->getChildNodes();
    for (int i=0; i<childern->getLength(); i++) {
      DOMNode* currentNode = childern->item(i) ;
      if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
        char* tagname = XMLString::transcode(currentNode->getNodeName());
	// deal with: constants, named_times, macros
        if (std::string(tagname) == "constants") {
          DOMNodeList* kids = currentNode->getChildNodes();
          for (int j=0; j<kids->getLength(); j++) {
            DOMNode* currentKid = kids->item(j) ;
            if (DOMNode::ELEMENT_NODE == currentKid->getNodeType()) {
              char* kidname = XMLString::transcode(currentKid->getNodeName());
              // types: constant, time_constant, numerical_constant, area_mass_flux_constant
              if (std::string(kidname) == "constant") {
	        attrMap = currentKid->getAttributes();
	        namedNode = attrMap->getNamedItem(XMLString::transcode("name"));
		if (namedNode) {
	            name = XMLString::transcode(namedNode->getNodeValue());
		}
                else {
                  throw_error_illformed("definitions", "name", "constant");
		}
	        namedNode = attrMap->getNamedItem(XMLString::transcode("type"));
		if (namedNode) {
	          type = XMLString::transcode(namedNode->getNodeValue());
		}
                else {
                  throw_error_illformed("definitions", "type", name);
		}
	        namedNode = attrMap->getNamedItem(XMLString::transcode("value"));
		if (namedNode) {
	          value = XMLString::transcode(namedNode->getNodeValue());
		}
                else {
                  throw_error_illformed("definitions", "value", name);
		}
		if (std::string(type) == "time") {
		  // check if time and convert to seconds - year = 365.25
		  // TODO: EIB - verify this works with spaces
		  // TODO: EIB - expect Akuna to move to no deliminator, need to test for this
                  time = convert_time_value(value);
                }
                else {
		  time = atof(value);
		}
		// add to list
		Teuchos::ParameterList tmp;
		tmp.set<std::string>("type",trim_string(type));
		tmp.set<double>("value",time);
		list.sublist("constants").sublist(name) = tmp;
                XMLString::release(&name);
                XMLString::release(&type);
                XMLString::release(&value);
	      }
              else if (std::string(kidname) == "time_constant") {
	        attrMap = currentKid->getAttributes();
	        namedNode = attrMap->getNamedItem(XMLString::transcode("name"));
		if (namedNode) {
	          name = XMLString::transcode(namedNode->getNodeValue());
		}
                else {
                  throw_error_illformed("definitions", "name", "time_constant");
		}
	        namedNode = attrMap->getNamedItem(XMLString::transcode("value"));
		if (namedNode) {
	          value = XMLString::transcode(namedNode->getNodeValue());
		}
                else {
                  throw_error_illformed("definitions", "value", name);
		}
		// check if time and convert to seconds - year = 365.25
		// TODO: EIB - verify this works with spaces
		// TODO: EIB - expect Akuna to move to no deliminator, need to test for this
		char_array = strtok(value,";, ");
		time = atof(char_array);
		char_array = strtok(NULL,";,");
		if (strcmp(char_array,"y")==0) { time = time*365.25*24.0*60.0*60.0; }
		else if (strcmp(char_array,"d")==0) { time = time*24.0*60.0*60.0; }
		else if (strcmp(char_array,"h")==0) { time = time*60.0*60.0; }
		// add to list
		Teuchos::ParameterList tmp;
		tmp.set<double>("value",time);
		list.sublist("time_constants").sublist(name) = tmp;
                XMLString::release(&name);
                XMLString::release(&value);
	      }
              else if (std::string(kidname) == "numerical_constant") {
	        attrMap = currentKid->getAttributes();
	        namedNode = attrMap->getNamedItem(XMLString::transcode("name"));
		if (namedNode) {
	          name = XMLString::transcode(namedNode->getNodeValue());
		}
                else {
                  throw_error_illformed("definitions", "name", "numerical_constant");
		}
	        namedNode = attrMap->getNamedItem(XMLString::transcode("value"));
		if (namedNode) {
	          value = XMLString::transcode(namedNode->getNodeValue());
		}
                else {
                  throw_error_illformed("definitions", "value", name);
		}
		// add to list
		Teuchos::ParameterList tmp;
		tmp.set<double>("value",atof(value));
		list.sublist("numerical_constant").sublist(name) = tmp;
                XMLString::release(&name);
                XMLString::release(&value);
	      }
              else if (std::string(kidname) == "area_mass_flux_constant") {
	        attrMap = currentKid->getAttributes();
	        namedNode = attrMap->getNamedItem(XMLString::transcode("name"));
		if (namedNode) {
	          name = XMLString::transcode(namedNode->getNodeValue());
		}
                else {
                  throw_error_illformed("definitions", "name", "area_mass_flux_constant");
		}
	        namedNode = attrMap->getNamedItem(XMLString::transcode("value"));
		if (namedNode) {
	          value = XMLString::transcode(namedNode->getNodeValue());
		}
                else {
                  throw_error_illformed("definitions", "value", name);
		}
		// add to list
		Teuchos::ParameterList tmp;
		tmp.set<double>("value",atof(value));
		list.sublist("area_mass_flux_constant").sublist(name) = tmp;
                XMLString::release(&name);
                XMLString::release(&value);
	      }
              XMLString::release(&kidname);
	    }
	  }
	}
        else if (std::string(tagname) == "named_times") {
	  //TODO: EIB - deal with named times
          DOMNodeList* kids = currentNode->getChildNodes();
          for (int j=0; j<kids->getLength(); j++) {
            DOMNode* currentKid = kids->item(j) ;
            if (DOMNode::ELEMENT_NODE == currentKid->getNodeType()) {
              char* kidname = XMLString::transcode(currentKid->getNodeName());
              // types: time
              if (std::string(kidname) == "time") {
	      }
              XMLString::release(&kidname);
	    }
	  }
	//} else if (strcmp(tagname,"macros")==0) {
	  //TODO: EIB - move macros here from outputs
        }
        XMLString::release(&tagname);
      }
    }
  }

  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
std::string get_amanzi_version(DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {
  std::stringstream old_version;
  
  //XMLCh* tag = XMLString::transcode("amanzi_input");

  DOMNodeList* nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("amanzi_input"));
  //XMLString::release(&tag);

  const XMLSize_t nodeCount = nodeList->getLength();  
  if (nodeList->getLength() > 0) {
    DOMNode* nodeGD = nodeList->item(0);
    DOMElement* elementGD = static_cast<DOMElement*>(nodeGD);
    std::string version(XMLString::transcode(elementGD->getAttribute(XMLString::transcode("version"))));
    
    int major, minor, micro;
    
    std::stringstream ss;
    ss << version;
    std::string ver;
    
    try {
      getline(ss,ver,'.');
      major = boost::lexical_cast<int>(ver);
      
      getline(ss,ver,'.');
      minor = boost::lexical_cast<int>(ver);
      
      getline(ss,ver);
      micro = boost::lexical_cast<int>(ver);
    }
    catch (...) {
      Exceptions::amanzi_throw(Errors::Message("The version string in the input file '"+version+"' has the wrong format, please use X.Y.Z, where X, Y, and Z are integers."));
    }

    if ( (major == AMANZI_INPUT_VERSION_MAJOR) && (minor == AMANZI_INPUT_VERSION_MINOR) && (micro == AMANZI_INPUT_VERSION_MICRO) ) {
      // now we can proceed, we translate to a v1.2.2 parameterlist
      old_version << AMANZI_OLD_INPUT_VERSION_MAJOR <<"."<< AMANZI_OLD_INPUT_VERSION_MINOR <<"."<< AMANZI_OLD_INPUT_VERSION_MICRO; 
    }
    else {
      std::stringstream ver;
      ver << AMANZI_INPUT_VERSION_MAJOR << "." << AMANZI_INPUT_VERSION_MINOR << "." << AMANZI_INPUT_VERSION_MICRO;
      Errors::Message msg;
      msg << "The input version " << version << " specified in the input file is not supported. This version of amanzi supports version "<< ver.str() << ".\n";
      if ((major == 2) && (minor == 1) && (micro == 0)) {
        msg << "  The python script UpdateSpec_210-211.py will update a 2.1.0 input file to 2.1.1.";
	msg <<"   The script is located in the source tree in tools/input and is installed in $INSTALL/bin";
      }
      Exceptions::amanzi_throw(msg);
    }
  }
  else {
    // amanzi inpurt description did not exist, error
    Exceptions::amanzi_throw(Errors::Message("Amanzi input description does not exist <amanzi_input version=...>"));
  }

  return old_version.str();
}

/* ****************************************************************
 * Empty
 ******************************************************************
 */
void get_sim_type(DOMDocument* xmlDoc, Teuchos::ParameterList* def_list) {
  std::stringstream old_version;
  
  Errors::Message msg;
  DOMNodeList* nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("amanzi_input"));

  const XMLSize_t nodeCount = nodeList->getLength();  
  if (nodeList->getLength() > 0) {
    DOMNode* nodeGD = nodeList->item(0);
    DOMElement* elementGD = static_cast<DOMElement*>(nodeGD);
    std::string sim_type = (XMLString::transcode(elementGD->getAttribute(XMLString::transcode("type"))));
    if (sim_type.length() > 0) {
      def_list->set<std::string>("sim_type",sim_type);
      if (strcmp(sim_type.c_str(),"structured")==0) {
	  isUnstr_ = false;
      }
      else {
	  isUnstr_ = true;
      }
    }
    else {
      throw_error_illformed("amanzi_input", "attribute", "type", "structured or unstructured");
    }
    
  }
  else {
    // amanzi input description did not exist, error
    throw_error_illformed("amanzi_input", "attribute", "type", "structured or unstructured");
  }

}
    
/* 
 ******************************************************************
 * Empty
 ****************************************************************** 
 */
Teuchos::ParameterList get_model_description(DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  Teuchos::ParameterList list;

  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *voI_->os() << "Getting Model Description" << std::endl;
  }
  
  // read in new stuff
  XMLCh* tag = XMLString::transcode("model_description");
  DOMNodeList* nodeList = xmlDoc->getElementsByTagName(tag);
  XMLString::release(&tag);

  // write to old format, mostly won't be read so just write out as strings
  const XMLSize_t nodeCount = nodeList->getLength() ;
  if (nodeList->getLength() > 0) {
    DOMNode* nodeGD = nodeList->item(0);
    DOMElement* elementGD = static_cast<DOMElement*>(nodeGD);
    char* model_id = XMLString::transcode(elementGD->getAttribute(XMLString::transcode("name")));
    list.set<std::string>("model_id",trim_string(model_id));

    DOMNodeList* childern = nodeGD->getChildNodes();
    for (int i=0; i<childern->getLength(); i++) {
      DOMNode* currentNode = childern->item(i) ;
      if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
        char* tagname = XMLString::transcode(currentNode->getNodeName());
        DOMNode::NodeType type = currentNode->getNodeType();
        if (strcmp(tagname,"units")!=0) {
          char* textContent = XMLString::transcode(currentNode->getTextContent());
          list.set<std::string>(tagname,trim_string(textContent));
          XMLString::release(&textContent);
        }
        XMLString::release(&tagname);
        }
    }
  }
  else {
	  // model_description didn't exist, report an error
  }

  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ****************************************************************** 
 */
Teuchos::ParameterList get_echo_translated(DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *voI_->os() << "Getting Echo" << std::endl;
  }
  
  XMLCh* TAG_echo = XMLString::transcode("echo_translated_input");
  XMLCh* ATTR_filename = XMLString::transcode("file_name");
  XMLCh* ATTR_format = XMLString::transcode("format");

  // Defaults
  std::string translated_filename = "translated_input.xml"; // NOTE: default changes for v1.  See below
  std::string translated_format = "v1";
  bool echo_translated = false;

  DOMNodeList* nodeList = xmlDoc->getElementsByTagName(TAG_echo);

  if (nodeList->getLength() > 0 ) {

    echo_translated = true;

    DOMElement* tmpElement = static_cast<DOMElement*>(nodeList->item(0));
    if (tmpElement->hasAttribute(ATTR_format)) {
      char* textContent = XMLString::transcode(tmpElement->getAttribute(ATTR_format));
      translated_format = trim_string(textContent);
      XMLString::release(&textContent);
    }

    // For v1, set default filename, which can be overridden
    if (translated_format == "v1") {
      translated_filename = def_list.get<std::string>("xmlfilename");
      std::string new_extension("_oldspec.xml");
      size_t pos = translated_filename.find(".xml");
      translated_filename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)12);
    }

    if (tmpElement->hasAttribute(ATTR_filename)) {
      char* textContent = XMLString::transcode(tmpElement->getAttribute(ATTR_filename));
      translated_filename = trim_string(textContent);
      XMLString::release(&textContent);
    }
  }

  XMLString::release( &TAG_echo );
  XMLString::release( &ATTR_filename );
  XMLString::release( &ATTR_format );

  Teuchos::ParameterList list;

  if (echo_translated) {
    list.set<std::string>("File Name",translated_filename);
    list.set<std::string>("Format",translated_format);
  }

  return list;
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_Mesh(DOMDocument* xmlDoc, Teuchos::ParameterList def_list ) {

  //TODO: EIB - see if mesh list ending up in right order/structure

  Teuchos::ParameterList list;

  bool generate = true;
  bool read = false;
  char *framework;
  Teuchos::ParameterList mesh_list;
  bool all_good = false;
  Errors::Message msg;
  std::stringstream helper;
    
  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Mesh" << std::endl;
  }

  Teuchos::RCP<Teuchos::StringValidator> meshfileValidator = rcp (new Teuchos::StringValidator(meshfileStrings));
  
  // read in new stuff
  XMLCh* tag = XMLString::transcode("mesh");
  DOMNodeList* nodeList = xmlDoc->getElementsByTagName(tag);
  XMLString::release(&tag);

  // read the attribute to set the framework sublist
  const XMLSize_t nodeCount = nodeList->getLength() ;
  if (nodeList->getLength() > 0) {
    DOMNode* nodeMesh = nodeList->item(0);
    DOMElement* elementMesh = static_cast<DOMElement*>(nodeMesh);
    // if unstructure, look for framework attribute
    if (isUnstr_) { 
      if(elementMesh->hasAttribute(XMLString::transcode("framework"))) {
        framework = XMLString::transcode(elementMesh->getAttribute(XMLString::transcode("framework")));
      } else { 
        msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing mesh - " ;
        msg << "framework was missing or ill-formed. \n  Use default framework='mstk' if unsure.  Please correct and try again \n" ;
        Exceptions::amanzi_throw(msg);
      }
    }

    // loop over child nodes
    DOMNodeList* children = nodeMesh->getChildNodes();
    // first figure out what the dimension is
    all_good = false;
    for (int i=0; i<children->getLength(); i++) {
      DOMNode* currentNode = children->item(i) ;
      if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
	char* tagname = XMLString::transcode(currentNode->getNodeName());
	if (strcmp(tagname,"dimension")==0) {
	  char* temp = XMLString::transcode(currentNode->getTextContent());
	  if (strlen(temp) > 0) {
	      dimension_ = get_int_constant(temp,def_list);
	      all_good = true;
	  }
          XMLString::release(&temp);
	}
        XMLString::release(&tagname);
      }
    }
    if (!all_good) {
      throw_error_illformed("mesh", "element", "dimension");
    }

    // now we can properly parse the generate/read list
    all_good = false;
    for (int i=0; i<children->getLength(); i++) {
      DOMNode* currentNode = children->item(i) ;
      if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
	char* tagname = XMLString::transcode(currentNode->getNodeName());   

	if (strcmp(tagname,"generate")==0) {
	  all_good = true;
	  generate = true;
	  read = false;
	  DOMElement* elementGen = static_cast<DOMElement*>(currentNode);

	  // get Number of Cells
	  Teuchos::Array<int> ncells; 
	  DOMNodeList* nodeList = elementGen->getElementsByTagName( XMLString::transcode("number_of_cells"));
	  DOMNode* node = nodeList->item(0);
	  DOMElement* elementNode = static_cast<DOMElement*>(node);
	  DOMNamedNodeMap *attrMap = node->getAttributes();
          DOMNode* nodeAttr;
          char* attrName;
	  char* temp;
	  // make sure number of attributes equals dimension
	  if ( attrMap->getLength() == dimension_) {
	    // loop over attributes to get nx, ny, nz as needed
	    for (int j=0; j<attrMap->getLength(); j++) {
	      nodeAttr = attrMap->item(j);
	      attrName =XMLString::transcode(nodeAttr->getNodeName());
	      if (attrName) {
	        temp = XMLString::transcode(nodeAttr->getNodeValue());
	        if (strlen(temp) > 0) {
	        ncells.append(get_int_constant(temp,def_list));
		} else {
		  all_good = false;
	          helper << "number_of_cells "<<attrName;
	        }
                XMLString::release(&temp);
	      } else {
		all_good = false;
	        helper << "number_of_cells "<<attrName;
	      }
	      XMLString::release(&attrName);
	    }
            mesh_list.set<Teuchos::Array<int> >("Number of Cells",ncells);
	  } else {
	    helper << "number_of_cells";
	    all_good = false;
	  }

	  // get Box - generalize
	  //char* char_array;
	  nodeList = elementGen->getElementsByTagName( XMLString::transcode("box"));
	  node = nodeList->item(0);
	  elementNode = static_cast<DOMElement*>(node);
	  temp = XMLString::transcode(elementNode->getAttribute( XMLString::transcode("low_coordinates")));
	  if (strlen(temp) > 0) {
	    // translate to array
	    Teuchos::Array<double> low = make_coordinates(temp, def_list);
            mesh_list.set<Teuchos::Array<double> >("Domain Low Coordinate",low);
	    if (low.length() != dimension_) {
	      helper << "low_coordinates";
	      all_good = false;
	    }
	  }
          else {
	    helper << "low_coordinates";
	    all_good = false;
	  }
	  XMLString::release(&temp);
	  temp = XMLString::transcode(elementNode->getAttribute( XMLString::transcode("high_coordinates")));
	  if (strlen(temp) > 0) {
	    // translate to array
	    Teuchos::Array<double> high = make_coordinates(temp, def_list);
            mesh_list.set<Teuchos::Array<double> >("Domain High Coordinate",high);
	    if (high.length() != dimension_) {
	      helper << "high_coordinates";
	      all_good = false;
	    }
	  }
          else {
	    helper << "high_coordinates";
	    all_good = false;
	  }
	  XMLString::release(&temp);
	}

	else if (strcmp(tagname,"read")==0) {
	  read = true;
	  generate = false;
	  bool goodtype = false;
	  bool goodname = false;
	  DOMElement* elementRead = static_cast<DOMElement*>(currentNode);

	  char* value = XMLString::transcode(elementRead->getElementsByTagName(
				  XMLString::transcode("format"))->item(0)->getTextContent());
          std::string format(trim_string(value));
          if (boost::iequals(format, "exodus ii")) {
            mesh_list.set<std::string>("Format","Exodus II");
            goodtype = true;
	  }
          else if (boost::iequals(format, "h5m")) {
            mesh_list.set<std::string>("Format","H5M");
            goodtype = true;
          }
          else {
            mesh_list.set<std::string>("Format",format,"Format of meshfile",meshfileValidator);
          }
	  char* filename = XMLString::transcode(elementRead->getElementsByTagName(
				  XMLString::transcode("file"))->item(0)->getTextContent());
	  if (strlen(filename) > 0) {
              mesh_list.set<std::string>("File",trim_string(filename));
	      goodname = true;
	  }
	  XMLString::release(&value);
	  XMLString::release(&filename);
	  if (goodtype && goodname) all_good = true;

	}
	//EIB - handles legacy files
	//TODO:: EIB - remove to conform to updated schema
	else if (strcmp(tagname,"file")==0) {
	  read = true;
	  generate = false;
	  char* filename = XMLString::transcode(currentNode->getTextContent());
	  if (strlen(filename) > 0) {
            mesh_list.set<std::string>("File",trim_string(filename));
	    mesh_list.set<std::string>("Format","Exodus II");
	    all_good = true;
	    std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
	    std::cout << "    Please note - the XML Schema for specifing a mesh file to read has been updated." << std::endl;
	    std::cout << "    See the Amanzi User Guide for the latest information." << std::endl;
	    std::cout <<"     The legacy format is being handled for now.  Please update input files for future versions." << std::endl;
	  }
	  XMLString::release(&filename);
	}

      }
    }
    if (!all_good) {
      throw_error_illformed("mesh", helper.str(), "generate/read");
    }
    
    if (generate || read) {
      if (isUnstr_) {
        if (strcmp(framework,"mstk")==0) {
          list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MSTK");
        } else if (strcmp(framework,"moab")==0) {
          list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MOAB");
        } else if (strcmp(framework,"simple")==0) {
          list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","Simple");
        } else if (strcmp(framework,"stk::mesh")==0) {
          list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","stk::mesh");
        } else {
          //list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MSTK");
          msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing mesh - " ;
          msg << "unknown framework=" << framework << ". See the schema for acceptable types. \n  Please correct and try again \n" ;
          Exceptions::amanzi_throw(msg);
	}
      }
    }

    if (generate) {
      if (isUnstr_) {
        list.sublist("Unstructured").sublist("Generate Mesh").sublist("Uniform Structured") = mesh_list;
      } else {
        list.sublist("Structured") = mesh_list;
      }
    }
    else if (read) {
      if (isUnstr_) {
	list.sublist("Unstructured").sublist("Read Mesh File") = mesh_list;
      }
    }
    else {
      // bad mesh, again if validated shouldn't need this
      msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing mesh - " ;
      msg << "\n  Please correct and try again \n" ;
      Exceptions::amanzi_throw(msg);
    }
  }
  else {
    // no mesh sub-elements
    throw_error_illformed("mesh", "element", "framework");
  }

  //XMLString::release(&framework);

  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_execution_controls(DOMDocument* xmlDoc, Teuchos::ParameterList* def_list ) {

  Teuchos::ParameterList list;

  DOMNodeList* nodeList;
  DOMNode* nodeTmp;
  DOMNode* nodeAttr;
  DOMNamedNodeMap* attrMap;
  char* tagName;
  char* attrName;
  char* textContent;
  char* elemContent;
  std::string meshbase;
  std::string value;
  Errors::Message msg;

  // This actually includes: process kernels, execution controls, and numerical controls
  // all three map back to the old exection controls
  
  if (isUnstr_) {
    meshbase = std::string("Unstructured Algorithm");
  } else {
    meshbase = std::string("Structured Algorithm");
  }

    // set so we don't have to reread transport and chemisty list later
  Teuchos::ParameterList fpkPL, tpkPL, cpkPL;
  bool flowON=false;
  bool staticflowON=false;
  bool transportON=false;
  bool chemistryON=false;

  // get process kernels node
  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Process Kernels" << std::endl;
  }
  
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("process_kernels"));
  for (int i=0; i<nodeList->getLength(); i++) {
    DOMNode* pkNode = nodeList->item(i);
    if (DOMNode::ELEMENT_NODE == pkNode->getNodeType()) {
      DOMElement* pkElement = static_cast<DOMElement*>(pkNode);
      // get flow
      DOMNodeList* tmpList = pkElement->getElementsByTagName(XMLString::transcode("flow"));
      if (tmpList->getLength() > 0 ) {
        DOMElement* flowElement = static_cast<DOMElement*>(tmpList->item(0));
        if (flowElement->hasAttribute((XMLString::transcode("state")))) {
	  textContent = XMLString::transcode(flowElement->getAttribute(XMLString::transcode("state")));
	  value = trim_string(textContent);
	  value[0] = std::tolower(value[0]);
          if (value == "off"){
            list.set<std::string>("Flow Model","Off");
          }
	  else if (value == "on"){
            flowON = true;
            char* textContent2 = XMLString::transcode(flowElement->getAttribute(XMLString::transcode("model")));
	    std::string value2 = trim_string(textContent2);
	    value2[0] = std::tolower(value2[0]);
            def_list->set<std::string>("flow",value2);
            if (value2 == "saturated") {
              list.set<std::string>("Flow Model","Single Phase");
            }
            else if (value2 =="richards") {
              list.set<std::string>("Flow Model","Richards");
            }
            else if (value2 == "constant") {
              // EIB - also need to set integration mode = transient with static flow
              list.set<std::string>("Flow Model","Single Phase");
              staticflowON = true;
            }
            else {
              throw_error_illformed("process_kernels", "model", "flow");
            }
            XMLString::release(&textContent2);
          }
          else {
            throw_error_illformed("process_kernels", "state", "flow");
          }
          XMLString::release(&textContent);
        }
        else {
          throw_error_illformed("process_kernels", "state", "flow");
        }
      }
      else {
        throw_error_illformed("process_kernels", "element", "flow");
      }

      // get transport
      tmpList = pkElement->getElementsByTagName(XMLString::transcode("transport"));
      DOMElement* transElement = static_cast<DOMElement*>(tmpList->item(0));
      if (transElement->hasAttribute((XMLString::transcode("state")))) {
        textContent = XMLString::transcode(transElement->getAttribute(XMLString::transcode("state")));
	value = trim_string(textContent);
	value[0] = std::tolower(value[0]);
        def_list->set<bool>("transport",false);
        if (value == "off"){
          list.set<std::string>("Transport Model","Off");
        }
        else {
          list.set<std::string>("Transport Model","On");
          transportON=true;
          def_list->set<bool>("transport",true);
        }
        XMLString::release(&textContent);
      }

      // get chemisty - TODO: EIB - assuming this will be set to OFF!!!!!
      // NOTE: EIB - old spec options seem to be ON/OFF, algorithm option goes under Numerical Control Parameters
      tmpList = pkElement->getElementsByTagName(XMLString::transcode("chemistry"));
      attrMap = tmpList->item(0)->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("state"));
      if (nodeAttr) {
        attrName = XMLString::transcode(nodeAttr->getNodeValue());
      }
      else {
        throw_error_missattr("process_kernels", "attribute", "state", "chemistry");
      }
      value = trim_string(attrName);
      value[0] = std::tolower(value[0]);
      if (value == "off"){
        list.set<std::string>("Chemistry Model","Off");
      }
      else {
    	chemistryON=true;
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("engine"));
	if (nodeAttr) {
          textContent = XMLString::transcode(nodeAttr->getNodeValue());
        }
        else {
          throw_error_missattr("process_kernels", "attribute", "engine", "chemistry");
        }

	if (strcmp(textContent,"amanzi")==0) {
          list.set<std::string>("Chemistry Model","Amanzi");
          def_list->set<std::string>("chemistry_engine","amanzi");
	}
        else if (strcmp(textContent,"pflotran")==0) {
          list.set<std::string>("Chemistry Model","Alquimia");
          def_list->set<std::string>("chemistry_engine","pflotran");
          
          // start chemistry list in def_list for later use - set engine
          Teuchos::ParameterList chemPL;
          if (def_list->isSublist("chemistry_PL")) {
            chemPL = def_list->sublist("chemistry_PL");
          }
          chemPL.set<std::string>("Engine","PFloTran");
          // get input file name
          DOMNode* nodeAttr1 = attrMap->getNamedItem(XMLString::transcode("input_filename"));
          if (nodeAttr1) {
            char* attrName1 = XMLString::transcode(nodeAttr1->getNodeValue());
            chemPL.set<std::string>("Engine Input File",attrName1);
          }
          // attach list if it didn't exist
          def_list->sublist("chemistry_PL") = chemPL;
	}
        else {
          //TODO: EIB - error handle here!!!
        }
        def_list->set<std::string>("chemistry_engine",textContent);
        XMLString::release(&textContent);
      }
      XMLString::release(&attrName);
    }
  }
 
  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Execution Controls" << std::endl;
  }

  // get execution contorls node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("execution_controls"));
  Teuchos::ParameterList ecsPL;
  Teuchos::ParameterList defPL;
  bool hasSteady = false;
  bool hasTrans = false;
  bool hasRestart = false;
  int numControlPeriods = 0;
  Teuchos::Array<double> start_times;
  double sim_start=-1.;
  double sim_end=-1.;
  Teuchos::ParameterList simPL;

  Teuchos::RCP<Teuchos::StringValidator> verbosityValidator = rcp (new Teuchos::StringValidator(verbosityStrings));
  
  for (int i=0; i<nodeList->getLength(); i++) {
    DOMNode* ecNode = nodeList->item(i);
    if (DOMNode::ELEMENT_NODE == ecNode->getNodeType()) {
      //loop over children
      DOMNodeList* children = ecNode->getChildNodes();
      for (int j=0; j<children->getLength(); j++) {
        DOMNode* currentNode = children->item(j) ;
        if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
    	  char* tagname = XMLString::transcode(currentNode->getNodeName());
          if (strcmp(tagname,"verbosity")==0) {
            attrMap = currentNode->getAttributes();
            nodeAttr = attrMap->getNamedItem(XMLString::transcode("level"));
            if (nodeAttr) {
              textContent = XMLString::transcode(nodeAttr->getNodeValue());
            }
            else {
              throw_error_missattr("execution_controls", "attribute", "level", "verbosity");
            }
            // This is a hack since structured requires uppercase first character
            std::string value(trim_string(textContent));
            value[0] = std::toupper(value[0]);
            if (value == "Default") {
              value = VERBOSITY_DEFAULT;
            }
            list.set<std::string>("Verbosity",value,"Verbosity level",verbosityValidator);
            simPL.set<std::string>("verbosity",value);
            XMLString::release(&textContent);

	  }
          else if (strcmp(tagname,"execution_control_defaults")==0) {
            attrMap = currentNode->getAttributes();
            for (int k=0; k<attrMap->getLength(); k++) {
              nodeAttr = attrMap->item(k);
              attrName =XMLString::transcode(nodeAttr->getNodeName());
              textContent = XMLString::transcode(nodeAttr->getNodeValue());
              defPL.set<std::string>(attrName,trim_string(textContent));
              // EIB:: only include default mode if > 1 EC
              //if (strcmp(attrName,"mode")==0) {
              //  if (strcmp(textContent,"steady")==0) {
              //    hasSteady = true;
              //  } else {
              //    hasTrans = true;
              //  }
              //}
            }
	  }
          else if (strcmp(tagname,"execution_control")==0) {
            Teuchos::ParameterList ecPL;
            numControlPeriods++;
            attrMap = currentNode->getAttributes();
            char* name;
            bool saveName=true;
            for (int k=0; k<attrMap->getLength(); k++) {
              nodeAttr = attrMap->item(k);
              attrName =XMLString::transcode(nodeAttr->getNodeName());
              textContent = XMLString::transcode(nodeAttr->getNodeValue());
              ecPL.set<std::string>(attrName,trim_string(textContent));
              if (strcmp(attrName,"start")==0) {
                if (saveName) name=textContent;
                double time = get_time_value(textContent, *def_list);
                if (start_times.length() == 0) {  // if first time through
                  start_times.append(time);
                  sim_start = time;
                }
                else {
                  if (time < sim_start) sim_start = time;           // check for simulation start time
                  if (time >= start_times[start_times.length()-1]) { // if already sorted
                    start_times.append(time);
                  }
                  else {                                          // else, sort
                    int idx = start_times.length()-1;
                    Teuchos::Array<double> hold;
                    hold.append(start_times[idx]);
                    start_times.remove(idx);
                    idx--;
                    while (time < start_times[idx]) {
                      hold.append(start_times[idx]);
                      start_times.remove(idx);
                      idx--;
                    }
                    start_times.append(time);
                    for (int i=0; i<hold.length(); i++) {
                      idx = hold.length()-1-i;
                      start_times.append(hold[hold.length()-idx]);
                    }
                  }
                }
              }
              if (strcmp(attrName,"end")==0) {
                double time = get_time_value(textContent, *def_list);
                if (time > sim_end) sim_end = time;                   // check for simulation end time
              }
              if (strcmp(attrName,"mode")==0) {
                if (strcmp(textContent,"steady")==0) {
                  hasSteady = true;
                  saveName = false;
                  name = textContent;
                }
                else {
                  hasTrans = true;
                }
              }
              if (strcmp(attrName,"restart")==0) {
                hasRestart = true;
                //name = attrName;
                //ecsPL.sublist("restart") = ecPL;
              }
            }
            if (hasRestart && ecPL.isParameter("restart")) ecsPL.sublist("restart") = ecPL;
            ecsPL.sublist(name) = ecPL;
	  }
	}
      }
    }
  }
  // do an end time check 
  if (sim_end == -1.) {
    sim_end = start_times[start_times.length()-1];
  }
  simPL.set<double>("simulation_start",sim_start);
  simPL.set<double>("simulation_end",sim_end);
  simPL.set<Teuchos::Array<double> >("simulation_start_times",start_times);
  //ecsPL.sublist("simulation") = simPL;
  def_list->sublist("simulation") = simPL;

  // If > 1 EC, then include default mode || mode wasn't set in any execution control statement
  if (numControlPeriods > 1 || (!hasTrans && !hasSteady)) {
    if (defPL.isParameter("mode")) {
      if (defPL.get<std::string>("mode") == "steady") hasSteady = true;
      if (defPL.get<std::string>("mode") == "transient") hasTrans = true;
    }
  }

  // Now, go back and sort things out
  bool haveSSF = false; // have steady-steady incr/red factors for later
  bool haveTF = false;  // have transient incr/red factors for later
  
  // Restart
  if (hasRestart) {
    std::string value = ecsPL.sublist("restart").get<std::string>("restart");
    Teuchos::ParameterList restartPL;
    /* EIB: proposed v1.2.2 update - Change Restart name */
    //restartPL.set<std::string>("Checkpoint Data File Name",value);
    //list.sublist("Restart from Checkpoint Data File") = restartPL;
    restartPL.set<std::string>("File Name",value);
    list.sublist("Restart") = restartPL;
  }

  // Translate to "Time Period Control"
  
  Teuchos::Array<double> start_array;
  Teuchos::Array<double> init_array;
  Teuchos::Array<double> max_array;
  
  // loop over non-restart entries to get list of time step mins/maxs
  for (int idx = 0; idx < start_times.length(); idx++) {
    for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
      if (it->first != "restart") {
        double time = get_time_value(it->first, *def_list);
        if (time == start_times[idx]) {
          // look for init_dt in current execution control
          if ( ecsPL.sublist(it->first).isParameter("init_dt") ) {
            init_array.append(get_time_value(ecsPL.sublist(it->first).get<std::string>("init_dt"), *def_list));
          }
          // if not there, look for default exectution control
          else if ( defPL.isParameter("init_dt") ){
            init_array.append(get_time_value(defPL.get<std::string>("init_dt"), *def_list));
          }
          // deremine mode and get defaults from InputParserIS_Def.hh
          else {
            // look for mode in current execution control
            if ( ecsPL.sublist(it->first).isParameter("mode") ) {
              if (ecsPL.sublist(it->first).get<std::string>("mode") == "steady") {
                init_array.append(ST_MIN_TS);
              }
              else {
                init_array.append(TR_MIN_TS);
              }
            }
            // else look in default execution contorl
            else if (defPL.isParameter("mode") ){
              if (defPL.get<std::string>("mode") == "steady") {
                init_array.append(ST_MIN_TS);
              }
              else {
                init_array.append(TR_MIN_TS);
              }
            }
            // finally assume min of steady and transient (shouldn't get to here)
            else {
              init_array.append(ST_MIN_TS);
            }
          }
          // repeat for max_dt
          if ( ecsPL.sublist(it->first).isParameter("max_dt") ) {
            max_array.append(get_time_value(ecsPL.sublist(it->first).get<std::string>("max_dt"), *def_list));
          }
          // if not there, look for default exectution control
          else if ( defPL.isParameter("max_dt") ){
            max_array.append(get_time_value(defPL.get<std::string>("max_dt"), *def_list));
          }
          // deremine mode and get defaults from InputParserIS_Def.hh
          else {
            // look for mode in current execution control
            if ( ecsPL.sublist(it->first).isParameter("mode") ) {
              if (ecsPL.sublist(it->first).get<std::string>("mode") == "steady") {
                max_array.append(ST_MAX_TS);
              }
              else {
                max_array.append(TR_MAX_TS);
              }
            }
            // else look in default execution contorl
            else if (defPL.isParameter("mode") ){
              if (defPL.get<std::string>("mode") == "steady") {
                max_array.append(ST_MAX_TS);
              }
              else {
                max_array.append(TR_MAX_TS);
              }
            }
            // finally assume min of steady and transient (shouldn't get to here)
            else {
              max_array.append(ST_MAX_TS);
            }
          }
        }
      }
    }
  }
  // add to "Time Period Control" in order
  list.sublist("Time Period Control").set<Teuchos::Array<double> >("Start Times",start_times);
  list.sublist("Time Period Control").set<Teuchos::Array<double> >("Initial Time Step",init_array);
  list.sublist("Time Period Control").set<Teuchos::Array<double> >("Maximum Time Step",max_array);
  
  // add default entry, if it exists
  if (defPL.isParameter("init_dt")) {
    list.sublist("Time Period Control").set<double>("Default Initial Time Step",get_time_value(defPL.get<std::string>("init_dt"), *def_list));
  }
  
  
  // Steady case
  if (hasSteady && !hasTrans) {
    if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
        *voI_->os() << "Creating Steady State Execution Control" << std::endl;
    }
    
    // NOTE:: if you edit the Steady case, it is repeated under the Initialize to Steady case, so edit there too!!!!
    Teuchos::ParameterList steadyPL;
    // look for values from default list
    // if not there, grab from ec list
    std::string value;
    std::string method;
    bool gotValue;
    if (defPL.isParameter("start")) {
      value = defPL.get<std::string>("start");
      gotValue = true;
    }
    else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	if (it->first != "restart") {
          if (ecsPL.sublist(it->first).isParameter("start")) {
            value = ecsPL.sublist(it->first).get<std::string>("start");
            gotValue = true;
	  }
	}
      }
    }
    if (gotValue) {
      double time = get_time_value(value, *def_list);
      steadyPL.set<double>("Start",time);
      gotValue = false;
    }
    else {
      // ERROR - for unstructured, optional for structured;
    }
    if (defPL.isParameter("end")) {
      value = defPL.get<std::string>("end");
      gotValue = true;
    }
    else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	if (it->first != "restart") {
          if (ecsPL.sublist(it->first).isParameter("end")) {
            value = ecsPL.sublist(it->first).get<std::string>("end");
            gotValue = true;
	  }
	}
      }
    }
    if (gotValue) {
      double time = get_time_value(value, *def_list);
      steadyPL.set<double>("End",time);
      gotValue = false;
    }
    else {
      // ERROR ;
    }
    if (defPL.isParameter("init_dt")) {
      value = defPL.get<std::string>("init_dt");
      gotValue = true;
    }
    else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	if (it->first != "restart") {
          if (ecsPL.sublist(it->first).isParameter("init_dt")) {
            value = ecsPL.sublist(it->first).get<std::string>("init_dt");
            gotValue = true;
	  }
	}
      }
    }
    if (gotValue) {
      steadyPL.set<double>("Initial Time Step",get_time_value(value,*def_list));
      gotValue = false;
    }
    else {
      // default value to 0.0
      steadyPL.set<double>("Initial Time Step",0.0);
    }
    if (defPL.isParameter("method")) {
      method = defPL.get<std::string>("method");
      gotValue = true;
    }
    else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	if (it->first != "restart") {
          if (ecsPL.sublist(it->first).isParameter("method")) {
            method = ecsPL.sublist(it->first).get<std::string>("method");
            gotValue = true;
	  }
	}
      }
    }
    if (gotValue && strcmp(method.c_str(),"picard")==0) {
      /* EIB - moved 'Use Picard' as 1.2.2 update */
      //steadyPL.set<bool>("Use Picard","true");
      fpkPL.set<bool>("Use Picard","true");
      gotValue = false;
    }
    else {
      // ERROR ;
    }
    if (defPL.isParameter("reduction_factor")) {
      haveSSF = true;
    }
    else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	if (it->first != "restart") {
          if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
            value = ecsPL.sublist(it->first).get<std::string>("reduction_factor");
            haveSSF = true;
	  }
	}
      }
    }
    if (defPL.isParameter("increase_factor")) {
      haveSSF = true;
    }
    else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	if (it->first != "restart") {
          if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
            haveSSF = true;
	  }
	}
      }
    }
    list.sublist("Time Integration Mode").sublist("Steady") = steadyPL;

  }
  else {
    if (!hasSteady) {
      // Transient case
      if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
          *voI_->os() << "Creating Transient Execution Control" << std::endl;
      }
      
      Teuchos::ParameterList transPL;
      // loop over ecs to set up, TPC lists
      Teuchos::Array<double> start_times;
      Teuchos::Array<double> init_steps;
      Teuchos::Array<double> max_steps;
      int max_cycles;
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
        // skip this list, if labeled "restart" it is a duplicate list
        if (it->first != "restart") {
          bool gotValue;
          std::string Value;
          if (ecsPL.sublist(it->first).isParameter("start")) {
            Value = ecsPL.sublist(it->first).get<std::string>("start");
            gotValue = true;
          }
          else {
            Value = defPL.get<std::string>("start");
            gotValue = true;
          }
          if (gotValue) {
            double time = get_time_value(Value,*def_list);
            start_times.append(time);
            gotValue = false;
          }
          if (ecsPL.sublist(it->first).isParameter("end")) {
            Value = ecsPL.sublist(it->first).get<std::string>("end");
            gotValue = true;
          }
          else {
            Value = defPL.get<std::string>("end");
            gotValue = true;
          }
          if (gotValue) {
            double time = get_time_value(Value,*def_list);
            transPL.set<double>("End",time);
            gotValue = false;
          }
          if (ecsPL.sublist(it->first).isParameter("init_dt")) {
            Value = ecsPL.sublist(it->first).get<std::string>("init_dt");
            gotValue = true;
          }
          else {
            if (defPL.isParameter("init_dt")) {
              Value = defPL.get<std::string>("init_dt");
              gotValue = true;
            }
          }
          if (gotValue) {
            init_steps.append(get_time_value(Value,*def_list));
            gotValue = false;
          }
          if (ecsPL.sublist(it->first).isParameter("max_dt")) {
            Value = ecsPL.sublist(it->first).get<std::string>("max_dt");
            gotValue = true;
          }
          else {
            if (defPL.isParameter("max_dt")) {
              Value = defPL.get<std::string>("max_dt");
              gotValue = true;
            }
          }
          if (gotValue) {
            max_steps.append(get_time_value(Value,*def_list));
            gotValue = false;
          }
          if (ecsPL.sublist(it->first).isParameter("reduction_factor") || defPL.isParameter("reduction_factor")) {
            haveTF = true;
          }
          if (ecsPL.sublist(it->first).isParameter("increase_factor") || defPL.isParameter("increase_factor")) {
            haveTF = true;
          }
          
          /* EIB: proposed v1.2.2 update - Add Maximum Cycle Number for Transient modes */
          if (ecsPL.sublist(it->first).isParameter("max_cycles")) {
            Value = ecsPL.sublist(it->first).get<std::string>("max_cycles");
            transPL.set<int>("Maximum Cycle Number",get_int_constant(Value,*def_list));
          }
        }
      }
      transPL.set<double>("Start",start_times[0]);
      if ( init_steps.length() > 0 ) transPL.set<double>("Initial Time Step",init_steps[0]);
      if ( max_steps.length() > 0 ) transPL.set<double>("Maximum Time Step Size",max_steps[0]);
      if  (!staticflowON) {
        list.sublist("Time Integration Mode").sublist("Transient") = transPL;
      }
      else {
        list.sublist("Time Integration Mode").sublist("Transient with Static Flow") = transPL;
      }
      if (start_times.length() > 1) {
        // to include "Time Period Control" list
        Teuchos::ParameterList tpcPL;
        tpcPL.set<Teuchos::Array<double> >("Start Times",start_times);
        if ( init_steps.length() > 0 ) tpcPL.set<Teuchos::Array<double> >("Initial Time Step",init_steps);
        if ( max_steps.length() > 0 ) tpcPL.set<Teuchos::Array<double> >("Maximum Time Step",max_steps);
        list.sublist("Time Period Control") = tpcPL;
      }
    }
    else {
      // Initialize to Steady case
      if (numControlPeriods==1) {
      // user screwed up and really meant Steady case
      if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
          *voI_->os() << "Creating Steady State Execution Control" << std::endl;
      }
        
      Teuchos::ParameterList steadyPL;
      // look for values from default list
      // if not there, grab from ec list
      std::string value;
      std::string method;
      bool gotValue;
      if (defPL.isParameter("start")) {
        value = defPL.get<std::string>("start");
        gotValue = true;
      }
      else {
        for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	    if (it->first != "restart") {
              if (ecsPL.sublist(it->first).isParameter("start")) {
                value = ecsPL.sublist(it->first).get<std::string>("start");
                gotValue = true;
	      }
	    }
          }
      }
      if (gotValue) {
        double time = get_time_value(value, *def_list);
        steadyPL.set<double>("Start",time);
        gotValue = false;
      }
      else {
        // ERROR - for unstructured, optional for structured;
      }
      if (defPL.isParameter("end")) {
        value = defPL.get<std::string>("end");
        gotValue = true;
      }
      else {
        for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	    if (it->first != "restart") {
              if (ecsPL.sublist(it->first).isParameter("end")) {
                value = ecsPL.sublist(it->first).get<std::string>("end");
                gotValue = true;
	      }
	    }
          }
      }
      if (gotValue) {
        double time = get_time_value(value, *def_list);
        steadyPL.set<double>("End",time);
        gotValue = false;
      }
      else {
        // ERROR ;
      }
      if (defPL.isParameter("init_dt")) {
        value = defPL.get<std::string>("init_dt");
        gotValue = true;
      }
      else {
        for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	    if (it->first != "restart") {
              if (ecsPL.sublist(it->first).isParameter("init_dt")) {
                value = ecsPL.sublist(it->first).get<std::string>("init_dt");
                gotValue = true;
	      }
	    }
          }
      }
      if (gotValue) {
        steadyPL.set<double>("Initial Time Step",get_time_value(value,*def_list));
        gotValue = false;
      }
      else {
        // default value to 0.0
        steadyPL.set<double>("Initial Time Step",0.0);
      }
      if (defPL.isParameter("method")) {
        method = defPL.get<std::string>("method");
        gotValue = true;
      }
      else {
        for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
          if (it->first != "restart") {
            if (ecsPL.sublist(it->first).isParameter("method")) {
              method = ecsPL.sublist(it->first).get<std::string>("method");
              gotValue = true;
            }
          }
        }
      }
      if (gotValue && strcmp(method.c_str(),"picard")==0) {
        /* EIB - moved 'Use Picard' as 1.2.2 update */
        //steadyPL.set<bool>("Use Picard","true");
        fpkPL.set<bool>("Use Picard","true");
        gotValue = false;
      }
      else {
        // ERROR ;
      }
      if (defPL.isParameter("reduction_factor")) {
        haveSSF = true;
      }
      else {
        for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
          if (it->first != "restart") {
            if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
              value = ecsPL.sublist(it->first).get<std::string>("reduction_factor");
              haveSSF = true;
            }
          }
        }
      }
      if (defPL.isParameter("increase_factor")) {
        haveSSF = true;
      }
      else {
        for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	    if (it->first != "restart") {
              if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
                haveSSF = true;
    	      }
    	    }
          }
      }
      list.sublist("Time Integration Mode").sublist("Steady") = steadyPL;
    }
      else {
      // proceed, user really meant Initialize to Steady case
      if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
          *voI_->os() << "Creating Initialize to Steady Execution Control" << std::endl;
      }
      Teuchos::Array<double> start_times;
      Teuchos::Array<double> init_steps;
      Teuchos::Array<double> max_steps;
      Teuchos::ParameterList timesPL;
      Teuchos::ParameterList initPL;
      std::string value;
      bool gotValue = false;
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
        if (it->first != "restart") {
	    std::string mode("none");
	    if (defPL.isParameter("mode")) mode = defPL.get<std::string>("mode");
            if (ecsPL.sublist(it->first).isParameter("mode")) mode = ecsPL.sublist(it->first).get<std::string>("mode");
	    if (strcmp(mode.c_str(),"steady")==0) {
	      if (ecsPL.sublist(it->first).isParameter("start")) {
                value = ecsPL.sublist(it->first).get<std::string>("start");
	        double time = get_time_value(value,*def_list);
	        initPL.set<double>("Start",time);
	      }
	      if (ecsPL.sublist(it->first).isParameter("end")) {
                value = ecsPL.sublist(it->first).get<std::string>("end");
	        double time = get_time_value(value,*def_list);
	        initPL.set<double>("Switch",time);
	      }
	      if (ecsPL.sublist(it->first).isParameter("init_dt")) {
                value = ecsPL.sublist(it->first).get<std::string>("init_dt");
	        initPL.set<double>("Steady Initial Time Step",get_time_value(value,*def_list));
	      } 
	      if (ecsPL.sublist(it->first).isParameter("method")) {
                value = ecsPL.sublist(it->first).get<std::string>("method");
                if (strcmp(value.c_str(),"true")==0) {
                  fpkPL.set<bool>("Use Picard","true");
                }
                /* EIB - moved 'Use Picard' as 1.2.2 update */
	          //initPL.set<bool>("Use Picard",true);
      	      }
	      if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
                haveSSF = true;
	      } 
	      if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
                haveSSF = true;
	      } 
	    }
            else {
	      if (ecsPL.sublist(it->first).isParameter("start")) {
                value = ecsPL.sublist(it->first).get<std::string>("start");
	        double time = get_time_value(value,*def_list);
	        start_times.append(time);
	      }
	      if (ecsPL.sublist(it->first).isParameter("end")) {
                value = ecsPL.sublist(it->first).get<std::string>("end");
	        double time = get_time_value(value,*def_list);
	        initPL.set<double>("End",time);
	      }
	      if (ecsPL.sublist(it->first).isParameter("init_dt")) {
                value = ecsPL.sublist(it->first).get<std::string>("init_dt");
	        gotValue = true;
	      }
              else {
	        if (defPL.isParameter("init_dt")) {
                  value = defPL.get<std::string>("init_dt");
	          gotValue = true;
	        }
	      }
	      if (gotValue) {
	        init_steps.append(get_time_value(value,*def_list));
	        gotValue = false;
	      }
	    }
	  }
      }
      initPL.set<double>("Transient Initial Time Step",init_steps[0]);
      list.sublist("Time Integration Mode").sublist("Initialize To Steady") = initPL;
      if (start_times.length() > 1) {
        timesPL.set<Teuchos::Array<double> >("Start Times",start_times);
        timesPL.set<Teuchos::Array<double> >("Initial Time Step",init_steps);
        list.sublist("Time Period Control") = timesPL;
      }
    }
    }
  }

  // get numerical controls node
  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Numerical Controls" << std::endl;
  }
  
  Teuchos::ParameterList ncPL, commonPL, algoPL;
  // grab some variables saved in the defPL and ecsPL
  Teuchos::ParameterList tcPL, ssPL;
  // first check the default list
  if ( defPL.isParameter("mode") ) {
    if (defPL.get<std::string>("mode") == "transient") {
      if (defPL.isParameter("increase_factor")) {
        tcPL.set<double>("transient time step increase factor",
                         get_double_constant(defPL.get<std::string>("increase_factor"),*def_list));
      }
      if (defPL.isParameter("reduction_factor")) {
        tcPL.set<double>("transient time step reduction factor",
                         get_double_constant(defPL.get<std::string>("reduction_factor"),*def_list));
      }
      if (defPL.isParameter("max_dt")) {
        tcPL.set<double>("transient max time step",
                         get_time_value(defPL.get<std::string>("max_dt"),*def_list));
      }
    }
    if (defPL.get<std::string>("mode") == "steady") {
      if (defPL.isParameter("increase_factor")) {
        ssPL.set<double>("steady time step increase factor",
                         get_double_constant(defPL.get<std::string>("increase_factor"),*def_list));
      }
      if (defPL.isParameter("reduction_factor")) {
        ssPL.set<double>("steady time step reduction factor",
                         get_double_constant(defPL.get<std::string>("reduction_factor"),*def_list));
      }
      if (defPL.isParameter("max_dt")) {
        ssPL.set<double>("steady max time step",
                         get_time_value(defPL.get<std::string>("max_dt"),*def_list));
      }
    }
  }
  // next check the individual execution controls, this will overwrite default values
  for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
    if (ecsPL.sublist(it->first).isParameter("mode")) {
      if (ecsPL.sublist(it->first).get<std::string>("mode") == "transient") {
        if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
          tcPL.set<double>("transient time step increase factor",
                           get_double_constant(ecsPL.sublist(it->first).get<std::string>("increase_factor"),*def_list));
        }
        if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
          tcPL.set<double>("transient time step reduction factor",
                           get_double_constant(ecsPL.sublist(it->first).get<std::string>("reduction_factor"),*def_list));
        }
        if (ecsPL.sublist(it->first).isParameter("max_dt")) {
          tcPL.set<double>("transient max time step",
                           get_time_value(ecsPL.sublist(it->first).get<std::string>("max_dt"),*def_list));
        }
      }
      if (ecsPL.sublist(it->first).get<std::string>("mode") == "steady") {
        if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
          ssPL.set<double>("steady time step increase factor",
                           get_double_constant(ecsPL.sublist(it->first).get<std::string>("increase_factor"),*def_list));
        }
        if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
          ssPL.set<double>("steady time step reduction factor",
                           get_double_constant(ecsPL.sublist(it->first).get<std::string>("reduction_factor"),*def_list));
        }
        if (ecsPL.sublist(it->first).isParameter("max_dt")) {
          ssPL.set<double>("steady max time step",
                           get_time_value(ecsPL.sublist(it->first).get<std::string>("max_dt"),*def_list));
        }
      }
    }
  }
  
  // The following is for the Unstructured Algorithm ONLY
  DOMNodeList* topList = xmlDoc->getElementsByTagName(XMLString::transcode("numerical_controls"));
  std::string algo_str_name, not_algo_str_name;
  if (isUnstr_) {
    algo_str_name = "unstructured_controls";
    not_algo_str_name = "structured_controls";
  }
  else {
    algo_str_name = "structured_controls";
    not_algo_str_name = "unstructured_controls";
  }
  bool found = false;
  bool found_other = false;
  DOMNodeList* algoList;
  // loop over children to see if "unstructured_controls" exists; if so, process it
  for (int i=0; i<topList->getLength(); i++) {
    DOMNode* topNode = topList->item(i);
    if (DOMNode::ELEMENT_NODE == topNode->getNodeType()) {
      DOMNodeList* childList = topNode->getChildNodes();
      for(int j=0; j<childList->getLength(); j++) {
        DOMNode* tmpNode = childList->item(j) ;
        if (DOMNode::ELEMENT_NODE == tmpNode->getNodeType()) {
          char* algoName = XMLString::transcode(tmpNode->getNodeName());
          if (strcmp(algoName,algo_str_name.c_str())==0) {
            algoList = xmlDoc->getElementsByTagName(XMLString::transcode(algo_str_name.c_str()));
            found = true;
          }
          else if (strcmp(algoName,not_algo_str_name.c_str())==0) {
            found_other = true;
          }
        }
      }
    }
  }
  
  // have unstructured framework and found unstructured_controls
  if (isUnstr_ && found) {
    for (int i=0; i<algoList->getLength(); i++) {
      DOMNode* ncNode = algoList->item(i);
      if (DOMNode::ELEMENT_NODE == ncNode->getNodeType()) {
        DOMNodeList* childList = ncNode->getChildNodes();
        for(int j=0; j<childList->getLength(); j++) {
          DOMNode* tmpNode = childList->item(j) ;
          if (DOMNode::ELEMENT_NODE == tmpNode->getNodeType()) {
            std::string nodeName(std::string(XMLString::transcode(tmpNode->getNodeName())));
            if (nodeName == "unstr_flow_controls") {
              DOMNodeList* children = tmpNode->getChildNodes();
              for (int k=0; k<children->getLength(); k++) {
                DOMNode* currentNode = children->item(k) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                  std::string elemName(std::string(XMLString::transcode(currentNode->getNodeName())));
                  if (elemName == "discretization_method") {
                    std::string textValue(std::string(XMLString::transcode(currentNode->getTextContent())));
                    std::string discr_method;
                    if (textValue == "fv-default") {
                      discr_method.append("FV: Default");
                    }
                    else if (textValue == "fv-monotone") {
                      discr_method.append("FV: Monotone");
                    }
                    
                    else if (textValue == "fv-multi_point_flux_approximation") {
                      discr_method.append("FV: Multi-Point Flux Approximation");
                    }
                    
                    else if (textValue == "fv-extended_to_boundary_edges") {
                      discr_method.append("FV: Extended to Boundary Edges");
                    }
                    else if (textValue == "mfd-default") {
                      discr_method.append("MFD: Default");
                    }
                    else if (textValue == "mfd-optimized_for_sparsity") {
                      discr_method.append("MFD: Optimized for Sparsity");
                    }
                    else if (textValue == "mfd-support_operator") {
                      discr_method.append("MFD: Support Operator");
                    }
                    else if (textValue == "mfd-optimized_for_monotonicity") {
                      discr_method.append("MFD: Optimized for Monotonicity");
                    }
                    else if (textValue == "mfd-two_point_flux_approximation") {
                      discr_method.append("MFD: Two-Point Flux Approximation");
                    }
                    fpkPL.set<std::string>("Discretization Method",discr_method.c_str());
                    XMLString::release(&textContent);
                  }
                  else if (elemName == "rel_perm_method") {
                    std::string textValue(std::string(XMLString::transcode(currentNode->getTextContent())));
                    std::string rel_perm_method("upwind-darcy_velocity");
                    if (textValue == "upwind-gravity") {
                      rel_perm_method = "Upwind: Gravity";
                    }
                    else if (textValue == "upwind-darcy_velocity") {
                      rel_perm_method = "Upwind: Darcy Velocity";
                    }
                    else if (textValue == "upwind-amanzi") {
                      rel_perm_method = "Upwind: Amanzi";
                    }
                    else if (textValue == "other-arithmetic_average") {
                      rel_perm_method = "Other: Arithmetic Average";
                    }
                    
                    else if (textValue == "other-harmonic_average") {
                      rel_perm_method = "Other: Harmonic Average";
                    }
                    fpkPL.set<std::string>("Relative Permeability",rel_perm_method.c_str());
                    XMLString::release(&textContent);

                  }
                  else if (elemName == "atmospheric_pressure") {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    fpkPL.set<double>("atmospheric pressure",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (elemName == "preconditioning_strategy") {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    std::string textValue(std::string(textContent));
                    std::string op("linearized_operator");
                    if (strcmp(textContent,"diffusion_operator")==0) {
                      op = "Diffusion Operator";
                    }
                    else if (strcmp(textContent,"linearized_operator")==0) {
                      op = "Linearized Operator";
                    }
                    fpkPL.set<std::string>("Preconditioning Strategy",op.c_str());
                    XMLString::release(&textContent);
                  }
                  else if (elemName != "comments") {
                    // warn about unrecognized element
                    std::stringstream elem_name;
                    elem_name << nodeName <<  "->"  << elemName;
                    throw_warning_skip(elem_name.str());
                  }
                }
              }
            }
            else if (nodeName == "unstr_transport_controls") {
              DOMNodeList* children = tmpNode->getChildNodes();
              for (int k=0; k<children->getLength(); k++) {
                DOMNode* currentNode = children->item(k) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                  std::string elemName(std::string(XMLString::transcode(currentNode->getNodeName())));
                  if (elemName == "algorithm") {
                    std::string textValue(std::string(XMLString::transcode(currentNode->getTextContent())));
                    std::string algo("explicit first-order");
                    if (textValue == "explicit first-order") {
                      algo = "Explicit First-Order";
                    }
                    else if (textValue == "explicit second-order") {
                      algo = "Explicit Second-Order";
                    }
                    else if (textValue == "none") {
                      algo = "None";
                    }
                    tpkPL.set<std::string>("Transport Integration Algorithm",algo.c_str());
                  }
                  else if (elemName == "sub_cycling") {
                    std::string textValue(std::string(XMLString::transcode(currentNode->getTextContent())));
                    tpkPL.set<bool>("transport subcycling",false);
                    if (textValue == "on") {
                      tpkPL.set<bool>("transport subcycling",true);
                    }
                    else {
                      tpkPL.set<bool>("transport subcycling",false);
                    }
                  }
                  else if (elemName == "cfl") {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    tpkPL.set<double>("CFL",atof(textContent));
                    XMLString::release(&textContent);
                  }
                  else if (elemName != "comments") {
                    // warn about unrecognized element
                    std::stringstream elem_name;
                    elem_name << nodeName <<  "->"  << elemName;
                    throw_warning_skip(elem_name.str());
                  }
                }
              }
            }
            else if (nodeName == "unstr_steady-state_controls") {
            // loop through children and deal with them
            DOMNodeList* children = tmpNode->getChildNodes();
            for (int k=0; k<children->getLength(); k++) {
              DOMNode* currentNode = children->item(k) ;
              if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
    	        char* tagname = XMLString::transcode(currentNode->getNodeName());
                if (strcmp(tagname,"min_iterations")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<int>("steady min iterations",get_int_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"max_iterations")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<int>("steady max iterations",get_int_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"limit_iterations")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<int>("steady limit iterations",get_int_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"max_preconditioner_lag_iterations")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<int>("steady max preconditioner lag iterations",get_int_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"nonlinear_tolerance")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<double>("steady nonlinear tolerance",get_double_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"nonlinear_iteration_divergence_factor")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<double>("steady nonlinear iteration divergence factor",get_double_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"nonlinear_iteration_damping_factor")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<double>("steady nonlinear iteration damping factor",get_double_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"max_divergent_iterations")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<int>("steady max divergent iterations",get_int_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"error_control_options")==0) { //default = 'pressure' options = 'pressure','residual'
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  Teuchos::Array<std::string> err_opts = make_regions_list(textContent);
                  ssPL.set<Teuchos::Array<std::string> >("steady error control options",err_opts);
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"error_rel_tol")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<double>("steady error rel tol",get_double_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"error_abs_tol")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<double>("steady error abs tol",get_double_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"restart_tolerance_relaxation_factor")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<double>("steady restart tolerance relaxation factor",get_double_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"restart_tolerance_relaxation_factor_damping")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<double>("steady restart tolerance relaxation factor damping",get_double_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"initialize_with_darcy")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  bool iwd(false);
                  std::string(textContent) == "true" ? iwd = true : iwd = false;
                  if (!iwd)
                    std::string(textContent) == "1" ? iwd = true : iwd = false;
                  ssPL.set<bool>("steady initialize with darcy",iwd);
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"preconditioner")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  
                  if (strcmp(textContent,"hypre_amg")==0) {
                    ssPL.set<std::string>("steady preconditioner","Hypre AMG");
                  }
                  else if (strcmp(textContent,"trilinos_ml")==0) {
                    ssPL.set<std::string>("steady preconditioner","Trilinos ML");
                  }
                  else if (strcmp(textContent,"block_ilu")==0) {
                    ssPL.set<std::string>("steady preconditioner","Block ILU");
                  }
                  else {
                    throw_error_illformed("steady preconditioner", "value", "preconditioner", "trilinos_ml, hypre_amg, block_ilu");
                  }
                  
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"nonlinear_iteration_initial_guess_extrapolation_order")==0) {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  ssPL.set<int>("steady nonlinear iteration initial guess extrapolation order",get_int_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (strcmp(tagname,"unstr_initialization")==0) {
                  Teuchos::ParameterList ptiPL;
                  DOMNodeList* kids = currentNode->getChildNodes();
                  for (int l=0; l<kids->getLength(); l++) {
                    DOMNode* curNode = kids->item(l) ;
                    if (DOMNode::ELEMENT_NODE == curNode->getNodeType()) {
                      std::string tag = std::string(XMLString::transcode(curNode->getNodeName()));
                      if (tag == "method") {
                        textContent = XMLString::transcode(curNode->getTextContent());
                        if (strcmp(textContent,"picard")==0) {
                          ptiPL.set<std::string>("time integration method","Picard");
                        }
                        else if (strcmp(textContent,"darcy_solver")==0) {
                          ptiPL.set<std::string>("time integration method","darcy solver");
                        }
                        else {
                          //TODO: EIB - warn, unrecognized method name
                        }
                        XMLString::release(&textContent);
                      }
                      else if (tag == "preconditioner") {
                        textContent = XMLString::transcode(curNode->getTextContent());
                        if (strcmp(textContent,"trilinos_ml")==0) {
                          ptiPL.set<std::string>("preconditioner","Trilinos ML");
                        }
                        else if (strcmp(textContent,"hypre_amg")==0) {
                          ptiPL.set<std::string>("preconditioner","Hypre AMG");
                        }
                        else if (strcmp(textContent,"block_ilu")==0) {
                          ptiPL.set<std::string>("preconditioner","Block ILU");
                        }
                        else {
                          throw_error_illformed("unstr_initialization", "value", "preconditioner", "trilinos_ml, hypre_amg, block_ilu");
                        }
                        XMLString::release(&textContent);
                      }
                      else if (tag == "linear_solver") {
                        textContent = XMLString::transcode(curNode->getTextContent());
                        if (strcmp(textContent,"aztec00")==0) {
                          ptiPL.set<std::string>("linear solver","AztecOO");
                        }
                        XMLString::release(&textContent);
                      }
                      else if (tag == "error_control_options") { //default = 'pressure'
                        textContent = XMLString::transcode(curNode->getTextContent());
                        Teuchos::Array<std::string> err_opts = make_regions_list(textContent);
                        ptiPL.set<Teuchos::Array<std::string> >("error control options",err_opts);
                        XMLString::release(&textContent);
                      }
                      else if (tag == "max_iterations") {
                        textContent = XMLString::transcode(curNode->getTextContent());
                        ptiPL.set<int>("picard maximum number of iterations",
                                       get_int_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (tag == "clipping_saturation") {
                        textContent = XMLString::transcode(curNode->getTextContent());
                        ptiPL.set<double>("clipping saturation value",
                                          get_double_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (tag == "clipping_pressure") {
                        textContent = XMLString::transcode(curNode->getTextContent());
                        ptiPL.set<double>("clipping pressure value",
                                          get_double_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (tag == "convergence_tolerance") {
                        textContent = XMLString::transcode(curNode->getTextContent());
                        ptiPL.set<double>("picard convergence tolerance",
                                          get_double_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                    }
                  }
                  //ssPL.sublist("Steady-State Pseudo-Time Implicit Solver") = ptiPL;
                  list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Initialization") = ptiPL;
                }
                else if (strcmp(tagname,"comments")!=0) {
                  // warn about unrecognized element
                  std::stringstream elem_name;
                  elem_name << nodeName <<  "->"  << tagname;
                  throw_warning_skip(elem_name.str());
                }
              }
	    }
	    // get items from steadyPL, from execution_controls
	    if (ecsPL.isSublist("steady")) {
	      if (ecsPL.sublist("steady").isParameter("max_dt")) {
		ssPL.set<double>("steady max time step",get_time_value(ecsPL.sublist("steady").get<std::string>("max_dt"),*def_list));
	      }
	    }

            //list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Steady-State Implicit Time Integration") = ssPL;
	  }
            else if (nodeName == "unstr_transient_controls") {
	    // check for incr/red factors from execution_controls first
	    // grab integration method, then loop through it's attributes
            DOMElement* tcElement = static_cast<DOMElement*>(tmpNode);
            
            DOMNodeList* children = tmpNode->getChildNodes();
            for (int k=0; k<children->getLength(); k++) {
              DOMNode* currentNode = children->item(k) ;
              if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                std::string tagname(XMLString::transcode(currentNode->getNodeName()));
                
                if (tagname == "bdf1_integration_method") {
                  DOMNodeList* gkids = currentNode->getChildNodes();
                  for (int l=0; l<gkids->getLength(); l++) {
                    DOMNode* kidNode = gkids->item(l) ;
                    if (DOMNode::ELEMENT_NODE == kidNode->getNodeType()) {
                      std::string bdf1_tagname(XMLString::transcode(kidNode->getNodeName()));
                      if (bdf1_tagname == "min_iterations") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<int>("transient min iterations",get_int_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "max_iterations") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<int>("transient max iterations",get_int_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "limit_iterations") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<int>("transient limit iterations",get_int_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "nonlinear_tolerance") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<double>("transient nonlinear tolerance",get_double_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "max_divergent_iterations") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<int>("transient max divergent iterations",get_int_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "max_preconditioner_lag_iterations") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<int>("transient max preconditioner lag iterations",get_int_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "nonlinear_iteration_damping_factor") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<double>("transient nonlinear iteration damping factor",get_double_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "nonlinear_iteration_divergence_factor") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<double>("transient nonlinear iteration divergence factor",get_double_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "restart_tolerance_relaxation_factor") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<double>("transient restart tolerance relaxation factor",get_double_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "restart_tolerance_relaxation_factor_damping") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        tcPL.set<double>("transient restart tolerance relaxation factor damping",get_double_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "initialize_with_darcy") {
                        textContent = XMLString::transcode(kidNode->getTextContent());bool iwd(false);
                        std::string(textContent) == "true" ? iwd = true : iwd = false;
                        if (!iwd)
                          std::string(textContent) == "1" ? iwd = true : iwd = false;
                        tcPL.set<bool>("transient initialize with darcy",iwd);
                        XMLString::release(&textContent);
                      }
                      else if (bdf1_tagname == "error_control_options") {
                        textContent = XMLString::transcode(kidNode->getTextContent());
                        Teuchos::Array<std::string> err_opts = make_regions_list(textContent);
                        tcPL.set<Teuchos::Array<std::string> >("transient error control options",err_opts);
                        XMLString::release(&textContent);
                      }
                      else {
                        if (bdf1_tagname != "comments") {
                          // warn about unrecognized element
                          std::stringstream elem_name;
                          elem_name << nodeName <<  "->"  << bdf1_tagname;
                          throw_warning_skip(elem_name.str());
                        }
                      }
                    }
                  }
                }
                else if (tagname == "preconditioner") {
                  std::string value(trim_string(XMLString::transcode(currentNode->getTextContent())));
                  if (value == "trilinos_ml") {
                    tcPL.set<std::string>("transient preconditioner","Trilinos ML");
                  }
                  else if (value == "hypre_amg") {
                    tcPL.set<std::string>("transient preconditioner","Hypre AMG");
                  }
                  else if (value == "block_ilu") {
                    tcPL.set<std::string>("transient preconditioner","Block ILU");
                  }
                  else {
                    throw_error_illformed("transient preconditioner", "value", "preconditioner", "trilinos_ml, hypre_amg, block_ilu");
                  }
                }
                else if (tagname == "initialize_with_darcy") {
                  std::string value(trim_string(XMLString::transcode(currentNode->getTextContent())));
                  bool iwd(false);
                  value == "true" ? iwd = true : iwd = false;
                  if (!iwd)
                    value == "1" ? iwd = true : iwd = false;
                  tcPL.set<bool>("transient initialize with darcy",iwd);
                  XMLString::release(&textContent);
                }
                else if (tagname == "nonlinear_iteration_initial_guess_extrapolation_order") {
                  textContent = XMLString::transcode(currentNode->getTextContent());
                  tcPL.set<int>("transient nonlinear iteration initial guess extrapolation order",get_int_constant(textContent,*def_list));
                  XMLString::release(&textContent);
                }
                else if (tagname != "comments") {
                  // warn about unrecognized element
                  std::stringstream elem_name;
                  elem_name << nodeName <<  "->"  << tagname;
                  throw_warning_skip(elem_name.str());
                }
              }
            }
          }
            else if (nodeName == "unstr_nonlinear_solver") {
              Teuchos::ParameterList nlPL;
              attrMap = tmpNode->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
              if (nodeAttr) {
                textContent = XMLString::transcode(nodeAttr->getNodeValue());
              }
              else {
                throw_error_missattr("numerical_controls","attribute","name","nonlinear_solver");
              }

              if (strcmp(textContent,"nka")==0) {
                nlPL.set<std::string>("Nonlinear Solver Type","NKA");
              }
              else if (strcmp(textContent,"newton")==0) {
                nlPL.set<std::string>("Nonlinear Solver Type","Newton");
              }
              else if (strcmp(textContent,"jfnk")==0) {
                nlPL.set<std::string>("Nonlinear Solver Type","JFNK");
              }
              else if (strcmp(textContent,"newton_picard")==0) {
                nlPL.set<std::string>("Nonlinear Solver Type","Newton-Picard");
              }
              XMLString::release(&textContent);
              
              // loop through children and deal with them
              DOMNodeList* children = tmpNode->getChildNodes();
              for (int k=0; k<children->getLength(); k++) {
                DOMNode* currentNode = children->item(k) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                  std::string tagname(std::string(XMLString::transcode(currentNode->getNodeName())));
                  if (tagname == "modify_correction") {
                    std::string textContent(trim_string(XMLString::transcode(currentNode->getTextContent())));
                    bool iwd(false);
                    textContent == "true" ? iwd = true : iwd = false;
                    if (!iwd)
                      textContent == "1" ? iwd = true : iwd = false;
                    nlPL.set<bool>("modify correction",iwd);
                  }
                  else if (tagname == "update_upwind_frequency") {
                    // which precondition is stored in attribute, options are: trilinos_ml, hypre_amg, block_ilu
                    std::string textContent(std::string(XMLString::transcode(currentNode->getTextContent())));
                    //usePCPL = true;
                    
                    if (textContent == "every_timestep") {
                      nlPL.set<std::string>("update upwind frequency","every timestep");
                    }
                    else if (textContent == "every_nonlinear_iteration") {
                      nlPL.set<std::string>("update upwind frequency","every nonlinear iteration");
                    }
                    else {
                      throw_error_illformed(nodeName, "value", "update_upwind_frequency", "every_timestep, every_nonlinear_iteration");
                    }
                  }
                  else if (tagname == "comments") {
                    // warn about unrecognized element
                    std::stringstream elem_name;
                    elem_name << nodeName <<  "->"  << tagname;
                    throw_warning_skip(elem_name.str());
                  }
                }
              }
              list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Nonlinear Solver") = nlPL;
            }
            else if (nodeName == "unstr_linear_solver") {
              Teuchos::ParameterList lsPL;
              Teuchos::ParameterList pcPL;
              bool usePCPL=false;
              // loop through children and deal with them
              DOMNodeList* children = tmpNode->getChildNodes();
              for (int k=0; k<children->getLength(); k++) {
                DOMNode* currentNode = children->item(k) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                  char* tagname = XMLString::transcode(currentNode->getNodeName());
                  if (strcmp(tagname,"method")==0) {
                    std::string value(trim_string(XMLString::transcode(currentNode->getTextContent())));
                    if (value.length() > 0) {
                    lsPL.set<std::string>("linear solver iterative method",value);
                    }
                    else { // set default value
                      lsPL.set<std::string>("linear solver iterative method","gmres");
                    }
                  }
                  else if (strcmp(tagname,"max_iterations")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    lsPL.set<int>("linear solver maximum iterations",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"tolerance")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    lsPL.set<double>("linear solver tolerance",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"cfl")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    tpkPL.set<double>("CFL",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"preconditioner")==0) {
                    // which precondition is stored in attribute, options are: trilinos_ml, hypre_amg, block_ilu
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    //usePCPL = true;
                    
                    if (strcmp(textContent,"hypre_amg")==0) {
                      lsPL.set<std::string>("linear solver preconditioner","Hypre AMG");
                    }
                    else if (strcmp(textContent,"trilinos_ml")==0) {
                      lsPL.set<std::string>("linear solver preconditioner","Trilinos ML");
		    }
                    else if (strcmp(textContent,"block_ilu")==0) {
                      lsPL.set<std::string>("linear solver preconditioner","Block ILU");
                    }
                    else {
                      throw_error_illformed("linear solver preconditioner", "value", "preconditioner", "trilinos_ml, hypre_amg, block_ilu");
                    }
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"comments")!=0) {
                    // warn about unrecognized element
                    std::stringstream elem_name;
                    elem_name << nodeName <<  "->"  << tagname;
                    throw_warning_skip(elem_name.str());
                  }
                }
              }
              list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Linear Solver") = lsPL;
              //if (usePCPL)
              //  list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Preconditioners") = pcPL;
            }
            else if (nodeName == "nonlinear_solver") {
            // EIB: creating sub for section that doesn't actually exist yet in the New Schema, but does in the Input Spec
            Teuchos::ParameterList nlsPL;
            DOMNodeList* children = tmpNode->getChildNodes();
            for (int k=0; k<children->getLength(); k++) {
                DOMNode* currentNode = children->item(k) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
    	          char* tagname = XMLString::transcode(currentNode->getNodeName());
                  if (strcmp(tagname,"nonlinear_solver_type")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    nlsPL.set<std::string>("Nonlinear Solver Type",trim_string(textContent));
                    XMLString::release(&textContent);
		  }
                }
	      }
              list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Nonlinear Solver") = nlsPL;
            }
            else if (nodeName == "unstr_preconditioners") {
              // loop over known preconditioner types
              DOMNodeList* children = tmpNode->getChildNodes();
              for (int k=0; k<children->getLength(); k++) {
                DOMNode* currentNode = children->item(k) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                  std::string elemName(std::string(XMLString::transcode(currentNode->getNodeName())));
                  if (elemName == "hypre_amg") {
                    // loop over children to get options
                    Teuchos::ParameterList hypre;
                    DOMNodeList* gkids = currentNode->getChildNodes();
                    for (int l=0; l<gkids->getLength(); l++) {
                      DOMNode* currentGKid = gkids->item(l) ;
                      if (DOMNode::ELEMENT_NODE == currentGKid->getNodeType()) {
                        std::string gkidName(std::string(XMLString::transcode(currentGKid->getNodeName())));
                        if (gkidName == "hypre_cycle_applications") {
                          hypre.set<int>("Hypre AMG cycle applications",
                                         get_int_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                        else if (gkidName == "hypre_smoother_sweeps") {
                          hypre.set<int>("Hypre AMG smoother sweeps",
                                         get_int_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                        else if (gkidName == "hypre_tolerance") {
                          hypre.set<double>("Hypre AMG tolerance",
                                            get_double_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                        else if (gkidName == "hypre_strong_threshold") {
                          hypre.set<double>("Hypre AMG strong threshold",
                                            get_double_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                      }
                    }
                    list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Preconditioners").sublist("Hypre AMG") = hypre;
                  }
                  else if (elemName == "trilinos_ml") {
                    // loop over children to get options
                    Teuchos::ParameterList trilinos;
                    DOMNodeList* gkids = currentNode->getChildNodes();
                    for (int l=0; l<gkids->getLength(); l++) {
                      DOMNode* currentGKid = gkids->item(l) ;
                      if (DOMNode::ELEMENT_NODE == currentGKid->getNodeType()) {
                        std::string gkidName(std::string(XMLString::transcode(currentGKid->getNodeName())));
                        if (gkidName == "trilinos_smoother_type") {
                          std::string stype(XMLString::transcode(currentNode->getTextContent()));
                          if (stype == "jacobi") {
                            trilinos.set<std::string>("ML smoother type","Jacobi");
                          }
                          else if (stype == "gauss_seidel") {
                            trilinos.set<std::string>("ML smoother type","Gauss-Seidel");
                          }
                          else if (stype == "ilu") {
                            trilinos.set<std::string>("ML smoother type","ILU");
                          }
                        
                        }
                        else if (gkidName == "trilinos_threshold") {
                          trilinos.set<double>("ML aggregation threshold",
                                         get_double_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                        else if (gkidName == "trilinos_smoother_sweeps") {
                          trilinos.set<int>("ML smoother sweeps",
                                            get_int_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                        else if (gkidName == "trilinos_cycle_applications") {
                          trilinos.set<int>("ML cycle applications",
                                            get_int_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                      }
                    }
                    list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Preconditioners").sublist("Trilinos ML") = trilinos;
                  }
                  else if (elemName == "block_ilu") {
                    // loop over children to get options
                    Teuchos::ParameterList ilu;
                    DOMNodeList* gkids = currentNode->getChildNodes();
                    for (int l=0; l<gkids->getLength(); l++) {
                      DOMNode* currentGKid = gkids->item(l) ;
                      if (DOMNode::ELEMENT_NODE == currentGKid->getNodeType()) {
                        std::string gkidName(std::string(XMLString::transcode(currentGKid->getNodeName())));
                        if (gkidName == "ilu_overlap") {
                          ilu.set<int>("Block ILU overlap",
                                         get_int_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                        else if (gkidName == "ilu_relax") {
                          ilu.set<double>("Block ILU relax value",
                                         get_double_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                        else if (gkidName == "ilu_rel_threshold") {
                          ilu.set<double>("Block ILU relative threshold",
                                            get_double_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                        else if (gkidName == "ilu_abs_threshold") {
                          ilu.set<double>("Block ILU absolute threshold",
                                            get_double_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                        else if (gkidName == "ilu_level_of_fill") {
                          ilu.set<int>("Block ILU level of fill",
                                       get_int_constant(XMLString::transcode(currentGKid->getTextContent()),*def_list));
                        }
                      }
                    }
                    list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Preconditioners").sublist("Block ILU") = ilu;
                  }
                  else {
                    throw_error_illformed("preconditioners", "subelement", "preconditioner type", "trilinos_ml, hypre_amg, block_ilu");
                  }
                }
              }
              
            }
            else if (nodeName == "unstr_chemistry_controls") {
              Teuchos::ParameterList chemistryPL;
              // chemistry options are difference based on engine
              // deal with common items
              Teuchos::Array<std::string> verb;
              if (voI_->getVerbLevel() == Teuchos::VERB_EXTREME) {
                  verb.append("error");
                  chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
              }
              else if (voI_->getVerbLevel() == Teuchos::VERB_HIGH) {
                  verb.append("warning");
                  chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
              }
              else if (voI_->getVerbLevel() == Teuchos::VERB_MEDIUM) {
                  verb.append("verbose");
                  chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
              }
              else if (voI_->getVerbLevel() == Teuchos::VERB_LOW ) {
                  verb.append("terse");
                  chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
              }
              else {
                  verb.append("silent");
                  chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
              }
              
              // determine engine and proceed accordingly
              if (def_list->isParameter("chemistry_engine")) {
                if (def_list->get<std::string>("chemistry_engine") == "amanzi") {
                  // go ahead and add bdg file to PL
                  // build bgd filename
                  std::string bgdfilename;
                  if (def_list->isParameter("xmlfilename") ) {
                    bgdfilename = def_list->get<std::string>("xmlfilename") ;
                    std::string new_extension(".bgd");
                    size_t pos = bgdfilename.find(".xml");
                    bgdfilename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)4);
                  }
                  else {
                    // defaulting to hardcoded name
                    bgdfilename = "isotherms.bgd" ;
                  }
                  // add bgd file and parameters to list
                  Teuchos::ParameterList bgdPL;
                  bgdPL.set<std::string>("Format","simple");
                  bgdPL.set<std::string>("File",bgdfilename);
                  chemistryPL.sublist("Thermodynamic Database") = bgdPL;
                  chemistryPL.set<std::string>("Activity Model","unit");
                  
                  // loop over chemistry controls to get other options to add to PL
                  DOMNodeList* children = tmpNode->getChildNodes();
                  for (int k=0; k<children->getLength(); k++) {
                    DOMNode* currentNode = children->item(k) ;
                    if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                      std::string tagname(std::string(XMLString::transcode(currentNode->getNodeName())));
                      if (tagname == "chem_tolerance") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("Tolerance",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "chem_tolerance", nodeName);
                        }
                      }
                      else if (tagname == "chem_max_newton_iterations") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<int>("Maximum Newton Iterations",get_int_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "chem_max_newton_iterations", nodeName);
                        }
                        // TODO:: EIB - this need to be added to schema!!
                      }
                      else if (tagname == "chem_max_time_step") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("Max Time Step (s)",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "chem_max_time_step", nodeName);
                        }
                      }
                      else if (tagname == "max_chemistry_transport_timestep_ratio") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("max chemistry to transport timestep ratio",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "max_chemistry_transport_timestep_ratio", nodeName);
                        }
                      }
                      else if (tagname == "process_model") {
                        // TODO: EIB - not sure where this goes anymore
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("Tolerance",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        // TODO: EIB - removed error message until I figure out where this went
                        /*
                         else {
                         throw_error_illformed(algo_str_name, "process_model", nodeName);
                         }
                         */
                      }
                      else if (tagname != "comments") {
                        // warn about unrecognized element
                        std::stringstream elem_name;
                        elem_name << nodeName <<  "->"  << tagname;
                        throw_warning_skip(elem_name.str());
                      }
                    }
                  }
                }
                else if (def_list->get<std::string>("chemistry_engine") == "pflotran") {
                  chemistryPL.set<std::string>("Engine","PFloTran");
                  
                  // check file *.in filename in def_list
                  if (def_list->isSublist("chemistry_PL")) {
                    if (def_list->sublist("chemistry_PL").isParameter("Engine Input File")) {
                      chemistryPL.set<std::string>("Engine Input File",def_list->sublist("chemistry_PL").get<std::string>("Engine Input File"));
                    }
                  }
                  
                  // loop over children
                  DOMNodeList* children = tmpNode->getChildNodes();
                  for (int k=0; k<children->getLength(); k++) {
                    DOMNode* currentNode = children->item(k) ;
                    if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                      std::string tagname(std::string(XMLString::transcode(currentNode->getNodeName())));
                      if (tagname == "chem_max_time_step") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("Max Time Step (s)",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "chem_max_time_step", nodeName);
                        }
                      }
                      else if (tagname == "generate_chemistry_engine_inputfile") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<std::string>("Engine Input File",textContent);
                          XMLString::release(&textContent);
                        }
                      }
                      else if (tagname == "read_chemistry_engine_inputfile") {
                        // TODO: EIB - not sure where this goes anymore
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<std::string>("Engine Input File",textContent);
                          XMLString::release(&textContent);
                        }
                      }
                      else if (tagname != "comments") {
                        // warn about unrecognized element
                        std::stringstream elem_name;
                        elem_name << nodeName <<  "->"  << tagname;
                        throw_warning_skip(elem_name.str());
                      }
                    }
                  }
                }
              }
            
              
              // now add chemistry list
              def_list->sublist("Chemistry") = chemistryPL;
            }
          }
        }
      }
    }
    if (flowON)
      list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Flow Process Kernel") = fpkPL;
    if (transportON)
      list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Transport Process Kernel") = tpkPL;
    if (chemistryON)
      list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Chemistry Process Kernel") = cpkPL;
    list.sublist("Numerical Control Parameters").sublist("Common Controls") = commonPL;
    list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Steady-State Implicit Time Integration") = ssPL;
    list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Transient Implicit Time Integration") = tcPL;
  }
  if (isUnstr_ && !found) {
    Errors::Message msg;
    msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing numerical_controls - " ;
    msg << "No element 'unstructured_controls' found. \n" ;
    msg << "The element is required for unstructured framework. \n";
    msg << "  Please correct and try again \n" ;
    Exceptions::amanzi_throw(msg);
  }
  if (isUnstr_ && found_other) {
    //TODO: EIB - ignoring message
  }
  
  // have structured framework and found structured_controls
  Teuchos::ParameterList expertPL, amrPL;
  if (!isUnstr_ && found) {
    for (int i=0; i<algoList->getLength(); i++) {
      DOMNode* ncNode = algoList->item(i);
      if (DOMNode::ELEMENT_NODE == ncNode->getNodeType()) {
        DOMNodeList* childList = ncNode->getChildNodes();
        for(int j=0; j<childList->getLength(); j++) {
          DOMNode* tmpNode = childList->item(j) ;
          if (DOMNode::ELEMENT_NODE == tmpNode->getNodeType()) {
            char* nodeName = XMLString::transcode(tmpNode->getNodeName());
            if (strcmp(nodeName,"str_steady-state_controls")==0) {
              // loop through children and deal with them
              DOMNodeList* children = tmpNode->getChildNodes();
              for (int k=0; k<children->getLength(); k++) {
                DOMNode* currentNode = children->item(k) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                  char* tagname = XMLString::transcode(currentNode->getNodeName());
                  if (strcmp(tagname,"max_pseudo_time")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("steady_max_psuedo_time",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"min_iterations")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("steady_min_iterations",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"limit_iterations")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("steady_limit_iterations",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"max_iterations")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("steady_max_iterations",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"min_iterations_2")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("steady_min_iterations_2",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"time_step_increase_factor_2")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("steady_time_step_increase_factor_2",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"time_step_increase_factor")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("steady_time_step_increase_factor",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"time_step_reduction_factor")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("steady_time_step_reduction_factor",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"max_consecutive_failures_1")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("steady_max_consecutive_failures_1",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"time_step_retry_factor_1")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("steady_time_step_retry_factor_1",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"max_consecutive_failures_2")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("steady_max_consecutive_failures_2",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"time_step_retry_factor_2")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("steady_time_step_retry_factor_2",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"time_step_retry_factor_f")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("steady_time_step_retry_factor_f",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"max_num_consecutive_success")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("steady_max_num_consecutive_success",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"extra_time_step_increase_factor")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("steady_extra_time_step_increase_factor",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"abort_on_pseudo_timestep_failure")==0) {
		    textContent = XMLString::transcode(currentNode->getTextContent());
                    if (strcmp(textContent,"true")==0) {
                      expertPL.set<bool>("steady_abort_on_pseudo_timestep_failure",true);
                    }
                    else  {
                      expertPL.set<bool>("steady_abort_on_pseudo_timestep_failure",false);
                    }
                    XMLString::release(&textContent);
                    //textContent = XMLString::transcode(currentNode->getTextContent());
                    //expertPL.set<int>("steady_abort_on_pseudo_timestep_failure",get_int_constant(textContent,*def_list));
                    //XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"limit_function_evals")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("steady_limit_function_evals",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"do_grid_sequence")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    if (strcmp(textContent,"true")==0) {
                      expertPL.set<bool>("steady_do_grid_sequence",true);
                    }
                    else  {
                      expertPL.set<bool>("steady_do_grid_sequence",false);
                    }
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"grid_sequence_new_level_dt_factor")==0) {
                    Teuchos::Array<double> factors;
                    DOMNodeList* intChildren = currentNode->getChildNodes();
                    for (int l=0; l<intChildren->getLength(); l++) {
                      DOMNode* currentIntKid = intChildren->item(l) ;
                      if (DOMNode::ELEMENT_NODE == currentIntKid->getNodeType()) {
                        char* textContent = XMLString::transcode(currentIntKid->getTextContent());
                        factors.append(get_double_constant(textContent,*def_list));
                        XMLString::release(&textContent);
                      }
                    }
                    expertPL.set<Teuchos::Array<double> >("steady_grid_sequence_new_level_dt_factor",factors);
                    //XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"dt_thresh_pure_steady")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("richard_dt_thresh_pure_steady",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
		  XMLString::release(&tagname);
                }
              }
            }
            else if (strcmp(nodeName,"str_transient_controls")==0) {
              // loop through children and deal with them
              DOMNodeList* children = tmpNode->getChildNodes();
              for (int k=0; k<children->getLength(); k++) {
                DOMNode* currentNode = children->item(k) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                  char* tagname = XMLString::transcode(currentNode->getNodeName());
                  if (strcmp(tagname,"max_ls_iterations")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("richard_max_ls_iterations",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"ls_reduction_factor")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("richard_ls_reduction_factor",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"min_ls_factor")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("richard_min_ls_factor",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"ls_acceptance_factor")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("richard_ls_acceptance_factor",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"monitor_line_search")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("richard_monitor_line_search",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"monitor_linear_solve")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("richard_monitor_linear_solve",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"use_fd_jac")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    if (strcmp(textContent,"true")==0) {
                      expertPL.set<bool>("richard_use_fd_jac",true);
                    }
                    else  {
                      expertPL.set<bool>("richard_use_fd_jac",false);
                    }
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"perturbation_scale_for_J")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("richard_perturbation_scale_for_J",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"use_dense_Jacobian")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    if (strcmp(textContent,"true")==0) {
                      expertPL.set<bool>("richard_use_dense_Jacobian",true);
                    }
                    else  {
                      expertPL.set<bool>("richard_use_dense_Jacobian",false);
                    }
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"upwind_krel")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    if (strcmp(textContent,"true")==0) {
                      expertPL.set<bool>("richard_upwind_krel",true);
                    }
                    else  {
                      expertPL.set<bool>("richard_upwind_krel",false);
                    }
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"pressure_maxorder")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<int>("richard_pressure_maxorder",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"scale_solution_before_solve")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    if (strcmp(textContent,"true")==0) {
                      expertPL.set<bool>("richard_scale_solution_before_solve",true);
                    }
                    else  {
                      expertPL.set<bool>("richard_scale_solution_before_solve",false);
                    }
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"semi_analytic_J")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    if (strcmp(textContent,"true")==0) {
                      expertPL.set<bool>("richard_semi_analytic_J",true);
                    }
                    else  {
                      expertPL.set<bool>("richard_semi_analytic_J",false);
                    }
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"cfl")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    expertPL.set<double>("cfl",get_double_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
		  else if (strcmp(tagName,"max_n_subcycle_transport")==0) {
		    textContent = XMLString::transcode(currentNode->getTextContent());
		    expertPL.set<int>("max_n_subcycle_transport",get_int_constant(textContent,*def_list));
		    XMLString::release(&textContent);
		  }
                }
              }
            }
            else if (strcmp(nodeName,"str_amr_controls")==0) {
              // loop through children and deal with them
	      Teuchos::Array<std::string> refinement_indicator_names;
              DOMNodeList* children = tmpNode->getChildNodes();
              for (int k=0; k<children->getLength(); k++) {
                DOMNode* currentNode = children->item(k) ;
                if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                  char* tagname = XMLString::transcode(currentNode->getNodeName());
                  if (strcmp(tagname,"amr_levels")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    amrPL.set<int>("Number Of AMR Levels",get_int_constant(textContent,*def_list));
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"refinement_ratio")==0) {
                    Teuchos::Array<int> factors;
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    factors = make_int_list(textContent);
                    amrPL.set<Teuchos::Array<int> >("Refinement Ratio",factors);
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"do_amr_subcycling")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    if (strcmp(textContent,"true")==0) {
                      expertPL.set<bool>("do_amr_subcycling",true);
                    }
                    else  {
                      expertPL.set<bool>("do_amr_subcycling",false);
                    }
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"regrid_interval")==0) {
                    Teuchos::Array<int> factors;
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    factors = make_int_list(textContent);
                    amrPL.set<Teuchos::Array<int> >("Regrid Interval",factors);
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"blocking_factor")==0) {
                    Teuchos::Array<int> factors;
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    factors = make_int_list(textContent);
                    amrPL.set<Teuchos::Array<int> >("Blocking Factor",factors);
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"number_error_buffer_cells")==0) {
                    Teuchos::Array<int> factors;
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    factors = make_int_list(textContent);
                    amrPL.set<Teuchos::Array<int> >("Number Error Buffer Cells",factors);
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"max_grid_size")==0) {
                    Teuchos::Array<int> factors;
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    factors = make_int_list(textContent);
                    amrPL.set<Teuchos::Array<int> >("Maximum Grid Size",factors);
                    XMLString::release(&textContent);
                  }
                  else if (strcmp(tagname,"refinement_indicator")==0) {
                    char* nameRefinement;
                    Teuchos::ParameterList refinePL;
                    bool foundOne = false;
                    bool foundName = false;
                    // loop over attributes to get name
                    attrMap = currentNode->getAttributes();
                    for (int l=0; l<attrMap->getLength(); l++) {
                      nodeAttr = attrMap->item(l);
                      attrName =XMLString::transcode(nodeAttr->getNodeName());
                      textContent = XMLString::transcode(nodeAttr->getNodeValue());
                      if (strcmp(attrName,"name")==0) {
                        nameRefinement = XMLString::transcode(nodeAttr->getNodeValue());
                        foundName = true;
                      }
                    }
                    if (!foundName) {
                      //error message: refinement_indicators must have name
                    }
		    refinement_indicator_names.push_back(std::string(nameRefinement));

                    // loop over children to get values
                    DOMNodeList* kids = currentNode->getChildNodes();
                    for (int l=0; l<kids->getLength(); l++) {
                      DOMNode* currentKid = kids->item(l) ;
                      if (DOMNode::ELEMENT_NODE == currentKid->getNodeType()) {
                        char* tagname = XMLString::transcode(currentKid->getNodeName());
                        if (strcmp(tagname,"field_name")==0) {
                          textContent = XMLString::transcode(currentKid->getTextContent());
                          refinePL.set<std::string>("Field Name",trim_string(textContent));
                          XMLString::release(&textContent);
                        }
                        else if (strcmp(tagname,"regions")==0) {
                          textContent = XMLString::transcode(currentKid->getTextContent());
                          Teuchos::Array<std::string> regs = make_regions_list(textContent);
                          refinePL.set<Teuchos::Array<std::string> >("Regions",regs);
                          XMLString::release(&textContent);
                        }
                        else if (strcmp(tagname,"max_refinement_level")==0) {
                          textContent = XMLString::transcode(currentKid->getTextContent());
                          refinePL.set<int>("Maximum Refinement Level",atoi(textContent));
                          XMLString::release(&textContent);
                        }
                        else if (strcmp(tagname,"start_time")==0) {
                          textContent = XMLString::transcode(currentKid->getTextContent());
                          refinePL.set<double>("Start Time",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else if (strcmp(tagname,"end_time")==0) {
                          textContent = XMLString::transcode(currentKid->getTextContent());
                          refinePL.set<double>("End Time",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else if (strcmp(tagname,"value_greater")==0) {
                          textContent = XMLString::transcode(currentKid->getTextContent());
                          if (!foundOne) {
                            refinePL.set<double>("Value Greater",get_double_constant(textContent,*def_list));
                            XMLString::release(&textContent);
                            foundOne = true;
                          }
                          else {
                            //TODO: EIB - error message, already have one
                          }
                          XMLString::release(&textContent);
                        }
                        else if (strcmp(tagname,"value_less")==0) {
                          textContent = XMLString::transcode(currentKid->getTextContent());
                          if (!foundOne) {
                            refinePL.set<double>("Value Less",get_double_constant(textContent,*def_list));
                            XMLString::release(&textContent);
                            foundOne = true;
                          }
                          else {
                            //TODO: EIB - error message, already have one
                          }
                          XMLString::release(&textContent);
                        }
                        else if (strcmp(tagname,"adjacent_difference_greater")==0) {
                          textContent = XMLString::transcode(currentKid->getTextContent());
                          if (!foundOne) {
                            refinePL.set<double>("Adjacent Difference Greater",get_double_constant(textContent,*def_list));
                            XMLString::release(&textContent);
                            foundOne = true;
                          }
                          else {
                            //TODO: EIB - error message, already have one
                          }
                          XMLString::release(&textContent);
                        }
                        else if (strcmp(tagname,"inside_region")==0) {
                          textContent = XMLString::transcode(currentKid->getTextContent());
                          if (!foundOne) {
                            if (strcmp(textContent,"on")==0) {
                              refinePL.set<bool>("Inside Region",true);
                            }
                            else  {
                              refinePL.set<bool>("Inside Region",false);
                            }
                            //XMLString::release(&textContent);
                            foundOne = true;
                          }
                          else {
                            //TODO: EIB - error message, already have one
                          }
                          XMLString::release(&textContent);
                        }
			XMLString::release(&tagname);
                      }
                    }
                    amrPL.sublist(nameRefinement) = refinePL;
                    XMLString::release(&nameRefinement);
                  }
                }
              }
	      if (refinement_indicator_names.size() > 0) {
		amrPL.set<Teuchos::Array<std::string> >("Refinement Indicators",refinement_indicator_names);
	      }
            }
            else if (strcmp(nodeName,"max_n_subcycle_transport")==0) {
              textContent = XMLString::transcode(tmpNode->getTextContent());
              expertPL.set<int>("max_n_subcycle_transport",get_int_constant(textContent,*def_list));
              XMLString::release(&textContent);
	    }
            else if (strcmp(nodeName,"str_chemistry_controls")==0) {
              Teuchos::ParameterList chemistryPL;
              // chemistry options are difference based on engine
              // deal with common items
              if (voI_->getVerbLevel() == Teuchos::VERB_EXTREME) {
                  chemistryPL.set<std::string>("Verbosity","error");
              }
              else if (voI_->getVerbLevel() == Teuchos::VERB_HIGH) {
                  chemistryPL.set<std::string>("Verbosity","warning");
              }
              else if (voI_->getVerbLevel() == Teuchos::VERB_MEDIUM) {
                  chemistryPL.set<std::string>("Verbosity","verbose");
              }
              else if (voI_->getVerbLevel() == Teuchos::VERB_LOW ) {
                  chemistryPL.set<std::string>("Verbosity","terse");
              }
              else {
                  chemistryPL.set<std::string>("Verbosity","silent");
              }
              
              // determine engine and proceed accordingly
              if (def_list->isParameter("chemistry_engine")) {
                if (def_list->get<std::string>("chemistry_engine") == "amanzi") {
                  // go ahead and add bdg file to PL
                  // build bgd filename
                  std::string bgdfilename;
                  if (def_list->isParameter("xmlfilename") ) {
                    bgdfilename = def_list->get<std::string>("xmlfilename") ;
                    std::string new_extension(".bgd");
                    size_t pos = bgdfilename.find(".xml");
                    bgdfilename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)4);
                  }
                  else {
                    // defaulting to hardcoded name
                    bgdfilename = "isotherms.bgd" ;
                  }
                  // add bgd file and parameters to list
                  Teuchos::ParameterList bgdPL;
                  bgdPL.set<std::string>("Format","simple");
                  bgdPL.set<std::string>("File",bgdfilename);
                  chemistryPL.sublist("Thermodynamic Database") = bgdPL;
                  chemistryPL.set<std::string>("Activity Model","unit");
                  
                  // loop over chemistry controls to get other options to add to PL
                  DOMNodeList* children = tmpNode->getChildNodes();
                  for (int k=0; k<children->getLength(); k++) {
                    DOMNode* currentNode = children->item(k) ;
                    if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                      std::string tagname(std::string(XMLString::transcode(currentNode->getNodeName())));
                      if (tagname == "chem_tolerance") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("Tolerance",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "chem_tolerance", nodeName);
                        }
                      }
                      else if (tagname == "chem_max_newton_iterations") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<int>("Maximum Newton Iterations",get_int_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "chem_max_newton_iterations", nodeName);
                        }
                        // TODO:: EIB - this need to be added to schema!!
                      }
                      else if (tagname == "chem_max_time_step") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("Max Time Step (s)",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "chem_max_time_step", nodeName);
                        }
                      }
                      else if (tagname == "max_chemistry_transport_timestep_ratio") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("max chemistry to transport timestep ratio",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "max_chemistry_transport_timestep_ratio", nodeName);
                        }
                      }
                      else if (tagname == "process_model") {
                        // TODO: EIB - not sure where this goes anymore
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("Tolerance",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        // TODO: EIB - removed error message until I figure out where this went
                        /*
                         else {
                         throw_error_illformed(algo_str_name, "process_model", nodeName);
                         }
                         */
                      }
                      else if (tagname != "comments") {
                        // warn about unrecognized element
                        std::stringstream elem_name;
                        elem_name << nodeName <<  "->"  << tagname;
                        throw_warning_skip(elem_name.str());
                      }
                    }
                  }
                }
                else if (def_list->get<std::string>("chemistry_engine") == "pflotran") {
                  chemistryPL.set<std::string>("Engine","PFloTran");
                  
                  // check file *.in filename in def_list
                  if (def_list->isSublist("chemistry_PL")) {
                    if (def_list->sublist("chemistry_PL").isParameter("Engine Input File")) {
                      chemistryPL.set<std::string>("Engine Input File",def_list->sublist("chemistry_PL").get<std::string>("Engine Input File"));
                    }
                  }
                  
                  // loop over children
                  DOMNodeList* children = tmpNode->getChildNodes();
                  for (int k=0; k<children->getLength(); k++) {
                    DOMNode* currentNode = children->item(k) ;
                    if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
                      std::string tagname(std::string(XMLString::transcode(currentNode->getNodeName())));
                      if (tagname == "chem_max_time_step") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<double>("Max Time Step (s)",get_double_constant(textContent,*def_list));
                          XMLString::release(&textContent);
                        }
                        else {
                          throw_error_illformed(algo_str_name, "chem_max_time_step", nodeName);
                        }
                      }
                      else if (tagname == "generate_chemistry_engine_inputfile") {
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<std::string>("Engine Input File",textContent);
                          XMLString::release(&textContent);
                        }
                      }
                      else if (tagname == "read_chemistry_engine_inputfile") {
                        // TODO: EIB - not sure where this goes anymore
                        if (currentNode) {
                          textContent = XMLString::transcode(currentNode->getTextContent());
                          chemistryPL.set<std::string>("Engine Input File",textContent);
                          XMLString::release(&textContent);
                        }
                      }
                      else if (tagname != "comments") {
                        // warn about unrecognized element
                        std::stringstream elem_name;
                        elem_name << nodeName <<  "->"  << tagname;
                        throw_warning_skip(elem_name.str());
                      }
                    }
                  }
                }
              }
            
              
              // now add chemistry list
              def_list->sublist("Chemistry") = chemistryPL;
            }
            else if (strcmp(nodeName,"petsc_options_file")==0) {
              // This specifies the name of the file containing petsc options
              // The default filename ".petsc" is found automatically
              // Utilize this option if the file uses an alternative name
              textContent = XMLString::transcode(tmpNode->getTextContent());
              def_list->set<std::string>("petsc_options_file",trim_string(textContent));
              XMLString::release(&textContent);
            }
	    XMLString::release(&nodeName);
          }
        }
      }
    }
    list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Expert Settings") = expertPL;
    list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Adaptive Mesh Refinement Control") = amrPL;
  }
  
  if (!isUnstr_ && !found) {
    Errors::Message msg;
    msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing numerical_controls - " ;
    msg << "No element 'structured_controls' found. \n" ;
    msg << "The element is required for structured framework. \n";
    msg << "  Please correct and try again \n" ;
    Exceptions::amanzi_throw(msg);
  }
  if (!isUnstr_ && found_other) {
    //TODO: EIB - ignoring message
  }
  
  return list;
}
  
/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_phases(DOMDocument* xmlDoc, Teuchos::ParameterList& def_list) {

  Teuchos::ParameterList list;

  DOMNodeList* nodeList;
  DOMNodeList* nodeList2;
  DOMNode* nodeTmp;
  DOMNode* nodeTmp2;
  DOMNode* nodeTmp3;
  DOMNode* nodeAttr;
  DOMNode* nodeAttr2;
  DOMNode* nodeAttr3;
  DOMNamedNodeMap* attrMap;
  DOMNamedNodeMap* attrMap2;
  std::string tagName;
  char* textContent;
  char* textContent2;

  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Phases" << std::endl;
  }

  // get phases node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("phases"));
  DOMNode* nodeEC = nodeList->item(0);
  DOMElement* elementEC = static_cast<DOMElement*>(nodeEC);

  // get comments node - EIB: removed set, old spec can't handle comments at this level
  nodeList = elementEC->getElementsByTagName(XMLString::transcode("comments"));
  if (nodeList->getLength() > 0) {
    nodeTmp = nodeList->item(0);
    textContent = XMLString::transcode(nodeTmp->getTextContent());
    //list.set<std::string>("comments",textContent);
    XMLString::release(&textContent);
  }

  // get liquid_phase node
  nodeList = elementEC->getElementsByTagName(XMLString::transcode("liquid_phase"));
  if (nodeList->getLength() > 0) {
    nodeTmp = nodeList->item(0);
    attrMap = nodeTmp->getAttributes();
    nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
    std::string phaseName;
    if (nodeAttr) {
      textContent = XMLString::transcode(nodeAttr->getNodeValue());
      phaseName = std::string(textContent);
      if (phaseName=="water") {
	phaseName="Water";
      }
      XMLString::release(&textContent);
    }
    else {
      throw_error_missattr("phases","attribute","name","liquid_phase");
    }

    DOMNodeList* childern = nodeTmp->getChildNodes();
    bool foundPC = false; // checking for phase components
    for (int i=0; i<childern->getLength(); i++) {
      DOMNode* cur = childern->item(i) ;
      if (DOMNode::ELEMENT_NODE == cur->getNodeType()) {
        DOMElement* propElem = static_cast<DOMElement*>(cur);
        std::string tagName(XMLString::transcode(cur->getNodeName()));
        textContent = XMLString::transcode(cur->getTextContent());
	//TODO: NOTE: EIB - skipping EOS, not currently supported
	if (tagName == "viscosity"){
          if (propElem->hasAttribute(XMLString::transcode("type"))) {
            Teuchos::ParameterList propertyPL;
            propertyPL = get_file_info(propertyPL, static_cast<DOMElement*>(cur), "viscosity", "liquid_phase");
            list.sublist("Aqueous").sublist("Phase Properties").sublist("Viscosity: File") = propertyPL;
            std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
            std::cout << "    Please note - the XML Schema allows for specifing Viscosity with a file read." << std::endl;
            std::cout << "    However, file read for this property has not been propagated through the input process."  << std::endl;
            std::cout << "    Please specify a value."  << std::endl;
          }
          else {
            list.sublist("Aqueous").sublist("Phase Properties").sublist("Viscosity: Uniform").set<double>("Viscosity",get_double_constant(textContent,def_list));
          }
	}
        else if (tagName == "density") {
          if (propElem->hasAttribute(XMLString::transcode("type"))) {
            Teuchos::ParameterList propertyPL;
            propertyPL = get_file_info(propertyPL, static_cast<DOMElement*>(cur), "density", "liquid_phase");
            list.sublist("Aqueous").sublist("Phase Properties").sublist("Density: File") = propertyPL;
            std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
            std::cout << "    Please note - the XML Schema allows for specifing Density with a file read." << std::endl;
            std::cout << "    However, file read for this property has not been propagated through the input process."  << std::endl;
            std::cout << "    Please specify a value."  << std::endl;
          }
          else {
            list.sublist("Aqueous").sublist("Phase Properties").sublist("Density: Uniform").set<double>("Density",get_double_constant(textContent,def_list));
          }
	}
	else if (tagName == "dissolved_components") {
	  Teuchos::ParameterList dcPL;
	  Teuchos::Array<double> diffusion;  // TODO: EIB - not using any diffusion constants right now, where do the go???
	  Teuchos::Array<std::string> solutes;
          
          // EIB: new stuff
          DOMNodeList* gkids = cur->getChildNodes();
          for (int j=0; j<gkids->getLength(); j++) {
            DOMNode* curGKid = gkids->item(j) ;
            if (DOMNode::ELEMENT_NODE == curGKid->getNodeType()) {
              std::string grpName  = XMLString::transcode(curGKid->getNodeName());
              /* TODO: EIB - finish implementing primaries and secondaries
              if (grpName == "primaries") {
                DOMNodeList* ggkids = curGKid->getChildNodes();
                for (int k=0; k<ggkids->getLength(); k++) {
                  DOMNode* curGGKid = ggkids->item(k) ;
                  if (DOMNode::ELEMENT_NODE == curGGKid->getNodeType()) {
                    std::string elemName  = XMLString::transcode(curGGKid->getNodeName());
                    if (elemName == "primary") {
                      // TODO: EIB - get primary name from element text
                      // TODO: EIB - loop over attributes
                      // TODO: EIB - if primary doubles as solute, then need diffusion; if not, skip it here
                    }
                    else {
                      // TODO: EIB - warn that skipping this because i don't recognize it
                    }
                  }
                }
              }
              else if (grpName == "secondaries") {
                DOMNodeList* ggkids = curGKid->getChildNodes();
                for (int k=0; k<ggkids->getLength(); k++) {
                  DOMNode* curGGKid = ggkids->item(k) ;
                  if (DOMNode::ELEMENT_NODE == curGGKid->getNodeType()) {
                    std::string elemName  = XMLString::transcode(curGGKid->getNodeName());
                    if (elemName == "secondary") {
                      // TODO: EIB - do something with this info
                    }
                    else {
                      // TODO: EIB - warn that skipping this because i don't recognize it
                    }
                  }
                }
              }
               */
              if (grpName == "solutes") {
                Teuchos::ParameterList solutesPL;
                DOMNodeList* ggkids = curGKid->getChildNodes();
                for (int k=0; k<ggkids->getLength(); k++) {
                  DOMNode* curGGKid = ggkids->item(k) ;
                  if (DOMNode::ELEMENT_NODE == curGGKid->getNodeType()) {
                    std::string elemName  = XMLString::transcode(curGGKid->getNodeName());
                    if (elemName == "solute") {
                      // TODO: EIB - do something with this info
                      Teuchos::ParameterList solPL;
                      std::string soluteName(XMLString::transcode(curGGKid->getTextContent()));
                      solutes.append(soluteName);
                      attrMap = curGGKid->getAttributes();
                      // put attribute - coefficient_of_diffusion in diffusion array
                      nodeAttr = attrMap->getNamedItem(XMLString::transcode("coefficient_of_diffusion"));
                      if (nodeAttr) {
                        textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
                        if (isUnstr_) {
                          solPL.set<double>("Molecular Diffusivity: Uniform",get_double_constant(textContent2,def_list));
                        }
                        else {
                          solPL.set<double>("Molecular Diffusivity",get_double_constant(textContent2,def_list));
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        if (isUnstr_) {
                          throw_error_missattr("dissolved_components", "attribute", "coefficient_of_diffusion", "solute");
                        }
                      }
                      // put attribute - first_order_decay_constant in diffusion array
                      nodeAttr = attrMap->getNamedItem(XMLString::transcode("first_order_decay_constant"));
                      if (nodeAttr) {
                        textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
                        solPL.set<double>("First Order Decay Constant",get_double_constant(textContent2,def_list));
                        XMLString::release(&textContent2);
                      }
                      // no error message, because this is optional
                      
                      // add parameters to solute PL
                      dcPL.sublist(soluteName) = solPL;
                    }
                    else {
                      // TODO: EIB - warn that skipping this because i don't recognize it
                    }
                  }
                }
                dcPL.set<Teuchos::Array<std::string> >("Component Solutes",solutes);
                list.sublist("Aqueous").sublist("Phase Components").sublist(phaseName ) = dcPL;
		def_list.set<Teuchos::Array<std::string> >("solutes",solutes);
                foundPC = true;
              }
            }
          }
          
          
	}
        XMLString::release(&textContent);
      }
    }
    if (!isUnstr_ && !foundPC) {
      Teuchos::ParameterList emptyPL;
      
      list.sublist("Aqueous").sublist("Phase Components").sublist("Water") = emptyPL;
    }
  }

  // get any solid_phase node
  nodeList = elementEC->getElementsByTagName(XMLString::transcode("solid_phase"));
  if (nodeList->getLength() > 0) {
    Teuchos::ParameterList spPL;
    Teuchos::Array<std::string> minerals;
    nodeTmp = nodeList->item(0);
    DOMElement* solidElem = static_cast<DOMElement*>(nodeTmp);
    nodeList2 = solidElem->getElementsByTagName(XMLString::transcode("minerals"));
    if (nodeList2->getLength() > 0) {
      nodeTmp2 = nodeList2->item(0);
      DOMNodeList* kids = nodeTmp2->getChildNodes();
      for (int i=0; i<kids->getLength(); i++) {
        DOMNode* curKid = kids->item(i) ;
        if (DOMNode::ELEMENT_NODE == curKid->getNodeType()) {
          tagName  = XMLString::transcode(curKid->getNodeName());
	  if (tagName == "mineral"){
            textContent2 = XMLString::transcode(curKid->getTextContent());
	    minerals.append(textContent2);
            XMLString::release(&textContent2);
	  }
	}
      }
    }
    spPL.set<Teuchos::Array<std::string> >("Minerals",minerals);
    list.sublist("Solid") = spPL;
  }

  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_regions(DOMDocument* xmlDoc, Teuchos::ParameterList* def_list) {

  Teuchos::ParameterList list;

  DOMNodeList* nodeList;
  DOMNode* nodeTmp;
  DOMNode* nodeAttr;
  DOMNamedNodeMap* attrMap;
  char* tagName;
  char* nodeName;
  char* regName;
  char* textContent;
  char* textContent2;
  char* char_array;

  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Regions" << std::endl;
  }
  Teuchos::ParameterList reg_names;

  // if structured, add default region names know inside
  if (!isUnstr_) {
    regionNames_string_.append("XLOBC");
    regionNames_string_.append("XHIBC");
    regionNames_string_.append("YLOBC");
    regionNames_string_.append("YHIBC");
    regionNames_string_.append("ZLOBC");
    regionNames_string_.append("ZHIBC");
    regionNames_string_.append("All");
  }
  
  // get regions node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("regions"));
  DOMNode* nodeRgn = nodeList->item(0);
  DOMElement* elementRgn = static_cast<DOMElement*>(nodeRgn);

  // just loop over the children and deal with them as they come
  // new options: comment, region, box, point
  DOMNodeList* childern = nodeRgn->getChildNodes();
  for (int i=0; i<childern->getLength(); i++) {
    DOMNode* cur = childern->item(i) ;
    if (DOMNode::ELEMENT_NODE == cur->getNodeType()) {
      tagName  = XMLString::transcode(cur->getNodeName());
      /* NOTE: EIB - geometry doesn't deal with extra comment node
      if (strcmp(tagName,"comments") == 0){
        textContent = XMLString::transcode(cur->getTextContent());
        list.set<std::string>("comments",textContent);
        XMLString::release(&textContent);
      } */
      bool haveName = false;
      
      /////////
      DOMElement* regElem;
      // if region is under a region tag,
      // get the name and set child element as region element
      if  (strcmp(tagName,"region") == 0){
        // get region name
        attrMap = cur->getAttributes();
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
        if (nodeAttr) {
          regName = XMLString::transcode(nodeAttr->getNodeValue());
          haveName = true;
        } else {
          throw_error_missattr("Regions","attribute","name","region");
        }
        // loop over children to get region element
        DOMNodeList* kids = cur->getChildNodes();
        for (int j=0; j<kids->getLength(); j++) {
          DOMNode* curKid = kids->item(j) ;
          if (DOMNode::ELEMENT_NODE == curKid->getNodeType()) {
            nodeName  = XMLString::transcode(curKid->getNodeName());
            if (strcmp(nodeName,"comments") != 0)
              regElem =static_cast<DOMElement*>(curKid);
            XMLString::release(&nodeName);
          }
        }
      }
      // else set the current element as region element
      else {
        regElem =static_cast<DOMElement*>(cur);
      }
      
      // get regElem type
      nodeName  = XMLString::transcode(regElem->getNodeName());
      
      // get name if needed
      if (!haveName) {
        if (regElem->hasAttribute(XMLString::transcode("name"))) {
          regName = XMLString::transcode(regElem->getAttribute(XMLString::transcode("name")));
        }
        else {
          if (strcmp(nodeName,"comments") != 0)
            throw_error_missattr("Regions","attribute","name",nodeName);
        }
      }
      
      // loop over attributes of region element
      
      if  (strcmp(nodeName,"box") == 0){
        
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"box");
          regionNames_string_.append(regName);
        }
        
        if (regElem->hasAttribute(XMLString::transcode("low_coordinates"))) {
          textContent2 = XMLString::transcode(regElem->getAttribute(XMLString::transcode("low_coordinates")));
        } else {
          throw_error_missattr("Regions","attribute","low_coordinates","box");
        }
        Teuchos::Array<double> low = make_coordinates(textContent2, *def_list);
        list.sublist(regName).sublist("Region: Box").set<Teuchos::Array<double> >("Low Coordinate",low);
        XMLString::release(&textContent2);
        
        if (regElem->hasAttribute(XMLString::transcode("high_coordinates"))) {
          textContent2 = XMLString::transcode(regElem->getAttribute(XMLString::transcode("high_coordinates")));
        } else {
          throw_error_missattr("Regions","attribute","high_coordinates","box");
        }
        Teuchos::Array<double> high = make_coordinates(textContent2,* def_list);
        list.sublist(regName).sublist("Region: Box").set<Teuchos::Array<double> >("High Coordinate",high);
        XMLString::release(&textContent2);
      }
      else if  (strcmp(nodeName,"plane") == 0){
        
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"plane");
          regionNames_string_.append(regName);
        }
        
        if (regElem->hasAttribute(XMLString::transcode("location"))) {
          textContent2 = XMLString::transcode(regElem->getAttribute(XMLString::transcode("location")));
        } else {
          throw_error_missattr("Regions","attribute","location","plane");
        }
        Teuchos::Array<double> loc = make_coordinates(textContent2, *def_list);
        list.sublist(regName).sublist("Region: Plane").set<Teuchos::Array<double> >("Location",loc);
        
        if (regElem->hasAttribute(XMLString::transcode("normal"))) {
          textContent2 = XMLString::transcode(regElem->getAttribute(XMLString::transcode("normal")));
        } else {
          throw_error_missattr("Regions","attribute","normal","plane");
        }
        Teuchos::Array<double> dir = make_coordinates(textContent2, *def_list);
        list.sublist(regName).sublist("Region: Plane").set<Teuchos::Array<double> >("Direction",dir);
        XMLString::release(&textContent2);
      }
      else if  (strcmp(nodeName,"region_file") == 0){
        //TODO: EIB - add file
        Teuchos::ParameterList rfPL;
        
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"region_file");
          regionNames_string_.append(regName);
        }
        
        if (regElem->hasAttribute(XMLString::transcode("name"))) {
          textContent2 = XMLString::transcode(regElem->getAttribute(XMLString::transcode("name")));
        } else {
          throw_error_missattr("Regions","attribute","name","region_file");
        }
        rfPL.set<std::string>("File",trim_string(textContent2));
        XMLString::release(&textContent2);
        
        if (regElem->hasAttribute(XMLString::transcode("type"))) {
          textContent2 = XMLString::transcode(regElem->getAttribute(XMLString::transcode("type")));
        } else {
          throw_error_missattr("Regions","attribute","type","region_file");
        }
        if  (strcmp(textContent2,"color") == 0){
          char* value;
          if (regElem->hasAttribute(XMLString::transcode("label"))) {
            value = XMLString::transcode(regElem->getAttribute(XMLString::transcode("label")));
          } else {
            throw_error_missattr("Regions","attribute","label","color");
          }
          rfPL.set<int>("Value",atoi(value));
          XMLString::release(&value);
          list.sublist(regName).sublist("Region: Color Function") = rfPL;
        }
        else if  (strcmp(textContent2,"labeled set") == 0){
          char* value ;
          if (regElem->hasAttribute(XMLString::transcode("label"))) {
            value = XMLString::transcode(regElem->getAttribute(XMLString::transcode("label")));
          } else {
            throw_error_missattr("Regions","attribute","label","labeled set");
          }
          rfPL.set<std::string>("Label",trim_string(value));
          XMLString::release(&value);
          
          if (regElem->hasAttribute(XMLString::transcode("format"))) {
            value = XMLString::transcode(regElem->getAttribute(XMLString::transcode("format")));
          } else {
            throw_error_missattr("Regions","attribute","format","labeled set");
          }
          if  (strcmp(value,"exodus ii") == 0){
            rfPL.set<std::string>("Format","Exodus II");
          }
          XMLString::release(&value);
          
          if (regElem->hasAttribute(XMLString::transcode("entity"))) {
            value = XMLString::transcode(regElem->getAttribute(XMLString::transcode("entity")));
          } else {
            throw_error_missattr("Regions","attribute","entity","labeled set");
          }
          rfPL.set<std::string>("Entity",trim_string(value));
          XMLString::release(&value);
          
          list.sublist(regName).sublist("Region: Labeled Set") = rfPL;
        }
        XMLString::release(&textContent2);
      }
      else if  (strcmp(nodeName,"point") == 0){
        
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"point");
          regionNames_string_.append(regName);
        }
        
        if (regElem->hasAttribute(XMLString::transcode("coordinate"))) {
          textContent2 = XMLString::transcode(regElem->getAttribute(XMLString::transcode("coordinate")));
        } else {
          throw_error_missattr("Regions","attribute","coordinate","point");
        }
        
        Teuchos::Array<double> coord = make_coordinates(textContent2, *def_list);
        list.sublist(regName).sublist("Region: Point").set<Teuchos::Array<double> >("Coordinate",coord);
        XMLString::release(&textContent2);
      }
      else if  (strcmp(nodeName,"polygonal_surface") == 0){
        if (!isUnstr_) {
          throw_error_str_ustr("Regions", tagName, "Unstructured");
        }
        
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"polygonal_surface");
          regionNames_string_.append(regName);
        }
        
         // if attribute 'num_points' exists, get it
        int num_points(-1);
        int pt_cnt(0);
        if (regElem->hasAttribute(XMLString::transcode("num_points"))) {
          textContent2 = XMLString::transcode(regElem->getAttribute(XMLString::transcode("num_points")));
          std::string str(textContent2);
          boost::algorithm::trim(str);
          num_points = atoi(textContent2);
          //list.sublist(regName).sublist("Region: Polygonal Surface").set<int>("Number of points",num_points);
          list.sublist(regName).sublist("Region: Polygon").set<int>("Number of points",num_points);
          XMLString::release(&textContent2);
        }
        if (regElem->hasAttribute(XMLString::transcode("tolerance"))) {
          textContent2 = XMLString::transcode(regElem->getAttribute(XMLString::transcode("tolerance")));
          list.sublist(regName).sublist("Region: Polygon").sublist("Expert Parameters").set<double>("Tolerance",get_double_constant(textContent2,*def_list));
          XMLString::release(&textContent2);
        }
        // get verticies (add count them)
        Teuchos::Array<double> points;
        DOMNodeList* gkids = regElem->getChildNodes();
        for (int j=0; j<gkids->getLength(); j++) {
          DOMNode* curGKid = gkids->item(j) ;
          if (DOMNode::ELEMENT_NODE == curGKid->getNodeType()) {
            nodeName  = XMLString::transcode(curGKid->getNodeName());
            if  (strcmp(nodeName,"point") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              Teuchos::Array<double> point = make_coordinates(textContent2, *def_list);
              for (Teuchos::Array<double>::iterator pt=point.begin(); pt!=point.end(); ++pt) {
                points.append(*pt);
              }
              pt_cnt++;
              XMLString::release(&textContent2);
            }
          }
        }
        //list.sublist(regName).sublist("Region: Polygonal Surface").set<Teuchos::Array<double> >("Points",points);
        list.sublist(regName).sublist("Region: Polygon").set<Teuchos::Array<double> >("Points",points);
        // check that 'num_points' with current count
        //if (!list.sublist(regName).sublist("Region: Polygonal Surface").isParameter("Number of points")) {
        //  list.sublist(regName).sublist("Region: Polygonal Surface").set<int>("Number of points",pt_cnt);
        //}
        if (!list.sublist(regName).sublist("Region: Polygon").isParameter("Number of points")) {
          list.sublist(regName).sublist("Region: Polygon").set<int>("Number of points",pt_cnt);
        }
      }
      else if  (strcmp(nodeName,"logical") == 0){
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"logical");
          regionNames_string_.append(regName);
        }
       
        // Loop over childern
	bool haveOp = false;
	bool haveRL = false;
	Teuchos::Array<Teuchos::Array<double> > points;
        DOMNodeList* gkids = regElem->getChildNodes();
        for (int j=0; j<gkids->getLength(); j++) {
          DOMNode* curGKid = gkids->item(j) ;
          if (DOMNode::ELEMENT_NODE == curGKid->getNodeType()) {
            nodeName  = XMLString::transcode(curGKid->getNodeName());
	    // deal with operation
            if  (strcmp(nodeName,"operation") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              if ( strcmp(textContent2,"union") == 0) {
                list.sublist(regName).sublist("Region: Logical").set<std::string>("Operation","Union");
              }
              else if (strcmp(textContent2,"intersection") == 0) {
                list.sublist(regName).sublist("Region: Logical").set<std::string>("Operation","Intersection");
              }
              else if (strcmp(textContent2,"subtraction") == 0) {
                list.sublist(regName).sublist("Region: Logical").set<std::string>("Operation","Subtraction");
              }
              else if (strcmp(textContent2,"complement") == 0) {
                list.sublist(regName).sublist("Region: Logical").set<std::string>("Operation","Complement");
              }
	      else {
		// error about missing or unsupported operation
		throw_error_illformed("Regions","element","operation","union, intersection, subtraction, or complement");
	      }
	      haveOp = true;
              XMLString::release(&textContent2);
            }
	    // deal with region list
	    else if  (strcmp(nodeName,"region_list") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
	      Teuchos::Array<std::string> regs = make_regions_list(textContent2);
              list.sublist(regName).sublist("Region: Logical").set<Teuchos::Array<std::string> >("Regions",regs);
	      haveRL = true;
              XMLString::release(&textContent2);
	    }
          }
        }
	if (!haveOp) {
	  throw_error_missattr("Regions","element","operation","logical");
	}
	if (!haveRL) {
	  throw_error_missattr("Regions","element","region_list","logical");
	}

      }
      else if  (strcmp(nodeName,"polygon") == 0){
        if (isUnstr_) {
          throw_error_str_ustr("Regions", tagName, "Structured");
        }
        if (dimension_ != 2) {
          Errors::Message msg;
          msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing Regions - " ;
          msg << "Regions type 'polygon' is only available for 2D. \n" ;
          msg << "  Please correct and try again \n" ;
          Exceptions::amanzi_throw(msg);
        }
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"polygon");
          regionNames_string_.append(regName);
        }
        
        // num_points for error checking
        char* numPoints;
        if (regElem->hasAttribute(XMLString::transcode("num_points"))) {
          numPoints = XMLString::transcode(regElem->getAttribute(XMLString::transcode("num_points")));
        } else {
          throw_error_missattr("Regions","attribute","num_points","polygon");
        }
        
        // get verticies (add count them)
        Teuchos::Array<double> pointsX;
        Teuchos::Array<double> pointsY;
        DOMNodeList* gkids = regElem->getChildNodes();
        for (int j=0; j<gkids->getLength(); j++) {
          DOMNode* curGKid = gkids->item(j) ;
          if (DOMNode::ELEMENT_NODE == curGKid->getNodeType()) {
            nodeName  = XMLString::transcode(curGKid->getNodeName());
            if  (strcmp(nodeName,"point") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              Teuchos::Array<double> point = make_coordinates(textContent2, *def_list);
              pointsX.append(point[0]);
              pointsY.append(point[1]);
              XMLString::release(&textContent2);
            }
          }
        }
        list.sublist(regName).sublist("Region: Polygon").set<Teuchos::Array<double> >("VerticesV1",pointsX);
        list.sublist(regName).sublist("Region: Polygon").set<Teuchos::Array<double> >("VerticesV2",pointsY);
      }
      else if  (strcmp(nodeName,"ellipse") == 0){
        if (isUnstr_) {
          throw_error_str_ustr("Regions", tagName, "Structured");
        }
        if (dimension_ != 2) {
          Errors::Message msg;
          msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing Regions - " ;
          msg << "Regions type 'ellipse' is only available for 2D. \n" ;
          msg << "  Please correct and try again \n" ;
          Exceptions::amanzi_throw(msg);
        }
        
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"ellipse");
          regionNames_string_.append(regName);
        }
        
        // get verticies (add count them)
        Teuchos::Array<Teuchos::Array<double> > points;
        DOMNodeList* gkids = regElem->getChildNodes();
        for (int j=0; j<gkids->getLength(); j++) {
          DOMNode* curGKid = gkids->item(j) ;
          if (DOMNode::ELEMENT_NODE == curGKid->getNodeType()) {
            nodeName  = XMLString::transcode(curGKid->getNodeName());
            if (strcmp(nodeName,"center") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              Teuchos::Array<double> center = make_coordinates(textContent2, *def_list);
              list.sublist(regName).sublist("Region: Ellipse").set<Teuchos::Array<double> >("Center",center);
              XMLString::release(&textContent2);
            }
            else if (strcmp(nodeName,"radius") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              Teuchos::Array<double> radius = make_coordinates(textContent2, *def_list);
              list.sublist(regName).sublist("Region: Ellipse").set<Teuchos::Array<double> >("Radius",radius);
              XMLString::release(&textContent2);
            }
          }
        }
      }
      else if  (strcmp(nodeName,"rotated_polygon") == 0){
        if (isUnstr_) {
          throw_error_str_ustr("Regions", tagName, "Structured");
        }
        if (dimension_ != 3) {
          Errors::Message msg;
          msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing Regions - " ;
          msg << "Regions type 'rotated_polygon' is only available for 3D. \n" ;
          msg << "  Please correct and try again \n" ;
          Exceptions::amanzi_throw(msg);
        }
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"rotated_polygon");
          regionNames_string_.append(regName);
        }
        
        // get verticies (add count them)
        Teuchos::Array<double> pointsX;
        Teuchos::Array<double> pointsY;
        Teuchos::Array<double> pointsZ;
        DOMNodeList* gkids = regElem->getChildNodes();
        for (int j=0; j<gkids->getLength(); j++) {
          DOMNode* curGKid = gkids->item(j) ;
          if (DOMNode::ELEMENT_NODE == curGKid->getNodeType()) {
            nodeName  = XMLString::transcode(curGKid->getNodeName());
            if  (strcmp(nodeName,"point") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              Teuchos::Array<double> point = make_coordinates(textContent2, *def_list);
              pointsX.append(point[0]);
              pointsY.append(point[1]);
              pointsZ.append(point[2]);
              XMLString::release(&textContent2);
            }
            else if  (strcmp(nodeName,"reference_point") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              Teuchos::Array<double> point = make_coordinates(textContent2, *def_list);
              list.sublist(regName).sublist("Region: Rotated Polygon").set<Teuchos::Array<double> >("Reference Point", point);
              XMLString::release(&textContent2);
            }
            else if  (strcmp(nodeName,"plane") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              if ( strcmp(textContent2,"xy") == 0 | strcmp(textContent2,"yx") == 0) {
                list.sublist(regName).sublist("Region: Rotated Polygon").set<std::string>("Plane","XY");
              }
              else if (strcmp(textContent2,"yz") == 0 | strcmp(textContent2,"zy") == 0) {
                list.sublist(regName).sublist("Region: Rotated Polygon").set<std::string>("Plane","YZ");
              }
              else if (strcmp(textContent2,"xz") == 0 | strcmp(textContent2,"zx") == 0) {
                list.sublist(regName).sublist("Region: Rotated Polygon").set<std::string>("Plane","XZ");
              }
              else {
                throw_error_illformed("Regions", "value", "plane", "xy, yx, xz");
              }
              XMLString::release(&textContent2);
            }
            else if  (strcmp(nodeName,"axis") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              if ( strcmp(textContent2,"x") == 0) {
                list.sublist(regName).sublist("Region: Rotated Polygon").set<std::string>("Plane","X");
              }
              else if (strcmp(textContent2,"y") == 0) {
                list.sublist(regName).sublist("Region: Rotated Polygon").set<std::string>("Plane","Y");
              }
              else if (strcmp(textContent2,"z") == 0) {
                list.sublist(regName).sublist("Region: Rotated Polygon").set<std::string>("Plane","Z");
              }
              else {
                throw_error_illformed("Regions", "value", "axis", "x, y, z");
              }
              XMLString::release(&textContent2);
            }
          }
        }
        list.sublist(regName).sublist("Region: Rotated Polygon").set<Teuchos::Array<double> >("VerticesV1",pointsX);
        list.sublist(regName).sublist("Region: Rotated Polygon").set<Teuchos::Array<double> >("VerticesV2",pointsY);
        list.sublist(regName).sublist("Region: Rotated Polygon").set<Teuchos::Array<double> >("VerticesV3",pointsZ);
      }
      else if  (strcmp(nodeName,"swept_polygon") == 0){
        if (isUnstr_) {
          throw_error_str_ustr("Regions", tagName, "Structured");
        }
        if (dimension_ != 3) {
          Errors::Message msg;
          msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing Regions - " ;
          msg << "Regions type 'swept_polygon' is only available for 3D. \n" ;
          msg << "  Please correct and try again \n" ;
          Exceptions::amanzi_throw(msg);
        }
        // add region name to array of region names
        if (reg_names.isParameter(regName)) {
          // warn, region of this name already exists, overwriting
        } else {
          reg_names.set<std::string>(regName,"swept_polygon");
          regionNames_string_.append(regName);
        }
        // get verticies (add count them)
        Teuchos::Array<double> pointsX;
        Teuchos::Array<double> pointsY;
        Teuchos::Array<double> pointsZ;
        Teuchos::Array<double> extent;
        extent.append(0.0);
        extent.append(0.0);
        DOMNodeList* gkids = regElem->getChildNodes();
        for (int j=0; j<gkids->getLength(); j++) {
          DOMNode* curGKid = gkids->item(j) ;
          if (DOMNode::ELEMENT_NODE == curGKid->getNodeType()) {
            nodeName  = XMLString::transcode(curGKid->getNodeName());
            if  (strcmp(nodeName,"point") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              Teuchos::Array<double> point = make_coordinates(textContent2, *def_list);
              pointsX.append(point[0]);
              pointsY.append(point[1]);
              pointsZ.append(point[2]);
              XMLString::release(&textContent2);
            }
            else if  (strcmp(nodeName,"extent_min") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              extent[0] = atof(textContent2);
              XMLString::release(&textContent2);
            }
            else if  (strcmp(nodeName,"extent_max") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              extent[0] = atof(textContent2);
              XMLString::release(&textContent2);
            }
            else if  (strcmp(nodeName,"plane") == 0){
              textContent2 = XMLString::transcode(curGKid->getTextContent());
              if ( strcmp(textContent2,"xy") == 0 | strcmp(textContent2,"yx") == 0) {
                list.sublist(regName).sublist("Region: Swept Polygon").set<std::string>("Plane","XY");
              }
              else if (strcmp(textContent2,"yz") == 0 | strcmp(textContent2,"zy") == 0) {
                list.sublist(regName).sublist("Region: Swept Polygon").set<std::string>("Plane","YZ");
              }
              else if (strcmp(textContent2,"xz") == 0 | strcmp(textContent2,"zx") == 0) {
                list.sublist(regName).sublist("Region: Swept Polygon").set<std::string>("Plane","XZ");
              }
              else {
                throw_error_illformed("Regions", "value", "plane", "xy, yx, xz");
              }
              XMLString::release(&textContent2);
            }
          }
        }
        list.sublist(regName).sublist("Region: Swept Polygon").set<Teuchos::Array<double> >("VerticesV1",pointsX);
        list.sublist(regName).sublist("Region: Swept Polygon").set<Teuchos::Array<double> >("VerticesV2",pointsY);
        list.sublist(regName).sublist("Region: Swept Polygon").set<Teuchos::Array<double> >("VerticesV3",pointsZ);
        list.sublist(regName).sublist("Region: Swept Polygon").set<Teuchos::Array<double> >("Extent",extent);
      }
      
      XMLString::release(&nodeName);
      
      XMLString::release(&tagName);
      if (haveName) XMLString::release(&regName);
    }
  }
  // add array of region names to def_list, use these names to check assigned_regions list against later
  def_list->sublist("regions") = reg_names; 
  
  return list;
  
}


/*
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_materials(DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  Teuchos::ParameterList list;

  DOMNodeList* nodeList;
  DOMNode* nodeTmp;
  DOMNode* nodeTmp2;
  DOMNode* nodeAttr;
  DOMNamedNodeMap* attrMap;
  char* tagName;
  char* propName;
  char* propValue;
  char* textContent;
  char* textContent2;
  char* char_array;
  char* attrName;
  char* attrValue;
  std::string matName;
  bool hasPerm = false;
  bool hasHC = false;

  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Materials" << std::endl;
  }

  // get regions node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("materials"));
  DOMNode* nodeMat = nodeList->item(0);
  DOMElement* elementMat = static_cast<DOMElement*>(nodeMat);

  // just loop over the children and deal with them as they come
  DOMNodeList* childern = nodeMat->getChildNodes();
  for (int i=0; i<childern->getLength(); i++) {
    bool cappressON = false;
    bool dispersionON = false;
    bool diffusionON = false;
    bool tortuosityON = false;
    Teuchos::ParameterList caplist;
    std::string capname;
    hasPerm = false;
    hasHC = false;
    DOMNode* cur = childern->item(i) ;
    if (DOMNode::ELEMENT_NODE == cur->getNodeType()) {
      attrMap = cur->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
      if (nodeAttr) {
        textContent = XMLString::transcode(nodeAttr->getNodeValue());
        matName = std::string(textContent);
        XMLString::release(&textContent);
      }
      else {
        throw_error_missattr("materials", "attribute", "name", "material");
      }

      Teuchos::ParameterList matlist(matName);
      DOMNodeList* kids = cur->getChildNodes();
      for (int j=0; j<kids->getLength(); j++) {
        DOMNode* curkid = kids->item(j) ;
        if (DOMNode::ELEMENT_NODE == curkid->getNodeType()) {
          tagName  = XMLString::transcode(curkid->getNodeName());
	  if (strcmp("assigned_regions",tagName)==0){
            //TODO: EIB - if this is more than 1 region -> assuming comma seperated list of strings????
            textContent2 = XMLString::transcode(curkid->getTextContent());
            Teuchos::Array<std::string> regs = make_regions_list(textContent2);
            matlist.set<Teuchos::Array<std::string> >("Assigned Regions",regs);
            XMLString::release(&textContent2);
            if (!compare_region_names(regs, def_list)) {
              Errors::Message msg;
              msg << "Amanzi::InputTranslator: ERROR - invalid region in Materials Section - " ;
              msg << "valid regions are: \n" ;
              for (int r=0; r<regionNames_string_.size(); r++) {
                msg << "    " << regionNames_string_[r] << "\n";
              }
              msg << "  Please correct and try again \n" ;
              Exceptions::amanzi_throw(msg);
            }
	  }
          else if  (strcmp("mechanical_properties",tagName)==0){
            DOMNodeList* list = curkid->getChildNodes();
            // loop over child: deal with porosity and density
            for (int k=0; k<list->getLength(); k++) {
              DOMNode* curkiddy = list->item(k) ;
              if (DOMNode::ELEMENT_NODE == curkiddy->getNodeType()) {
                propName  = XMLString::transcode(curkiddy->getNodeName());
                if  (strcmp("porosity",propName)==0){
                  DOMElement* propElem = static_cast<DOMElement*>(curkiddy);
                  Teuchos::ParameterList propertyPL;

		  const std::string type_str_DEF("nonfile");
		  const std::string type_str_gslib("gslib");
                  const std::string type_str_file("file");
                  const std::string type_str_exo("exodus ii");

		  std::string type_str(type_str_DEF);

		  if (propElem->hasAttribute(XMLString::transcode("type"))) {
		    type_str = std::string(XMLString::transcode(propElem->getAttribute(XMLString::transcode("type"))));
		  }

		  if (type_str != type_str_DEF) {

		    std::string porName;
		    if (type_str == type_str_gslib) {

		      const std::string dfile_DEF("porosity_data");
		      get_gslib_info(propertyPL, def_list, propElem, "porosity", "mechanical_properties", dfile_DEF);
		      porName = "Porosity: GSLib";

		    } else if (type_str == type_str_file || type_str == type_str_exo) {

		      propertyPL = get_file_info(propertyPL, propElem, "porosity", "mechanical_properties");
		      porName = "Porosity: File";

		      matlist.sublist("Porosity: Uniform") = propertyPL;
		      std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
		      std::cout << "    Please note - the XML Schema allows for specifing Porosity with a file read." << std::endl;
		      std::cout << "    However, file read for this property has not been propagated through the input process."  << std::endl;
		      std::cout << "    Please specify a value."  << std::endl;

		    } else {
                      
                      std::string section;
                      section = "Material (" + matName + ") - porosity";
                      throw_error_illformed(section, "attribute", "type", "gslib or file");

		    }
		    
		    matlist.sublist(porName) = propertyPL;

		  }
		  else {

                    if (propElem->hasAttribute(XMLString::transcode("value"))) {
                      textContent2 = XMLString::transcode(propElem->getAttribute(XMLString::transcode("value")));
                      matlist.sublist("Porosity: Uniform").set<double>("Value",get_double_constant(textContent2,def_list));
                      XMLString::release(&textContent2);
                    }else {
                      throw_error_missattr("material", "attribute", "value", "porosity");
                    }
                  }
                }
                else if  (strcmp("particle_density",propName)==0){
                  // TODO: EIB - should be check value >= 0.
                  DOMElement* propElem = static_cast<DOMElement*>(curkiddy);
                  Teuchos::ParameterList propertyPL;
                  if (propElem->hasAttribute(XMLString::transcode("type"))) {
                    propertyPL = get_file_info(propertyPL, propElem, "particle_density", "mechanical_properties");
                    
                    matlist.sublist("Particle Density: Uniform") = propertyPL;
                    std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
                    std::cout << "    Please note - the XML Schema allows for specifing Particle Density with a file read." << std::endl;
                    std::cout << "    However, file read for this property has not been propagated through the input process."  << std::endl;
                    std::cout << "    Please specify a value."  << std::endl;
                  }
                  else {
                    if (propElem->hasAttribute(XMLString::transcode("value"))) {
                      textContent2 = XMLString::transcode(propElem->getAttribute(XMLString::transcode("value")));
                      matlist.sublist("Particle Density: Uniform").set<double>("Value",get_double_constant(textContent2,def_list));
                      XMLString::release(&textContent2);
                    }else {
                      throw_error_missattr("material", "attribute", "value", "particle_density");
                    }
                  }
                }
                else if  (strcmp("specific_storage",propName)==0){
                  DOMElement* propElem = static_cast<DOMElement*>(curkiddy);
                  Teuchos::ParameterList propertyPL;
                  if (propElem->hasAttribute(XMLString::transcode("type"))) {
                    propertyPL = get_file_info(propertyPL, propElem, "specific_storage", "mechanical_properties");
                    
                    matlist.sublist("Specific Storage: Uniform") = propertyPL;
                    std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
                    std::cout << "    Please note - the XML Schema allows for specifing Specific Storage with a file read." << std::endl;
                    std::cout << "    However, file read for this property has not been propagated through the input process."  << std::endl;
                    std::cout << "    Please specify a value."  << std::endl;
                  }
                  else {
                    if (propElem->hasAttribute(XMLString::transcode("value"))) {
                      textContent2 = XMLString::transcode(propElem->getAttribute(XMLString::transcode("value")));
                      matlist.sublist("Specific Storage: Uniform").set<double>("Value",get_double_constant(textContent2,def_list));
                      XMLString::release(&textContent2);
                    }else {
                      throw_error_missattr("material", "attribute", "value", "specific_storage");
                    }
                  }
                }
                else if  (strcmp("specific_yield",propName)==0){
                  DOMElement* propElem = static_cast<DOMElement*>(curkiddy);
                  Teuchos::ParameterList propertyPL;
                  if (propElem->hasAttribute(XMLString::transcode("type"))) {
                    propertyPL = get_file_info(propertyPL, propElem, "specific_yield", "mechanical_properties");
                    matlist.sublist("Specific Yield: Uniform") = propertyPL;
                    std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
                    std::cout << "    Please note - the XML Schema allows for specifing Specific Yield with a file read." << std::endl;
                    std::cout << "    However, file read for this property has not been propagated through the input process."  << std::endl;
                    std::cout << "    Please specify a value."  << std::endl;
                  }
                  else {
                    if (propElem->hasAttribute(XMLString::transcode("value"))) {
                      textContent2 = XMLString::transcode(propElem->getAttribute(XMLString::transcode("value")));
                      matlist.sublist("Specific Yield: Uniform").set<double>("Value",get_double_constant(textContent2,def_list));
                      XMLString::release(&textContent2);
                    }else {
                      throw_error_missattr("material", "attribute", "value", "specific_yield");
                    }
                  }

                }
                else if  (strcmp("dispersion_tensor",propName)==0){
                  // TODO: EIB - not handling file case
                  DOMElement* dispElem = static_cast<DOMElement*>(curkiddy);
                  if (dispElem->hasAttribute(XMLString::transcode("type"))) {
                    textContent2 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("type")));
                    char* textContent3;
                    if (std::string(textContent2) == "uniform_isotropic") {
                      if (dispElem->hasAttribute(XMLString::transcode("alpha_l"))) {
                        textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_l")));
                        matlist.sublist("Dispersion Tensor: Uniform Isotropic").set<double>("alphaL",get_double_constant(textContent3,def_list));
                        XMLString::release(&textContent3);
                      }else {
                        throw_error_missattr("material", "attribute", "alpha_l", "dispersion_tensor");
                      }
                      if (dispElem->hasAttribute(XMLString::transcode("alpha_t"))) {
                        textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_t")));
                        matlist.sublist("Dispersion Tensor: Uniform Isotropic").set<double>("alphaT",get_double_constant(textContent3,def_list));
                        XMLString::release(&textContent3);
                      }else {
                        throw_error_missattr("material", "attribute", "alpha_t", "dispersion_tensor");
                      }
                    }
                    else if (std::string(textContent2) == "burnett_frind") {
                      if (dispElem->hasAttribute(XMLString::transcode("alpha_l"))) {
                        textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_l")));
                        matlist.sublist("Dispersion Tensor: Burnett-Frind").set<double>("alphaL",get_double_constant(textContent3,def_list));
                        XMLString::release(&textContent3);
                      }else {
                        throw_error_missattr("material", "attribute", "alpha_l", "dispersion_tensor");
                      }
                      if (dispElem->hasAttribute(XMLString::transcode("alpha_th"))) {
                        textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_th")));
                        matlist.sublist("Dispersion Tensor: Burnett-Frind").set<double>("alphaTH",get_double_constant(textContent3,def_list));
                        XMLString::release(&textContent3);
                      }else {
                        throw_error_missattr("material", "attribute", "alpha_th", "dispersion_tensor");
                      }
                      if (dispElem->hasAttribute(XMLString::transcode("alpha_tv"))) {
                        textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_tv")));
                        matlist.sublist("Dispersion Tensor: Burnett-Frind").set<double>("alphaTV",get_double_constant(textContent3,def_list));
                        XMLString::release(&textContent3);
                      }else {
                        throw_error_missattr("material", "attribute", "alpha_tv", "dispersion_tensor");
                      }
                    }
                    else if (std::string(textContent2) == "lichtner_kelkar_robinson") {
                      if (dispElem->hasAttribute(XMLString::transcode("alpha_lh"))) {
                        textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_lh")));
                        matlist.sublist("Dispersion Tensor: Lichtner-Kelkar-Robinson").set<double>("alphaLH",get_double_constant(textContent3,def_list));
                        XMLString::release(&textContent3);
                      }
                      else {
                        throw_error_missattr("material", "attribute", "alpha_lh", "dispersion_tensor");
                      }
                      if (dispElem->hasAttribute(XMLString::transcode("alpha_lv"))) {
                        textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_lv")));
                        matlist.sublist("Dispersion Tensor: Lichtner-Kelkar-Robinson").set<double>("alphaLV",get_double_constant(textContent3,def_list));
                        XMLString::release(&textContent3);
                      }
                      else {
                        throw_error_missattr("material", "attribute", "alpha_lv", "dispersion_tensor");
                      }
                      if (dispElem->hasAttribute(XMLString::transcode("alpha_th"))) {
                        textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_th")));
                        matlist.sublist("Dispersion Tensor: Lichtner-Kelkar-Robinson").set<double>("alphaTH",get_double_constant(textContent3,def_list));
                        XMLString::release(&textContent3);
                      }
                      else {
                        throw_error_missattr("material", "attribute", "alpha_th", "dispersion_tensor");
                      }
                      if (dispElem->hasAttribute(XMLString::transcode("alpha_tv"))) {
                        textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_tv")));
                        matlist.sublist("Dispersion Tensor: Lichtner-Kelkar-Robinson").set<double>("alphaTV",get_double_constant(textContent3,def_list));
                        XMLString::release(&textContent3);
                      }else {
                        throw_error_missattr("material", "attribute", "alpha_tv", "dispersion_tensor");
                      }
                    }
                    XMLString::release(&textContent2);
                    dispersionON = true;
                  }
                  else {
                    char* textContent3;
                    if (dispElem->hasAttribute(XMLString::transcode("alpha_l"))) {
                      textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_l")));
                      matlist.sublist("Dispersion Tensor: Uniform Isotropic").set<double>("alphaL",get_double_constant(textContent3,def_list));
                      XMLString::release(&textContent3);
                    }
                    else {
                      throw_error_missattr("material", "attribute", "alpha_l", "dispersion_tensor");
                    }
                    if (dispElem->hasAttribute(XMLString::transcode("alpha_t"))) {
                      textContent3 = XMLString::transcode(dispElem->getAttribute(XMLString::transcode("alpha_t")));
                      matlist.sublist("Dispersion Tensor: Uniform Isotropic").set<double>("alphaT",get_double_constant(textContent3,def_list));
                      XMLString::release(&textContent3);
                    }
                    else {
                      throw_error_missattr("material", "attribute", "alpha_t", "dispersion_tensor");
                    }
                    dispersionON = true;
                    std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
                    std::cout << "    Please note - the XML Schema for specifing dispersion_tensor has been updated." << std::endl;
                    std::cout << "    See the Amanzi User Guide for the latest information." << std::endl;
                    std::cout <<"     The legacy format is being handled for now.  Please update input files for future versions." << std::endl;
                    //TODO: EIB - once synced with Akuna, remove above and just have error below.
                    //throw_error_missattr("material", "attribute", "type", "dispersion_tensor");
                  }
                }
                else if  (strcmp("tortuosity",propName)==0){
                  DOMElement* propElem = static_cast<DOMElement*>(curkiddy);
                  Teuchos::ParameterList propertyPL;
                  if (propElem->hasAttribute(XMLString::transcode("type"))) {
                    propertyPL = get_file_info(propertyPL, propElem, "tortuosity", "mechanical_properties");
                    matlist.sublist("Tortuosity Aqueous: Uniform") = propertyPL;
                    std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
                    std::cout << "    Please note - the XML Schema allows for specifing Tortuosity with a file read." << std::endl;
                    std::cout << "    However, file read for this property has not been propagated through the input process."  << std::endl;
                    std::cout << "    Please specify a value."  << std::endl;
                  }
                  else {
                    if (propElem->hasAttribute(XMLString::transcode("value"))) {
                      textContent2 = XMLString::transcode(propElem->getAttribute(XMLString::transcode("value")));
                      matlist.sublist("Tortuosity Aqueous: Uniform").set<double>("Value",get_double_constant(textContent2,def_list));
                      XMLString::release(&textContent2);
                    }else {
                      throw_error_missattr("material", "attribute", "value", "tortuosity");
                    }
                  }
                }
              }
            }
	  }
          else if  (strcmp("permeability",tagName)==0) {
            DOMElement* permElem = static_cast<DOMElement*>(curkid);
            Teuchos::ParameterList perm;
            Teuchos::ParameterList permTmp;
            std::string permName("Intrinsic Permeability: Anisotropic Uniform");
            if (hasHC || hasPerm) {
              Errors::Message msg;
              msg << "Amanzi::InputTranslator: ERROR - can only specify one of ";
              msg << "'permeability' or 'hydraulic_conductivity' per material. Material - " << matName << "\n";
              msg << "  Please correct and try again \n" ;
              Exceptions::amanzi_throw(msg);
            }
            hasPerm = true;

	    const std::string type_str_DEF("nonfile");
	    const std::string type_str_gslib("gslib");
            const std::string type_str_file("file");
            const std::string type_str_exo("exodus ii");

	    std::string type_str(type_str_DEF);

	    if (permElem->hasAttribute(XMLString::transcode("type"))) {
	      type_str = std::string(XMLString::transcode(permElem->getAttribute(XMLString::transcode("type"))));
	    }

	    if (type_str != type_str_DEF) {

	      if (type_str == type_str_gslib) {

		const std::string dfile_DEF("permeability_data");
		get_gslib_info(perm, def_list, permElem, "permeability", "materials", dfile_DEF);
		permName = "Intrinsic Permeability: GSLib";

	      } else if (type_str == type_str_file || type_str == type_str_exo) {

		perm = get_file_info(perm, permElem, "permeability", "materials");
		permName = "Intrinsic Permeability: File";

	      } else {
                
                std::string section;
                section = "Material (" + matName + ") - permeability";
                throw_error_illformed(section, "attribute", "type", "gslib or file");

	      }
	    }
            else {
              // loop over attributes to get x,y,z
              char *x,*y,*z;
              attrMap = curkid->getAttributes();
              for (int k=0; k<attrMap->getLength(); k++) {
                DOMNode* attrNode = attrMap->item(k) ;
                if (DOMNode::ATTRIBUTE_NODE == attrNode->getNodeType()) {
                  char* attrName = XMLString::transcode(attrNode->getNodeName());
                  if (strcmp("x",attrName)==0){
                    x = XMLString::transcode(attrNode->getNodeValue());
                    permTmp.set<double>("x",get_double_constant(x,def_list));
                    XMLString::release(&x);
                  } else if (strcmp("y",attrName)==0){
                    y = XMLString::transcode(attrNode->getNodeValue());
                    permTmp.set<double>("y",get_double_constant(y,def_list));
                    XMLString::release(&y);
                  } else if (strcmp("z",attrName)==0){
                    z = XMLString::transcode(attrNode->getNodeValue());
                    permTmp.set<double>("z",get_double_constant(z,def_list));
                    XMLString::release(&z);
                  }
                }
              }
              bool checkY = true , checkZ = true;
              if (permTmp.isParameter("y")) {
                if (permTmp.get<double>("x") != permTmp.get<double>("y")) checkY = false;
              }
              if (permTmp.isParameter("z")) {
                if (permTmp.get<double>("x") != permTmp.get<double>("z")) checkZ = false;
              }
              if (checkY && checkZ) {
                perm.set<double>("Value",permTmp.get<double>("x"));
                //matlist.sublist("Intrinsic Permeability: Uniform") = perm;
                permName = "Intrinsic Permeability: Uniform";
              }
              else {
                perm.set<double>("x",permTmp.get<double>("x"));
                if (permTmp.isParameter("y")) {
                  perm.set<double>("y",permTmp.get<double>("y")) ;
                }
                if (permTmp.isParameter("z")) {
                  perm.set<double>("z",permTmp.get<double>("z")) ;
                }
              }
            }
            matlist.sublist(permName) = perm;
	  }
          else if  (strcmp("hydraulic_conductivity",tagName)==0) {
            DOMElement* hydcondElem = static_cast<DOMElement*>(curkid);
            Teuchos::ParameterList hydcond;
            Teuchos::ParameterList hydcondTmp;
            std::string hydcondName("Hydraulic Conductivity: Anisotropic Uniform");
            if (hasHC || hasPerm) {
              Errors::Message msg;
              msg << "Amanzi::InputTranslator: ERROR - can only specify one of ";
              msg << "'permeability' or 'hydraulic_conductivity' per material. Material - " << matName << "\n";
              msg << "  Please correct and try again \n" ;
              Exceptions::amanzi_throw(msg);
            }
            hasHC = true;

	    const std::string type_str_DEF("nonfile");
	    const std::string type_str_gslib("gslib");
	    const std::string type_str_file("file");

	    std::string type_str(type_str_DEF);

	    if (hydcondElem->hasAttribute(XMLString::transcode("type"))) {
	      type_str = std::string(XMLString::transcode(hydcondElem->getAttribute(XMLString::transcode("type"))));
	    }

	    if (type_str != type_str_DEF) {

	      if (type_str == type_str_gslib) {

		const std::string dfile_DEF("hydraulic_conductivity_data");
		get_gslib_info(hydcond, def_list, hydcondElem, "hydraulic_conductivity", "materials", dfile_DEF);
		hydcondName = "Hydraulic Conductivity: GSLib";

	      } else if (type_str == type_str_file) {

		hydcond = get_file_info(hydcond, hydcondElem, "hydraulic_conductivity", "materials");
		hydcondName = "Hydraulic Conductivity: File";

	      } else {

		Errors::Message msg;
		msg << "Amanzi::InputTranslator: ERROR - hydraulic_conductivity 'type' can only be either";
		msg << " 'gslib' or 'file' per material. Material - " << matName << "\n";
		msg << "  Please correct and try again \n" ;
		Exceptions::amanzi_throw(msg);

	      }
	    }
            else {
              // loop over attributes to get x,y,z
              char *x,*y,*z;
              attrMap = curkid->getAttributes();
              for (int k=0; k<attrMap->getLength(); k++) {
                DOMNode* attrNode = attrMap->item(k) ;
                if (DOMNode::ATTRIBUTE_NODE == attrNode->getNodeType()) {
                  char* attrName = XMLString::transcode(attrNode->getNodeName());
                  if (strcmp("x",attrName)==0){
                    x = XMLString::transcode(attrNode->getNodeValue());
                    hydcondTmp.set<double>("x",get_double_constant(x,def_list));
                    XMLString::release(&x);
                  } else if (strcmp("y",attrName)==0){
                    y = XMLString::transcode(attrNode->getNodeValue());
                    hydcondTmp.set<double>("y",get_double_constant(y,def_list));
                    XMLString::release(&y);
                  } else if (strcmp("z",attrName)==0){
                    z = XMLString::transcode(attrNode->getNodeValue());
                    hydcondTmp.set<double>("z",get_double_constant(z,def_list));
                    XMLString::release(&z);
                  }
                }
              }
              bool checkY = true , checkZ = true;
              if (hydcondTmp.isParameter("y")) {
                if (hydcondTmp.get<double>("x") != hydcondTmp.get<double>("y")) checkY = false;
              }
              if (hydcondTmp.isParameter("z")) {
                if (hydcondTmp.get<double>("x") != hydcondTmp.get<double>("z")) checkZ = false;
              }
              if (checkY && checkZ) {
                hydcond.set<double>("Value",hydcondTmp.get<double>("x"));
                hydcondName = "Hydraulic Conductivity: Uniform";
              }
              else {
                hydcond.set<double>("x",hydcondTmp.get<double>("x"));
                if (hydcondTmp.isParameter("y")) {
                  hydcond.set<double>("y",hydcondTmp.get<double>("y")) ;
                }
                if (hydcondTmp.isParameter("z")) {
                  hydcond.set<double>("z",hydcondTmp.get<double>("z")) ;
                }
              }
            }
            matlist.sublist(hydcondName) = hydcond;
	  }
          else if  (strcmp("cap_pressure",tagName)==0){
            attrMap = curkid->getAttributes();
            nodeAttr = attrMap->getNamedItem(XMLString::transcode("model"));
            if (nodeAttr) {
              textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
            }
            else {
              throw_error_missattr("material", "attribute", "model", "cap_pressure");
            }

            if  (strcmp("van_genuchten",textContent2)==0){
              // TODO: EIB - not handling file case
              cappressON = true;
              DOMNodeList* paramList= curkid->getChildNodes();
              for (int k=0; k<paramList->getLength(); k++) {
                DOMNode* paramNode = paramList->item(k) ;
                if (DOMNode::ELEMENT_NODE == paramNode->getNodeType()) {
                  propName  = XMLString::transcode(paramNode->getNodeName());
                  if  (strcmp("parameters",propName)==0){
                    attrMap = paramNode->getAttributes();
                    for (int l=0; l<attrMap->getLength(); l++) {
                      DOMNode* attr = attrMap->item(l) ;
                      attrName  = XMLString::transcode(attr->getNodeName());
                      attrValue  = XMLString::transcode(attr->getNodeValue());
                      if (strcmp(attrName,"sr")==0) {
                        caplist.set<double>("Sr",get_double_constant(attrValue,def_list));
                      }
                      else if (strcmp(attrName,"optional_krel_smoothing_interval")==0) {
                        caplist.set<double>("krel smoothing interval",get_double_constant(attrValue,def_list));
                      }
                      else {
                        caplist.set<double>(attrName,get_double_constant(attrValue,def_list));
                      }
                      XMLString::release(&attrName);
                      XMLString::release(&attrValue);
                    }
                  }
                  else if (strcmp("read",propName)==0) {
                    // TODO: EIB - implement read, how to seperate parameters?
                    std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
                    std::cout << "    Please note - the XML Schema allows for specifing cap_pressure:van_genuchten with a file read." << std::endl;
                    std::cout << "    However, file read for this property has not been propagated through the input process."  << std::endl;
                    std::cout << "    Please specify parameters."  << std::endl;
                  }
                  else {
                    // ill-formed error
                    throw_error_illformed("material","element","van_genuchten","parameters or read");
                  }
                }
              }
              capname = "Capillary Pressure: van Genuchten";
            }
            else if (strcmp("brooks_corey",textContent2)==0){
              // TODO: EIB - not handling file case
              cappressON = true;
              DOMNodeList* paramList= curkid->getChildNodes();
              for (int k=0; k<paramList->getLength(); k++) {
                DOMNode* paramNode = paramList->item(k) ;
                if (DOMNode::ELEMENT_NODE == paramNode->getNodeType()) {
                  propName  = XMLString::transcode(paramNode->getNodeName());
                  if  (strcmp("parameters",propName)==0){
                    attrMap = paramNode->getAttributes();
                    for (int l=0; l<attrMap->getLength(); l++) {
                      DOMNode* attr = attrMap->item(l) ;
                      attrName  = XMLString::transcode(attr->getNodeName());
                      attrValue  = XMLString::transcode(attr->getNodeValue());
                      if (strcmp(attrName,"sr")==0) {
                        caplist.set<double>("Sr",get_double_constant(attrValue,def_list));
                      }
                      else if (strcmp(attrName,"optional_krel_smoothing_interval")==0) {
                        caplist.set<double>("krel smoothing interval",get_double_constant(attrValue,def_list));
                      }
                      else {
                        caplist.set<double>(attrName,get_double_constant(attrValue,def_list));
                      }
                      XMLString::release(&attrName);
                      XMLString::release(&attrValue);
                    }
                  }
                }
              }
              capname = "Capillary Pressure: Brooks Corey";
            }
	    XMLString::release(&textContent2);
	  }
          else if  (strcmp("rel_perm",tagName)==0){
            // TODO: EIB - how to handle if cappress=false? ie, caplist not setup yet?
            attrMap = curkid->getAttributes();
            nodeAttr = attrMap->getNamedItem(XMLString::transcode("model"));
            if (nodeAttr) {
              textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
            }
            else {
              throw_error_missattr("material", "attribute", "model", "rel_perm");
            }

            if (strcmp(textContent2,"burdine")==0) {
              caplist.set<std::string>("Relative Permeability","Burdine");
              // TODO: EIB - should get child element exp here. Don't know where it goes in v1.2.1 (ELL_BURDINE=2??)
            }
            else if (strcmp(textContent2,"mualem")==0) {
              caplist.set<std::string>("Relative Permeability","Mualem");
            }
            if (strcmp(textContent2,"none")!=0) {
              DOMNodeList* paramList= curkid->getChildNodes();
              for (int k=0; k<paramList->getLength(); k++) {
                DOMNode* paramNode = paramList->item(k) ;
                if (DOMNode::ELEMENT_NODE == paramNode->getNodeType()) {
                  propName  = XMLString::transcode(paramNode->getNodeName());
                  if  (strcmp("optional_krel_smoothing_interval",propName)==0){
                    propValue  = XMLString::transcode(paramNode->getTextContent());
                    caplist.set<double>("krel smoothing interval",get_double_constant(propValue,def_list));
                  }
                }
              }
            }

	  }
          else if  (strcmp("sorption_isotherms",tagName)==0){
            // TODO: barker - write bgd file to go with this!!!!
	    // TODO: barker - should check that chemistry is on and set to amanzi native
            Teuchos::ParameterList sorptionPL;
	    // loop over child: deal with solutes
            DOMNodeList* list = curkid->getChildNodes();
            for (int k=0; k<list->getLength(); k++) {
              DOMNode* curkiddy = list->item(k) ;
              if (DOMNode::ELEMENT_NODE == curkiddy->getNodeType()) {
                propName  = XMLString::transcode(curkiddy->getNodeName());
		// error checking to make sure this is a solute element, all properties are in attributes
	        if  (strcmp("solute",propName)==0){
                  Teuchos::ParameterList solutePL;
		  char *soluteName;
		  char *modelName;
                  // get solute name (attribute)
                  attrMap = curkiddy->getAttributes();
                  nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
                  if (nodeAttr) {
                    soluteName = XMLString::transcode(nodeAttr->getNodeValue());
                  }
                  else {
                    //TODO: EIB - add error msg here
                  }
                  // loop over children
                  DOMNodeList* propElementList = curkiddy->getChildNodes();
                  for (int l=0; l<propElementList->getLength(); l++) {
                    DOMNode* propNode = propElementList->item(l) ;
                    if (DOMNode::ELEMENT_NODE == propNode->getNodeType()) {
                      char* propertyName  = XMLString::transcode(propNode->getNodeName());
                      DOMNamedNodeMap* propAttrMap;
                      // Look for molecular_diffusion element (moved to phases)
                      /*
                      if (strcmp("molecular_diffusion",propertyName)==0){
                        // loop over molecular_diffusion attributes
                        propAttrMap = propNode->getAttributes();
                        for (int m=0; m<propAttrMap->getLength(); m++) {
                          DOMNode* attrNode = propAttrMap->item(m) ;
                          if (DOMNode::ATTRIBUTE_NODE == attrNode->getNodeType()) {
                            char* attrName = XMLString::transcode(attrNode->getNodeName());
                            if (strcmp("type",attrName)==0) {
                              // should get type
                              // then get filename
                              Errors::Message msg;
                              msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing solutes - " ;
                              msg << "This translator does not handle molecular_diffusion file reads yet. \n" ;
                              msg << "  Please specify a value and try again \n" ;
                              Exceptions::amanzi_throw(msg);
                            }
                            if (strcmp("value",attrName)==0) {
                              //solutePL.set<double>("Molecular Diffusivity: Uniform 1",get_double_constant(XMLString::transcode(attrNode->getNodeValue()),def_list));
                              //TODO: EIB - this is a hack for now!!!!
                              diffusionON = true;
                              matlist.sublist("Molecular Diffusivity: Uniform").set<double>("Value",get_double_constant(XMLString::transcode(attrNode->getNodeValue()),def_list));
                            }
                          }
                        }
                      }
                       */
                      // Look for kd_model element (optional)
                      if (strcmp("kd_model",propertyName)==0){
                        // loop over kd_model attributes
                        propAttrMap = propNode->getAttributes();
                        for (int m=0; m<propAttrMap->getLength(); m++) {
                          DOMNode* attrNode = propAttrMap->item(m) ;
                          if (DOMNode::ATTRIBUTE_NODE == attrNode->getNodeType()) {
                            char* attrName = XMLString::transcode(attrNode->getNodeName());
                            if (strcmp("model",attrName)==0){
                              modelName = XMLString::transcode(attrNode->getNodeValue());
                              //TODO: EIB - add checking of modelName = linear | langmuir | freundlich
                            } else if (strcmp("kd",attrName)==0){
                              solutePL.set<double>("Kd",get_double_constant(XMLString::transcode(attrNode->getNodeValue()),def_list));
                            } else if (strcmp("b",attrName)==0){
                              solutePL.set<double>("Langmuir b",get_double_constant(XMLString::transcode(attrNode->getNodeValue()),def_list));
                            } else if (strcmp("n",attrName)==0){
                              solutePL.set<double>("Freundlich n",get_double_constant(XMLString::transcode(attrNode->getNodeValue()),def_list));
                            }
                          }
                        }
                      }
                    }
                  }
                  sorptionPL.sublist(soluteName) = solutePL;
		}
	      }
	    }
	    matlist.sublist("Sorption Isotherms") = sorptionPL;
	    // write BGD file
            if (def_list.isParameter("chemistry_engine")) {
              if (def_list.get<std::string>("chemistry_engine") == "amanzi") {
                write_BDG_file(sorptionPL, def_list);
              }
            }
            else {
              write_BDG_file(sorptionPL, def_list);
            }
	    // Chemistry list is also necessary - this is created under numerical controls section
	  }
 	  XMLString::release(&tagName);
          // If dispersion or diffusion is on, need the other.  This is a hack to correct for user forgetting one.
          // doesn't work right since moveing diffusion to phases
          /*
          if (dispersionON != diffusionON != tortuosityON)
          {
            // If dispersion is on: default Moleculare diffusion and Tortuosity
            if (dispersionON) {
              // TODO: EIB - add notification message here
              if (!diffusionON) {
                matlist.sublist("Molecular Diffusivity: Uniform").set<double>("Value",0.0);
              }
              if (!tortuosityON) {
                matlist.sublist("Tortuosity: Uniform").set<double>("Value",0.0);
              }
            }
            // else if dispersion is off, remove Molecular Diffusivity and Tortuosity
            else {
              // TODO: EIB - add notification message here
              if (matlist.isSublist("Molecular Diffusivity: Uniform")) {
                matlist.remove("Molecular Diffusivity: Uniform");
              }
              if (matlist.isSublist("Tortuosity: Uniform")) {
                matlist.remove("Tortuosity: Uniform");
              }
            }
          }
           */
	}
      }
      if(cappressON) matlist.sublist(capname) = caplist;
      list.sublist(matName) = matlist;
      if (!hasPerm and !hasHC){
        Teuchos::ParameterList perm;
        perm.set<double>("Value",0.0);
        list.sublist(matName).sublist("Intrinsic Permeability: Uniform") = perm;
      }
    }
  }

  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_initial_conditions(DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  Teuchos::ParameterList list;

  DOMNodeList* nodeList;
  DOMNode* nodeTmp;
  DOMNode* nodeAttr;
  DOMNamedNodeMap* attrMap;
  char* tagName;
  char* propName;
  char* textContent;
  char* textContent2;
  char* char_array;
  char* attrName;
  char* attrValue;
  std::string phaseName;

  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Initial Conditions" << std::endl;
  }

  // get regions node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("initial_conditions"));
  DOMNode* nodeIC = nodeList->item(0);
  DOMElement* elementIC = static_cast<DOMElement*>(nodeIC);

  // just loop over the children and deal with them as they come
  DOMNodeList* childern = nodeIC->getChildNodes();
  for (int i=0; i<childern->getLength(); i++) {
    DOMNode* cur = childern->item(i) ;
    if (DOMNode::ELEMENT_NODE == cur->getNodeType()) {
      // get name of IC, then loop over it's children to fill it in
      attrMap = cur->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
      if (nodeAttr) {
        textContent = XMLString::transcode(nodeAttr->getNodeValue());
      }
      else {
        throw_error_missattr("initial_conditions", "attribute", "name", "initial_condition");
      }

      Teuchos::ParameterList iclist(textContent);
      DOMNodeList* IC = cur->getChildNodes();
      for (int j=0; j<IC->getLength(); j++) {
        DOMNode* ICNode = IC->item(j) ;
        tagName  = XMLString::transcode(ICNode->getNodeName());
        //NOTE: EIB - ignoring comments for now
        if (strcmp(tagName,"assigned_regions")==0) {
	  //TODO: EIB - if this is more than 1 region -> assuming comma seperated list of strings????
          textContent2 = XMLString::transcode(ICNode->getTextContent());
	  Teuchos::Array<std::string> regs = make_regions_list(textContent2);
	  iclist.set<Teuchos::Array<std::string> >("Assigned Regions",regs);
	  XMLString::release(&textContent2);
	  if (!compare_region_names(regs, def_list)) {
            Errors::Message msg;
            msg << "Amanzi::InputTranslator: ERROR - invalid region in Initial Conditions Section - " ;
            msg << "valid regions are: \n" ;
            for (int r=0; r<regionNames_string_.size(); r++) {
              msg << "    " << regionNames_string_[r]<< "\n";
            }
            msg << "  Please correct and try again \n" ;
            Exceptions::amanzi_throw(msg);
	  }
        }
        else if (strcmp(tagName,"liquid_phase")==0) {
          //TODO: EIB - deal with liquid phase
          attrMap = ICNode->getAttributes();
          nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
	  if (nodeAttr) {
            textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
            phaseName = std::string(textContent2);
            if (phaseName=="water") {
              phaseName="Water";
            }
            XMLString::release(&textContent2);
	  }
          else {
            throw_error_missattr("initial_conditions", "attribute", "name", "liquid_phase");
	  }

	  // loop over children, deal with liquid_component, solute_component, geomchemistry
          DOMNodeList* compList = ICNode->getChildNodes();
          for (int k=0; k<compList->getLength(); k++) {
            //TODO: EIB - add ELEMENT check
            DOMNode* compNode = compList->item(k) ;
            char* compName  = XMLString::transcode(compNode->getNodeName());
            if (strcmp(compName,"liquid_component")==0) {
	      // loop over children to find pressure
              DOMNodeList* childList = compNode->getChildNodes();
              for (int l=0; l<childList->getLength(); l++) {
                DOMNode* pressure = childList->item(l) ;
                char* pressName  = XMLString::transcode(pressure->getNodeName());
	        Teuchos::ParameterList pressureList;
                if (isUnstr_) pressureList.set<std::string>("Phase","Aqueous");
                if (strcmp(pressName,"uniform_pressure")==0 || strcmp(pressName,"uniform_saturation")==0) {
	          // loop over attributes to get info
	          attrMap = pressure->getAttributes();
                  nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
		  if (nodeAttr) {
                    attrValue = XMLString::transcode(nodeAttr->getNodeValue());
	          } else {
                    throw_error_missattr("initial_conditions", "attribute", "value", "uniform_pressure or uniform_saturation");
	          }

	          pressureList.set<double>("Value",get_double_constant(attrValue,def_list));
	          XMLString::release(&attrValue);
                  if (strcmp(pressName,"uniform_pressure")==0 ) {
                    iclist.sublist("IC: Uniform Pressure") = pressureList;
		  } else {
                    iclist.sublist("IC: Uniform Saturation") = pressureList;
		  }
		}
		else if (strcmp(pressName,"linear_pressure")==0 || strcmp(pressName,"linear_saturation")==0) {
	            char* char_array;
		    //value
	            attrMap = pressure->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
		    if (nodeAttr) {
                      attrValue = XMLString::transcode(nodeAttr->getNodeValue());
                      pressureList.set<double>("Reference Value",get_double_constant(attrValue,def_list));
                      XMLString::release(&attrValue);
	            } else {
                      throw_error_missattr("initial_conditions", "attribute", "value", "linear_pressure or linear_saturation");
	            }

		    //reference_coord
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("reference_coord"));
		    if (nodeAttr) {
                      attrValue = XMLString::transcode(nodeAttr->getNodeValue());
                      Teuchos::Array<double> coord = make_coordinates(attrValue, def_list);
                      pressureList.set<Teuchos::Array<double> >("Reference Point",coord);
                      XMLString::release(&attrValue);
	            } else {
                      throw_error_missattr("initial_conditions", "attribute", "reference_coord", "linear_pressure or linear_saturation");
	            }
		    //gradient
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("gradient"));
		    if (nodeAttr) {
                      attrValue = XMLString::transcode(nodeAttr->getNodeValue());
                      Teuchos::Array<double> grad = make_coordinates(attrValue, def_list);
                      pressureList.set<Teuchos::Array<double> >("Gradient Value",grad);
                      XMLString::release(&attrValue);
	            } else {
                      throw_error_missattr("initial_conditions", "attribute", "gradient", "linear_pressure or linear_saturation");
	            }
                    if (strcmp(pressName,"linear_pressure")==0 ) {
		      iclist.sublist("IC: Linear Pressure") = pressureList;
		    } else {
		      iclist.sublist("IC: Linear Saturation") = pressureList;
		    }
		}
		else if (strcmp(pressName,"velocity")==0 ) {
		    Teuchos::Array<double> vel_vector;
	            attrMap = pressure->getAttributes();
		    // get velocity vector
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("x"));
		    if (nodeAttr) {
                      attrValue = XMLString::transcode(nodeAttr->getNodeValue());
		      vel_vector.append(get_double_constant(attrValue, def_list));
	              XMLString::release(&attrValue);
		    }
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("y"));
		    if (nodeAttr) {
                      attrValue = XMLString::transcode(nodeAttr->getNodeValue());
		      vel_vector.append(get_double_constant(attrValue, def_list));
	              XMLString::release(&attrValue);
		    }
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("z"));
		    if (nodeAttr) {
                      attrValue = XMLString::transcode(nodeAttr->getNodeValue());
		      vel_vector.append(get_double_constant(attrValue, def_list));
	              XMLString::release(&attrValue);
		    }
		    // check that vector dimension == problem dimension
		    if (vel_vector.length() != dimension_) {
		      Errors::Message msg;
		      msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing initial_conditions - " ;
		      msg << "velocity vectory size does not match the spatial dimension of the problem. \n" ;
		      msg << "  Please correct and try again \n" ;
		      Exceptions::amanzi_throw(msg);
		    }
		    // set up new element
		    // TODO:: EIB - does function="uniform | linear" translate to IC: Uniform Velocity and IC: Linear Velocity???
		    //              Linear Velocity does not exist => need error message as such
		    pressureList.set<Teuchos::Array<double> >("Velocity Vector",vel_vector);
		    iclist.sublist("IC: Uniform Velocity") = pressureList;

		}
	      }
			      
	    }
	    else if (strcmp(compName,"solute_component")==0) {
	      char* solName;
	      char* funcType;
	      Teuchos::ParameterList sclist;
	      attrMap = compNode->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
	      if (nodeAttr) {
                solName = XMLString::transcode(nodeAttr->getNodeValue());
	      }
              else {
                throw_error_missattr("initial_conditions", "attribute", "name", "solute_component");
	      }

	      attrMap = compNode->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("function"));
	      if (nodeAttr) {
                funcType = XMLString::transcode(nodeAttr->getNodeValue());
	      }
              else {
                throw_error_missattr("initial_conditions", "attribute", "function", "solute_component");
	      }
	      if (strcmp(funcType,"uniform")==0){
                nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
		if (nodeAttr) {
                  textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
		  sclist.sublist("IC: Uniform Concentration").set<double>("Value",get_double_constant(textContent2,def_list));
		  XMLString::release(&textContent2);
	        }
                else {
                  throw_error_missattr("initial_conditions", "attribute", "value", "solute_component");
	        }
	      }
              else if (strcmp(funcType,"linear")==0){
		// TODO: EIB - currently can't handle this
	      }
	      //TODO: EIB - not added concerntation units, confused by what to add. grab from units?
	      iclist.sublist("Solute IC").sublist("Aqueous").sublist(phaseName).sublist(solName) = sclist;
	      XMLString::release(&solName);
	      XMLString::release(&funcType);
	    }
	    else if (strcmp(compName,"geochemistry")==0) {

	      const Teuchos::Array<std::string>& soluteNames = def_list.get<Teuchos::Array<std::string> >("solutes");

    	      if (soluteNames.size() > 0) {
    	    	char* funcType;
    	    	Teuchos::ParameterList sclist;
    	    	attrMap = compNode->getAttributes();
    	    	nodeAttr = attrMap->getNamedItem(XMLString::transcode("function"));
    	    	if (nodeAttr) {
    	    	  funcType = XMLString::transcode(nodeAttr->getNodeValue());
    	    	}
    	    	else {
    	    	  throw_error_missattr("initial_conditions", "attribute", "function", "geochemistry");
    	    	}
    	    	if (strcmp(funcType,"uniform")==0) {
    	    	  nodeAttr = attrMap->getNamedItem(XMLString::transcode("constraint"));
    	    	  if (nodeAttr) {
    	    	    textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
    	    	    for (int i=0; i<soluteNames.size(); ++i) {
    	    	      Teuchos::ParameterList sic, sic1;
		      sic.set<std::string>("Geochemical Condition",trim_string(textContent2));
    	    	      sic1.sublist("IC: Uniform Concentration") = sic;
    	    	      sclist.sublist(soluteNames[i]) = sic1;
    	    	    }
    	    	    XMLString::release(&textContent2);
    	    	  }
    	    	  else {
    	    	    throw_error_missattr("initial_conditions", "attribute", "value", "solute_component");
    	    	  }
    	    	}
    	    	else if (strcmp(funcType,"linear")==0) {
    	    	  // TODO: EIB - currently can't handle this
    	    	}
    	    	//TODO: EIB - not added concerntation units, confused by what to add. grab from units?
    	    	iclist.sublist("Solute IC").sublist("Aqueous").sublist(phaseName) = sclist;
    	    	XMLString::release(&funcType);
    	      }
	    }
	  }
	}
        else if (strcmp(tagName,"solid_phase")==0) {
          //TODO: EIB - deal with solid phase -> mineral, geochemisty
        }
        XMLString::release(&tagName);
      }
      list.sublist(textContent) = iclist;
      XMLString::release(&textContent);
    }
  }

  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_boundary_conditions(DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  Teuchos::ParameterList list;

  DOMNodeList* nodeList;
  DOMNodeList* BCList;
  DOMNode* nodeTmp;
  DOMNode* nodeAttr;
  DOMNamedNodeMap* attrMap;
  char* tagName;
  char* propName;
  std::string phaseName("Water");
  char* textContent;
  char* textContent2;
  char* char_array;
  char* attrName;
  char* attrValue;


  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Boundary Conditions" << std::endl;
  }

  // get BCs node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("boundary_conditions"));

  if (nodeList->getLength() > 0 ){ // boundary conditions tag does not have to exist

    DOMNode* nodeBC = nodeList->item(0);
    DOMElement* elementBC = static_cast<DOMElement*>(nodeBC);

    // get list of BCs
    BCList = elementBC->getElementsByTagName(XMLString::transcode("boundary_condition"));
    for (int i=0; i<BCList->getLength(); i++) {
      DOMNode* cur = BCList->item(i) ;
      if (DOMNode::ELEMENT_NODE == cur->getNodeType()) {
        // get name of BC, then loop over it's children to fill it in
        attrMap = cur->getAttributes();
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
        if (nodeAttr) {
          textContent = XMLString::transcode(nodeAttr->getNodeValue());
        }
        else {
          throw_error_missattr("boundary_conditions", "attribute", "name", "boundary_condition");
        }

        Teuchos::ParameterList bclist(textContent);
        DOMNodeList* BC = cur->getChildNodes();
        for (int j=0; j<BC->getLength(); j++) {
          DOMNode* BCNode = BC->item(j) ;
          tagName  = XMLString::transcode(BCNode->getNodeName());
          //NOTE: EIB - ignoring comments for now
          if (strcmp(tagName,"assigned_regions")==0) {
        //TODO: EIB - if this is more than 1 region -> assuming comma seperated list of strings????
            textContent2 = XMLString::transcode(BCNode->getTextContent());
	    Teuchos::Array<std::string> regs = make_regions_list(textContent2);
	    bclist.set<Teuchos::Array<std::string> >("Assigned Regions",regs);
	    XMLString::release(&textContent2);
	    if (!compare_region_names(regs, def_list)) {
              Errors::Message msg;
              msg << "Amanzi::InputTranslator: ERROR - invalid region in Boundary Conditions Section - " ;
              msg << "valid regions are: \n" ;
              for (int r=0; r<regionNames_string_.size(); r++) {
                msg << "    " << regionNames_string_[r]<< "\n";
              }
              msg << "  Please correct and try again \n" ;
              Exceptions::amanzi_throw(msg);
	    }
          }
          else if (strcmp(tagName,"liquid_phase")==0) {
            //TODO: EIB - deal with liquid phase
            attrMap = BCNode->getAttributes();
            nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
	    if (nodeAttr) {
              textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
              phaseName = std::string(textContent2);
              if (phaseName=="water") {
                phaseName="Water";
              }
              XMLString::release(&textContent2);
	    }
            else {
              throw_error_missattr("boundary_conditions", "attribute", "name", "liquid_phase");
	    }

	    // loop over children, deal with liquid_component, solute_component, geomchemistry
            DOMNodeList* compList = BCNode->getChildNodes();
            for (int k=0; k<compList->getLength(); k++) {
              DOMNode* compNode = compList->item(k) ;
              char* compName  = XMLString::transcode(compNode->getNodeName());

              if (strcmp(compName,"liquid_component")==0) {
	        Teuchos::Array<double> vals;
	        Teuchos::Array<double> times;
	        Teuchos::Array<std::string> funcs;
	        std::string bcname;
	        std::string valname;
	        bool hasCoordsys = false;
                std::string coordsys("Absolute");
                bool hasSubmodel = false;
                std::string submodel;
                Teuchos::Array<double> grad;
                Teuchos::Array<double> ref_pt;
                double ref_val;

	        // loop over children and deal with bc's
                DOMNodeList* bcChildList = compNode->getChildNodes();
                for (int l=0; l<bcChildList->getLength(); l++) {
                  DOMNode* bcChildNode = bcChildList->item(l) ;
                  if (DOMNode::ELEMENT_NODE == bcChildNode->getNodeType()) {
                    char* bcChildName  = XMLString::transcode(bcChildNode->getNodeName());
		    std::string function;
		    double value;
		    double time;
                    bool success = true;
                    std::string xmlBCName;
                    std::stringstream errmsg;

		    if (strcmp(bcChildName,"seepage_face")==0){
                      xmlBCName = "seepage_face";
		      bcname = "BC: Seepage";
		      valname = "Inward Mass Flux";
                      DOMElement* bcElem = static_cast<DOMElement*>(bcChildNode);
                      if (bcElem->hasAttribute(XMLString::transcode("function"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("function")));
                        if (strlen(textContent2) > 0 ) {
                          if (strcmp(textContent2,"linear")==0) {
                            function = "Linear";
                          }
                          else if (strcmp(textContent2,"constant")==0) {
                            function = "Constant";
                          }
                          else if (strcmp(textContent2,"uniform")==0) {
                            function = "Uniform";
                          }
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'function' for BC seepage_face " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'function' for BC seepage_face " << std::endl ;
                      }
                      if (bcElem->hasAttribute(XMLString::transcode("inward_mass_flux"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("inward_mass_flux")));
                        if (strlen(textContent2) > 0 ) {
                          value = get_time_value(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'inward_mass_flux' for BC seepage_face " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'inward_mass_flux' for BC seepage_face " << std::endl ;
                      }
                      if (bcElem->hasAttribute(XMLString::transcode("start"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("start")));
                        if (strlen(textContent2) > 0 ) {
                          time = get_time_value(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'start' for BC seepage_face \n" ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'start' for BC seepage_face " << std::endl ;
                      }
		    }
                    else if (strcmp(bcChildName,"no_flow")==0) {
                      xmlBCName = "no_flow";
		      bcname = "BC: Zero Flow";
                      DOMElement* bcElem = static_cast<DOMElement*>(bcChildNode);
                      if (bcElem->hasAttribute(XMLString::transcode("function"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("function")));
                        if (strlen(textContent2) > 0 ) {
                          if (strcmp(textContent2,"linear")==0) {
                            function = "Linear";
                          }
                          else if (strcmp(textContent2,"constant")==0) {
                            function = "Constant";
                          }
                          else if (strcmp(textContent2,"uniform")==0) {
                            function = "Uniform";
                          }
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'function' for BC no_flow " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'function' for BC no_flow " << std::endl ;
                      }
                      if (bcElem->hasAttribute(XMLString::transcode("start"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("start")));
                        if (strlen(textContent2) > 0 ) {
                          time = get_time_value(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'start' for BC no_flow " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'start' for BC no_flow " << std::endl ;
                      }
                    }
                    else if (strcmp(bcChildName,"linear_hydrostatic")==0) {
                      xmlBCName = "linear_hydrostatic";
                      bcname = "BC: Linear Hydrostatic";
                      DOMElement* bcElem = static_cast<DOMElement*>(bcChildNode);
                      if (bcElem->hasAttribute(XMLString::transcode("gradient_value"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("gradient_value")));
                        if (strlen(textContent2) > 0 ) {
                          // translate to array
                          grad = make_coordinates(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'gradient_value' for BC linear_hydrostatic " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'gradient_value' for BC linear_hydrostatic " << std::endl ;
                      }
                      if (bcElem->hasAttribute(XMLString::transcode("reference_point"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("reference_point")));
                        if (strlen(textContent2) > 0 ) {
                          // translate to array
                          ref_pt = make_coordinates(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'reference_point' for BC linear_hydrostatic " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'reference_point' for BC linear_hydrostatic " << std::endl ;
                      }
                      if (bcElem->hasAttribute(XMLString::transcode("reference_water_table_height"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("reference_water_table_height")));
                        if (strlen(textContent2) > 0 ) {
                          ref_val = get_double_constant(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'reference_water_table_height' for BC linear_hydrostatic " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'reference_water_table_height' for BC linear_hydrostatic " << std::endl ;
                      }
                      if (bcElem->hasAttribute(XMLString::transcode("submodel"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("submodel")));
                        hasSubmodel = true;
                        if (strcmp(textContent2,"no_flow_above_water_table")==0){
                          submodel = "No Flow Above Water Table";
                        }
                        XMLString::release(&textContent2);
                      }
                    }
                    else if (strcmp(bcChildName,"linear_pressure")==0) {
                      xmlBCName = "linear_pressure";
                      bcname = "BC: Linear Pressure";
                      DOMElement* bcElem = static_cast<DOMElement*>(bcChildNode);
                      if (bcElem->hasAttribute(XMLString::transcode("gradient_value"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("gradient_value")));
                        if (strlen(textContent2) > 0 ) {
                          // translate to array
                          grad = make_coordinates(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'gradient_value' for BC linear_pressure " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'gradient_value' for BC linear_pressure " << std::endl ;
                      }
                      if (bcElem->hasAttribute(XMLString::transcode("reference_point"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("reference_point")));
                        if (strlen(textContent2) > 0 ) {
                          // translate to array
                          ref_pt = make_coordinates(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'reference_point' for BC linear_pressure " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'reference_point' for BC linear_pressure " << std::endl ;
                      }
                      if (bcElem->hasAttribute(XMLString::transcode("reference_value"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("reference_value")));
                        if (strlen(textContent2) > 0 ) {
                          ref_val = get_double_constant(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'reference_value' for BC linear_pressure " << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'reference_value' for BC linear_pressure " << std::endl ;
                      }
                    }
                    else {
                      DOMElement* bcElem = static_cast<DOMElement*>(bcChildNode);
                      // translate boundary condition name here
		      if (strcmp(bcChildName,"inward_mass_flux")==0) {
		        bcname = "BC: Flux";
		        valname = "Inward Mass Flux";
		      }
                      else if (strcmp(bcChildName,"inward_volumetric_flux")==0) {
		        bcname = "BC: Flux";
		        valname = "Outward Volumetric Flux";
                      }
                      else if (strcmp(bcChildName,"outward_mass_flux")==0) {
		        bcname = "BC: Flux";
		        valname = "Inward Mass Flux";
		      }
                      else if (strcmp(bcChildName,"outward_volumetric_flux")==0) {
		        bcname = "BC: Flux";
		        valname = "Outward Volumetric Flux";
		      }
                      else if (strcmp(bcChildName,"uniform_pressure")==0) {
		        bcname = "BC: Uniform Pressure";
		        valname = "Values";
		      }
                      else if (strcmp(bcChildName,"hydrostatic")==0) {
		        bcname = "BC: Hydrostatic";
		        valname = "Water Table Height";
		        //TODO: EIB - update if necessary
		        if (bcElem->hasAttribute(XMLString::transcode("submodel"))) {
                          textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("submodel")));
		          hasSubmodel = true;
		          if (strcmp(textContent2,"no_flow_above_water_table")==0){
			    submodel = "No Flow Above Water Table";
		          }
                          XMLString::release(&textContent2);
		        }
                        if (bcElem->hasAttribute(XMLString::transcode("coordinate_system"))) {
                          textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("coordinate_system")));
                          hasCoordsys = true;
                          if (strcmp(textContent2,"relative to mesh top")==0){
                            coordsys = "Relative";
                          }
                          XMLString::release(&textContent2);
                        }
		      }
                      // get function type
                      if (bcElem->hasAttribute(XMLString::transcode("function"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("function")));
                        if (strlen(textContent2) > 0 ) {
                          if (strcmp(textContent2,"linear")==0) {
                            function = "Linear";
                          }
                          else if (strcmp(textContent2,"constant")==0) {
                            function = "Constant";
                          }
                          else if (strcmp(textContent2,"uniform")==0) {
                            function = "Uniform";
                          }
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'function' for BC " << bcChildName << std::endl ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'function' for BC " << bcChildName << std::endl ;
                      }
                      // get value
                      if (bcElem->hasAttribute(XMLString::transcode("value"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("value")));
                        if (strlen(textContent2) > 0 ) {
                          value = get_double_constant(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'value' for BC " << bcChildName << " \n" ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'value' for BC " << bcChildName << std::endl ;
                      }
                      // get start time
                      if (bcElem->hasAttribute(XMLString::transcode("start"))) {
                        textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("start")));
                        if (strlen(textContent2) > 0 ) {
                          time = get_time_value(textContent2, def_list);
                        }
                        else {
                          success = false;
                          errmsg << "  Ill-formed 'start' for BC " << bcChildName << " \n" ;
                        }
                        XMLString::release(&textContent2);
                      }
                      else {
                        success = false;
                        errmsg << "  Missing 'start' for BC " << bcChildName << std::endl ;
                      }
		    }
                    
                    // Throw accumulated exceptions
                    if (!success) {
                      Errors::Message msg;
                      msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing boundary_conditions - " ;
                      msg << errmsg.str().c_str();
                      msg << "  Please correct and try again \n" ;
                      Exceptions::amanzi_throw(msg);
                    }

		    // put time, function, value in sorted order
		    if (times.length() == 0) {               // if first time through
                      times.append(time);
		      funcs.append(function);
                      vals.append(value);
		    }
                    else {
		      if (time >= times[times.length()-1]) { // if already sorted
		        times.append(time);
		        funcs.append(function);
		        vals.append(value);
		      }
                      else {                              // otherwise, sort
                        int here(0);
                        int idx(0);
                        // figure out the correct position (times already in order from previous steps)
                        while (times[idx] < time) {
                          idx++;
                          here = idx;
                        }
                        // append last values to end
                        times.append(times[times.length()-1]);
                        vals.append(vals[vals.length()-1]);
                        funcs.append(funcs[funcs.length()-1]);
                        // copy entries over one
                        for (int i=times.length()-1; i>here; i--) {
                          times[i] = times[i-1];
                          vals[i] = vals[i-1];
                          funcs[i] = funcs[i-1];
                        }
                        // add new entries "here"
                        times[here] = time;
                        vals[here] = value;
                        funcs[here] = function;
                      }
		    }
		  }
	        }
	        // if len array == 1: add dummy vals to create and interval (only for unstructured)
	        if (times.length()==1 && bcname != "BC: Zero Flow" && isUnstr_){
		  // EIB: thought this would work, simple version doesn't in the steady case where end time = 0
		  //      overkill to do correctly
                  //if (def_list.sublist("simulation").isParameter("simulation_end")) {
		  //  double end_time = (def_list.sublist("simulation").get<double>("simulation_end");
		  //  if (end_time > (times[0])) {
	 	  //    times.append(end_time);
		  //  } else {
		  //    times.append(times[0]+1.);
		  //  }
		  //}
		  //times.append(times[0]+1.);
		  if (def_list.sublist("simulation").isParameter("simulation_end")) {
	  	      times.append(def_list.sublist("simulation").get<double>("simulation_end")+1.);
		  } else {
	  	      times.append(times[0]+1.);
		  }
		  vals.append(vals[0]);
	        }
	        //EIB - this is iffy!!! Talked with Ellen, this is consistent with her assumptions in Akuna, for now
	        if (times.length()==funcs.length()) funcs.remove(funcs.length()-1);

	        // create a BC new list here
                Teuchos::ParameterList newbclist;
                if (bcname == "BC: Linear Hydrostatic") {
                  newbclist.set<Teuchos::Array<double> >("Gradient Value",grad);
                  newbclist.set<Teuchos::Array<double> >("Reference Point",ref_pt);
                  newbclist.set<double>("Reference Water Table Height",ref_val);
                  if (hasSubmodel) newbclist.set<std::string>("Submodel",submodel);
                }
                else if (bcname == "BC: Linear Pressure") {
                  newbclist.set<Teuchos::Array<double> >("Gradient Value",grad);
                  newbclist.set<Teuchos::Array<double> >("Reference Point",ref_pt);
                  newbclist.set<double>("Reference Value",ref_val);
                  
                }
                else {
	          newbclist.set<Teuchos::Array<double> >("Times",times);
	          newbclist.set<Teuchos::Array<std::string> >("Time Functions",funcs);
	          if (bcname != "BC: Zero Flow") newbclist.set<Teuchos::Array<double> >(valname,vals);
	          if (bcname == "BC: Hydrostatic" && hasCoordsys) newbclist.set<std::string>("Coordinate System",coordsys);
                  if (bcname == "BC: Hydrostatic" && hasSubmodel) newbclist.set<std::string>("Submodel",submodel);

                }
                bclist.sublist(bcname) = newbclist;
	      }
              if (strcmp(compName,"solute_component")==0) {
	        char* solName;
	        //Teuchos::ParameterList sclist;
	        // loop over elements to build time series, add to list
	        Teuchos::Array<double> vals;
	        Teuchos::Array<double> times;
	        Teuchos::Array<std::string> funcs;

	        // EIB - loop over elements, aqueous_conc, to get all. Could have multiple time steps for multiple components
	        Teuchos::ParameterList sc_tmplist;
                DOMNodeList* acList = compNode->getChildNodes();
                for (int l=0; l<acList->getLength(); l++) {
                  DOMNode* cur = acList->item(l) ;
                  if (DOMNode::ELEMENT_NODE == cur->getNodeType()) {
		    std::string function;
		    double value;
		    double time;
                    DOMElement* bcElem = static_cast<DOMElement*>(cur);
		    if (strcmp(XMLString::transcode(bcElem->getTagName()),"aqueous_conc") != 0) {
                      throw_error_illformed("boundary_conditions->solute_component", "element", "aqueous_conc");
		    }
                    solName = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("name")));
		    if (solName[0] == '\0'){
                      throw_error_illformed("boundary_conditions->solute_component", "name", "aqueous_conc");
		    }
                    textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("function")));
		    if (textContent2[0] == '\0'){
                      throw_error_illformed("boundary_conditions->solute_component", "function", "aqueous_conc");
		    }
		    if (strcmp(textContent2,"linear")==0) {
		      function = "Linear";
		    }
                    else if (strcmp(textContent2,"constant")==0) {
		      function = "Constant";
		      // EIB - uniform is a space option, not a time option
		    }
                    textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("value")));
		    if (textContent2[0] == '\0'){
                      throw_error_illformed("boundary_conditions->solute_component", "value", "aqueous_conc");
		    }
                    value = get_time_value(textContent2, def_list);
                    textContent2 = XMLString::transcode(bcElem->getAttribute(XMLString::transcode("start")));
		    if (textContent2[0] == '\0'){
                      throw_error_illformed("boundary_conditions->solute_component", "start", "aqueous_conc");
		    }
                    time = get_time_value(textContent2, def_list);
		    if ( sc_tmplist.isParameter(solName)) {
		      Teuchos::Array<double> tmpV;
		      Teuchos::Array<double> tmpT;
		      Teuchos::Array<std::string> tmpF;
		      tmpV = sc_tmplist.sublist(solName).get<Teuchos::Array<double> >("value");
		      tmpV.append(value);
		      sc_tmplist.sublist(solName).set<Teuchos::Array<double> >("value",tmpV);
		      tmpT = sc_tmplist.sublist(solName).get<Teuchos::Array<double> >("time");
		      tmpT.append(time);
		      sc_tmplist.sublist(solName).set<Teuchos::Array<double> >("time",tmpT);
		      tmpF = sc_tmplist.sublist(solName).get<Teuchos::Array<std::string> >("function");
		      tmpF.append(function);
		      sc_tmplist.sublist(solName).set<Teuchos::Array<std::string> >("function",tmpF);
		    }
                    else {
		      Teuchos::ParameterList tmpsc_list;
		      Teuchos::Array<double> tmpV;
		      Teuchos::Array<double> tmpT;
		      Teuchos::Array<std::string> tmpF;
		      tmpV.append(value);
		      tmpsc_list.set<Teuchos::Array<double> >("value",tmpV);
		      tmpT.append(time);
		      tmpsc_list.set<Teuchos::Array<double> >("time",tmpT);
		      tmpF.append(function);
		      tmpsc_list.set<Teuchos::Array<std::string> >("function",tmpF);
		      sc_tmplist.sublist(solName) = tmpsc_list;
		    }
		  }
	        }
	        // EIB - now do time sorting, if have multiple time steps for each component
                for (Teuchos::ParameterList::ConstIterator i = sc_tmplist.begin(); i != sc_tmplist.end(); i++) {
                  Teuchos::ParameterList& curComp_list = sc_tmplist.sublist(sc_tmplist.name(i)) ;
                  Teuchos::Array<double> times = curComp_list.get<Teuchos::Array<double> >("time") ;
                  Teuchos::Array<double> values = curComp_list.get<Teuchos::Array<double> >("value") ;
                  Teuchos::Array<std::string> funcs = curComp_list.get<Teuchos::Array<std::string> >("function") ;
                  Teuchos::Array<double> sort_vals;
                  Teuchos::Array<double> sort_times;
		  Teuchos::Array<std::string> sort_func;
                  for (int j=0; j<times.length(); j++) {
                    if (j==0) {
                      sort_times.append(times[j]);
                      sort_func.append(funcs[j]);
                      sort_vals.append(values[j]);
                    }
                    else {
                      if (times[j] >= sort_times[sort_times.length()-1]) { //if already sorted
                        sort_times.append(times[j]);
                        sort_func.append(funcs[j]);
                        sort_vals.append(values[j]);
                      }
                      else {                                            // otherwise sort
                        // otherwise, sort
                        int here(0);
                        int idx(0);
                        // figure out the correct position (times already in order from previous steps)
                        while (sort_times[idx] < times[j]) {
                          idx++;
                          here = idx;
                        }
                        // append last values to end
                        sort_times.append(sort_times[sort_times.length()-1]);
                        sort_vals.append(sort_vals[sort_vals.length()-1]);
                        sort_func.append(sort_func[sort_func.length()-1]);
                        // copy entries over one
                        for (int i=sort_times.length()-1; i>here; i--) {
                          sort_times[i] = sort_times[i-1];
                          sort_vals[i] = sort_vals[i-1];
                          sort_func[i] = sort_func[i-1];
                        }
                        // add new entries "here"
                        sort_times[here] = times[j];
                        sort_vals[here] = values[j];
                        sort_func[here] = funcs[j];

                      }
                    }
		  }
                  // if len array == 1: add dummy vals to create and interval
                  if (sort_times.length()==1){
                    //sort_times.append(sort_times[0]+1.);
                    if (def_list.sublist("simulation").isParameter("simulation_end")) {
                      sort_times.append(def_list.sublist("simulation").get<double>("simulation_end")+1.);
                    }
                    else {
                      sort_times.append(times[0]+1.);
                    }
                    sort_vals.append(sort_vals[0]);
	          }
                  //EIB - this is iffy!!! Talked with Ellen, this is consistent with her assumptions in Akuna, for now
                  if (sort_times.length()==sort_func.length()) sort_func.remove(sort_func.length()-1);
                  // EIB - add these new sorted arrays back to PL
                  sc_tmplist.sublist(sc_tmplist.name(i)).set<Teuchos::Array<double> >("sorted_times",sort_times);
                  sc_tmplist.sublist(sc_tmplist.name(i)).set<Teuchos::Array<double> >("sort_values",sort_vals);
                  sc_tmplist.sublist(sc_tmplist.name(i)).set<Teuchos::Array<std::string> >("sorted_functions",sort_func);
                }
	        //TODO: EIB - not added concerntation units, need to grab from units
                // EIB - now add each solute BC to PL
                for (Teuchos::ParameterList::ConstIterator i = sc_tmplist.begin(); i != sc_tmplist.end(); i++) {
                  Teuchos::ParameterList sclist;
                  Teuchos::ParameterList& curComp_list = sc_tmplist.sublist(sc_tmplist.name(i)) ;
	          sclist.sublist("BC: Uniform Concentration").set<Teuchos::Array<double> >("Times",curComp_list.get<Teuchos::Array<double> >("sorted_times"));
	          sclist.sublist("BC: Uniform Concentration").set<Teuchos::Array<std::string> >("Time Functions",curComp_list.get<Teuchos::Array<std::string> >("sorted_functions"));
	          sclist.sublist("BC: Uniform Concentration").set<Teuchos::Array<double> >("Values",curComp_list.get<Teuchos::Array<double> >("sort_values"));
                  bclist.sublist("Solute BC").sublist("Aqueous").sublist(phaseName).sublist(sc_tmplist.name(i)) = sclist;
                }
	        XMLString::release(&solName);
	      }
              else if (strcmp(compName,"geochemistry")==0) {

		const Teuchos::Array<std::string>& soluteNames = def_list.get<Teuchos::Array<std::string> >("solutes");

		if (soluteNames.size() > 0) {
		  char* funcType;
		  Teuchos::ParameterList sbclist;
		  attrMap = compNode->getAttributes();
		  nodeAttr = attrMap->getNamedItem(XMLString::transcode("function"));
		  if (nodeAttr) {
		    funcType = XMLString::transcode(nodeAttr->getNodeValue());
		  }
		  else {
		    throw_error_missattr("boundary_conditions", "attribute", "function", "geochemistry");
		  }
		  if (strcmp(funcType,"uniform")==0) {
		    nodeAttr = attrMap->getNamedItem(XMLString::transcode("constraint"));
		    if (nodeAttr) {
		      textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
		      for (int i=0; i<soluteNames.size(); ++i) {
			Teuchos::ParameterList sbc, sbc1;
			sbc.set<std::string>("Geochemical Condition",trim_string(textContent2));
			sbc1.sublist("BC: Uniform Concentration") = sbc;
			sbclist.sublist(soluteNames[i]) = sbc1;
		      }
		      XMLString::release(&textContent2);
		    }
		    else {
		      throw_error_missattr("boundary_conditions", "attribute", "value", "solute_component");
		    }
		  }
		  else if (strcmp(funcType,"linear")==0) {
		    // TODO: EIB - currently can't handle this
		  }
		  //TODO: EIB - not added concerntation units, confused by what to add. grab from units?
		  bclist.sublist("Solute BC").sublist("Aqueous").sublist(phaseName) = sbclist;
		  XMLString::release(&funcType);
		}
	      }
	    }
          } else if (strcmp(tagName,"solid_phase")==0) {
            //TODO: EIB - deal with solid phase -> mineral, geochemisty
          }
          XMLString::release(&tagName);
        }
        list.sublist(textContent) = bclist;
        XMLString::release(&textContent);
      }
    }
  }
  return list;
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_sources(DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  Teuchos::ParameterList list;

  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Sources" << std::endl;
  }

  // get Sources node
  DOMNodeList* nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("sources"));
  if (nodeList->getLength() > 0) {
  DOMNode* nodeSC = nodeList->item(0);
  DOMElement* elementSC = static_cast<DOMElement*>(nodeSC);

  // get list of Source
  DOMNodeList* SCList = elementSC->getElementsByTagName(XMLString::transcode("source"));
  std::string phase("Aqueous"); //NOTE:: EIB: currently only support this, add checks later
  std::string component("Water");

  for (int i=0; i<SCList->getLength(); i++) {
    DOMNode* cur = SCList->item(i) ;

    if (DOMNode::ELEMENT_NODE == cur->getNodeType()) {
      char* textContent;
      DOMNamedNodeMap* attrMap = cur->getAttributes();
      DOMNode* nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
      if (nodeAttr) {
        textContent = XMLString::transcode(nodeAttr->getNodeValue());
      }
      else {
        throw_error_missattr("sources", "attribute", "name", "source");
      }

      Teuchos::ParameterList sclist;
      DOMNodeList* SC = cur->getChildNodes();
      for (int j=0; j<SC->getLength(); j++) {
        DOMNode* SCNode = SC->item(j) ;

        if (DOMNode::ELEMENT_NODE == SCNode->getNodeType()) {
          char* tagName  = XMLString::transcode(SCNode->getNodeName());
          if (strcmp(tagName,"assigned_regions")==0) {
            char* textContent2 = XMLString::transcode(SCNode->getTextContent());
	    Teuchos::Array<std::string> regs = make_regions_list(textContent2);
	    sclist.set<Teuchos::Array<std::string> >("Assigned Regions",regs);
	    XMLString::release(&textContent2);
	    if (!compare_region_names(regs, def_list)) {
              Errors::Message msg;
              msg << "Amanzi::InputTranslator: ERROR - invalid region in Sources Section - " ;
              msg << "valid regions are: \n" ;
              for (int r=0; r<regionNames_string_.size(); r++) {
                msg << "    " << regionNames_string_[r]<< "\n";
              }
              msg << "  Please correct and try again \n" ;
              Exceptions::amanzi_throw(msg);
	    }
          }
          else if (strcmp(tagName,"liquid_phase")==0) {
            attrMap = SCNode->getAttributes();
            nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
            std::string phaseName;
	    if (nodeAttr) {
              char* textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
              phaseName = std::string(textContent2);
              if (phaseName=="water") {
                phaseName="Water";
              }
              XMLString::release(&textContent2);
            }
            else {
              throw_error_missattr("sources", "attribute", "name", "liquid_phase");
            }

            DOMNodeList* compList = SCNode->getChildNodes();
            for (int k=0; k<compList->getLength(); k++) {
	      DOMNode* compNode = compList->item(k) ;
              char* compName  = XMLString::transcode(compNode->getNodeName());
              DOMNamedNodeMap* attrMap2 = compNode->getAttributes();
              if (strcmp(compName,"liquid_component")==0) {
	        Teuchos::Array<double> vals;
	        Teuchos::Array<double> times;
	        Teuchos::Array<std::string> funcs;
	        std::string scname;
		// get component name
                attrMap = compNode->getAttributes();
                nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
                char* compName2 ;
		if (nodeAttr) {
                  compName2 = XMLString::transcode(nodeAttr->getNodeValue());
                  //if (std::string(compName2) != "Water" || std::string(compName2) != "water") {
		  //  component = std::string(compName2);
                  //}
                }
                else {
                  throw_error_missattr("sources", "attribute", "name", "liquid_component");
                }

		// loop over children
                DOMNodeList* scChildList = compNode->getChildNodes();
                for (int l=0; l<scChildList->getLength(); l++) {
                  DOMNode* scChildNode = scChildList->item(l) ;
                  if (DOMNode::ELEMENT_NODE == scChildNode->getNodeType()) {
                    char* scChildName  = XMLString::transcode(scChildNode->getNodeName());
		    if (strcmp(scChildName,"volume_weighted")==0){
		     scname = "Source: Volume Weighted";
		    }
                    //else if (strcmp(scChildName,"permeability_weighted")==0){
                    else if (strcmp(scChildName,"perm_weighted")==0){
		     scname = "Source: Permeability Weighted";
		    }
                    
		    // loop over any attributes that may exist
                    DOMNamedNodeMap* attrMap2 = scChildNode->getAttributes();
                    std::string function;
                    double value;
                    double time;
                    for (int l=0; l<attrMap2->getLength(); l++) {
                      DOMNode* attrNode = attrMap2->item(l) ;
                      if (DOMNode::ATTRIBUTE_NODE == attrNode->getNodeType()) {
                        
                        char* attrName = XMLString::transcode(attrNode->getNodeName());
                        char* attrValue = XMLString::transcode(attrNode->getNodeValue());
                        if (strcmp(attrName,"function")==0) {
		          if (strcmp(attrValue,"linear")==0) {
		            function = "Linear";
		          }
                          else if (strcmp(attrValue,"constant")==0) {
		            function = "Constant";
		          }
                          else if (strcmp(attrValue,"uniform")==0) {
		            function = "Uniform";
		          }
			}
                        else if (strcmp(attrName,"start")==0) {
		          time = get_time_value(attrValue, def_list);
			}
                        else if (strcmp(attrName,"value")==0) {
		          value = get_double_constant(attrValue, def_list);
                        }
                        XMLString::release(&attrName);
                        XMLString::release(&attrValue);
                      }
		    }
                    
                    // put time, function, value in sorted order
                    if (times.length() == 0) {               // if first time through
                      times.append(time);
                      funcs.append(function);
                      vals.append(value);
                    }
                    else {
                      if (time >= times[times.length()-1]) { // if already sorted
                        times.append(time);
                        funcs.append(function);
                        vals.append(value);
                      }
                      else {                              // otherwise, sort
                        int here(0);
                        int idx(0);
                        // figure out the correct position (times already in order from previous steps)
                        while (times[idx] < time) {
                          idx++;
                          here = idx;
                        }
                        // append last values to end
                        times.append(times[times.length()-1]);
                        vals.append(vals[vals.length()-1]);
                        funcs.append(funcs[funcs.length()-1]);
                        // copy entries over one
                        for (int m=times.length()-1; m>here; m--) {
                          times[m] = times[m-1];
                          vals[m] = vals[m-1];
                          funcs[m] = funcs[m-1];
                        }
                        // add new entries "here"
                        times[here] = time;
                        vals[here] = value;
                        funcs[here] = function;
                      }
                    }
		  }
		}
	        if (times.length()==1 ){
	  	  //times.append(times[0]+1.);
		  if (def_list.sublist("simulation").isParameter("simulation_end")) {
	  	    times.append(def_list.sublist("simulation").get<double>("simulation_end")+1.);
		  }
                  else {
	  	    times.append(times[0]+1.);
		  }
		  vals.append(vals[0]);
	        }
	        //EIB - this is iffy!!! Talked with Ellen, this is consistent with her assumptions in Akuna, for now
	        if (times.length()==funcs.length() && funcs.length()>0) funcs.remove(funcs.length()-1); 
	        Teuchos::ParameterList newsclist;
		if (times.length() > 0) {
	          newsclist.set<Teuchos::Array<double> >("Times",times);
	          newsclist.set<Teuchos::Array<std::string> >("Time Functions",funcs);
	          newsclist.set<Teuchos::Array<double> >("Values",vals);
		}
	        sclist.sublist(scname) = newsclist;
	      }
              else if (strcmp(compName,"solute_component")==0) {
	        Teuchos::Array<double> vals;
	        Teuchos::Array<double> times;
	        Teuchos::Array<std::string> funcs;
	        std::string scname;
	        std::string soluteName;
                bool isReleaseModel = false;
                Teuchos::ParameterList newsclist;
                DOMNodeList* scChildList = compNode->getChildNodes();
                for (int l=0; l<scChildList->getLength(); l++) {
                  DOMNode* scChildNode = scChildList->item(l) ;
                  if (DOMNode::ELEMENT_NODE == scChildNode->getNodeType()) {
                    char* scChildName  = XMLString::transcode(scChildNode->getNodeName());
		    if (strcmp(scChildName,"flow_weighted_conc")==0){
		     scname = "Source: Flow Weighted Concentration";
		    }
                    else if (strcmp(scChildName,"uniform_conc")==0){
		     scname = "Source: Uniform Concentration";
                    }
                    else if (strcmp(scChildName,"diffusion_dominated_release")==0){
                      scname = "Source: Diffusion Dominated Release Model";
                      isReleaseModel = true;
                    }
		    
                    // loop over any attributes that may exist
                    DOMNamedNodeMap* attrMap2 = scChildNode->getAttributes();
                    std::string function;
                    double value;
                    double time;
                    for (int l=0; l<attrMap2->getLength(); l++) {
                      DOMNode* attrNode = attrMap2->item(l) ;
                      if (DOMNode::ATTRIBUTE_NODE == attrNode->getNodeType()) {
                        char* attrName = XMLString::transcode(attrNode->getNodeName());
                        char* attrValue = XMLString::transcode(attrNode->getNodeValue());
			if (strcmp(attrName,"function")==0) {
		          if (strcmp(attrValue,"linear")==0) {
		            function = "Linear";
		          }
                          else if (strcmp(attrValue,"constant")==0) {
		            function = "Constant";
		          }
                          else if (strcmp(attrValue,"uniform")==0) {
		            function = "Uniform";
		          }
			}
                        else if (strcmp(attrName,"start")==0) {
		          time = get_time_value(attrValue, def_list);
			}
                        else if (strcmp(attrName,"value")==0) {
		          value = get_double_constant(attrValue, def_list);
			}
                        else if (strcmp(attrName,"name")==0) {
		          soluteName = attrValue;
			}
                        else if (strcmp(attrName,"total_inventory")==0) {
                          newsclist.set<double>("Total Inventory",atof(attrValue));
                        }
                        else if (strcmp(attrName,"mixing_length")==0) {
                          newsclist.set<double>("Mixing Length",atof(attrValue));
                        }
                        else if (strcmp(attrName,"effective_diffusion_coefficient")==0) {
                          newsclist.set<double>("Effective Diffusion Coefficient",atof(attrValue));
                        }
                        XMLString::release(&attrName);
                        XMLString::release(&attrValue);
		      }
                    }
                    
                    // put time, function, value in sorted order
                    if (times.length() == 0) {               // if first time through
                      times.append(time);
                      funcs.append(function);
                      vals.append(value);
                    }
                    else {
                      if (time >= times[times.length()-1]) { // if already sorted
                        times.append(time);
                        funcs.append(function);
                        vals.append(value);
                      }
                      else {                              // otherwise, sort
                        int here(0);
                        int idx(0);
                        // figure out the correct position (times already in order from previous steps)
                        while (times[idx] < time) {
                          idx++;
                          here = idx;
                        }
                        // append last values to end
                        times.append(times[times.length()-1]);
                        vals.append(vals[vals.length()-1]);
                        funcs.append(funcs[funcs.length()-1]);
                        // copy entries over one
                        for (int m=times.length()-1; m>here; m--) {
                          times[m] = times[m-1];
                          vals[m] = vals[m-1];
                          funcs[m] = funcs[m-1];
                        }
                        // add new entries "here"
                        times[here] = time;
                        vals[here] = value;
                        funcs[here] = function;
                      }
                    }
		  }
		}
	        if (times.length()==1 ){
		  if (def_list.sublist("simulation").isParameter("simulation_end")) {
	  	    times.append(def_list.sublist("simulation").get<double>("simulation_end")+1.);
		  }
                  else {
	  	    times.append(times[0]+1.);
		  }
		  vals.append(vals[0]);
	        }
	        //EIB - this is iffy!!! Talked with Ellen, this is consistent with her assumptions in Akuna, for now
	        if (times.length()==funcs.length() && funcs.length()>0) funcs.remove(funcs.length()-1);
                newsclist.set<Teuchos::Array<double> >("Times",times);
                if (!isReleaseModel) {
                  newsclist.set<Teuchos::Array<std::string> >("Time Functions",funcs);
                  newsclist.set<Teuchos::Array<double> >("Values",vals);

                }
                //sclist.sublist(scname) = newsclist;
                sclist.sublist("Solute SOURCE").sublist(phase).sublist(component).sublist(soluteName).sublist(scname) = newsclist;
                //sclist.sublist("Solute SOURCE").sublist(phase).sublist(soluteName).sublist(scname) = newsclist;
	      }
	    }
	  }
	}
      }
      list.sublist(textContent) = sclist;
      XMLString::release(&textContent);
    }
  }
  }

  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_output(DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  Teuchos::ParameterList list;

  DOMNodeList* nodeList;
  DOMNodeList* visList;
  DOMNodeList* chkList;
  DOMNodeList* obsList;
  DOMNodeList* tmpList;
  DOMNamedNodeMap* attrMap;
  DOMNode* tmpNode;
  DOMNode* nodeAttr;
  char* textContent;
  char* textContent2;


  if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *voI_->os() << "Getting Macros" << std::endl;
  }

  // get definitions node - this node MAY exist ONCE
  // this contains any time macros and cycle macros
  // they are stored in the outputs of the old format
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("definitions"));
  if (nodeList->getLength() > 0) {
    DOMNode* defNode = nodeList->item(0);
    DOMElement* defElement = static_cast<DOMElement*>(defNode);
    DOMNodeList* macroList = xmlDoc->getElementsByTagName(XMLString::transcode("macros"));
    Teuchos::ParameterList tmPL;
    Teuchos::ParameterList cmPL;
    //loop over children
    DOMNodeList* children = macroList->item(0)->getChildNodes();
    for (int i=0; i<children->getLength(); i++) {
      DOMNode* currentNode = children->item(i) ;
      if (DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
	char* tagname = XMLString::transcode(currentNode->getNodeName());
	if (strcmp(tagname,"time_macro")==0) {
          Teuchos::ParameterList tm_parameter;
          attrMap = currentNode->getAttributes();
          nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
	  if (nodeAttr) {
            textContent = XMLString::transcode(nodeAttr->getNodeValue());
	  }
          else {
            throw_error_missattr("definitions", "attribute", "name", "time_macro");
	  }

	  // deal differently if "times" or "start-inter-stop"
          DOMNodeList* childList = currentNode->getChildNodes();
	  bool isTime = false;
          for (int j=0; j<childList->getLength(); j++) {
            DOMNode* timeNode = childList->item(j) ;
            if (DOMNode::ELEMENT_NODE == timeNode->getNodeType()) {
	      if (strcmp(XMLString::transcode(timeNode->getNodeName()),"time")==0)
		      isTime = true;
	    }
	  }
	  if ( isTime ) {
            Teuchos::Array<double> times;
            for (int j=0; j<childList->getLength(); j++) {
              DOMNode* timeNode = childList->item(j) ;
              if (DOMNode::ELEMENT_NODE == timeNode->getNodeType()) {
	        char* nodeTxt = XMLString::transcode(timeNode->getTextContent());
	        times.append(get_time_value(nodeTxt,def_list));
	        XMLString::release(&nodeTxt);
	      }
	    }
	    tm_parameter.set<Teuchos::Array<double> >("Values", times);
	  }
          else {
            DOMElement* curElement = static_cast<DOMElement*>(currentNode);
            DOMNodeList* curList = curElement->getElementsByTagName(XMLString::transcode("start"));
            tmpNode = curList->item(0);
	    char* nodeTxt = XMLString::transcode(tmpNode->getTextContent());
            Teuchos::Array<double> sps;
	    sps.append(get_time_value(nodeTxt,def_list));
	    XMLString::release(&nodeTxt);
            curList = curElement->getElementsByTagName(XMLString::transcode("timestep_interval"));
	    if (curList->getLength() >0) {
              tmpNode = curList->item(0);
	      nodeTxt = XMLString::transcode(tmpNode->getTextContent());
	      sps.append(get_time_value(nodeTxt,def_list));
	      XMLString::release(&nodeTxt);
              curList = curElement->getElementsByTagName(XMLString::transcode("stop"));
	      if (curList->getLength() >0) {
                tmpNode = curList->item(0);
	        nodeTxt = XMLString::transcode(tmpNode->getTextContent());
	        sps.append(get_time_value(nodeTxt,def_list));
	        XMLString::release(&nodeTxt);
	      }
              else {
	        sps.append(-1.0);
	      }
	      tm_parameter.set<Teuchos::Array<double> >("Start_Period_Stop", sps);
	    }
            else {
	      tm_parameter.set<Teuchos::Array<double> >("Values", sps);
	    }
	  }
	  tmPL.sublist(textContent) = tm_parameter;
	  XMLString::release(&textContent);
	}
        else if (strcmp(tagname,"cycle_macro")==0) {
          Teuchos::ParameterList cm_parameter;
          attrMap = currentNode->getAttributes();
          nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
	  if (nodeAttr) {
            textContent = XMLString::transcode(nodeAttr->getNodeValue());
	  }
          else {
            throw_error_missattr("definitions", "attribute", "name", "cycle_macro");
	  }

          DOMElement* curElement = static_cast<DOMElement*>(currentNode);
          DOMNodeList* curList = curElement->getElementsByTagName(XMLString::transcode("start"));
          tmpNode = curList->item(0);
	  char* nodeTxt = XMLString::transcode(tmpNode->getTextContent());
          Teuchos::Array<int> sps;
	  sps.append(get_int_constant(nodeTxt,def_list));
	  XMLString::release(&nodeTxt);
          curList = curElement->getElementsByTagName(XMLString::transcode("timestep_interval"));
	  if (curList->getLength() >0) {
            tmpNode = curList->item(0);
	    nodeTxt = XMLString::transcode(tmpNode->getTextContent());
	    sps.append(get_int_constant(nodeTxt,def_list));
	    XMLString::release(&nodeTxt);
            curList = curElement->getElementsByTagName(XMLString::transcode("stop"));
	    if (curList->getLength() >0) {
              tmpNode = curList->item(0);
	      nodeTxt = XMLString::transcode(tmpNode->getTextContent());
	      sps.append(get_int_constant(nodeTxt,def_list));
	      XMLString::release(&nodeTxt);
	    }
            else {
	      sps.append(-1.0);
	    }
	    cm_parameter.set<Teuchos::Array<int> >("Start_Period_Stop", sps);
	  }
          else {
	    cm_parameter.set<Teuchos::Array<int> >("Values", sps);
	  }
	  cmPL.sublist(textContent) = cm_parameter;
	  XMLString::release(&textContent);
	}
      }
    }
    list.sublist("Time Macros") = tmPL;
    list.sublist("Cycle Macros") = cmPL;
  }

  // get output node - this node must exist ONCE
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("output"));
  DOMNode* outNode = nodeList->item(0);
  DOMElement* outElement = static_cast<DOMElement*>(outNode);
  if (DOMNode::ELEMENT_NODE == outNode->getNodeType()) {
    DOMNodeList* outChildList = outNode->getChildNodes();
    for (int m=0; m<outChildList->getLength(); m++) {
      DOMNode* curoutNode = outChildList->item(m) ;
      if (DOMNode::ELEMENT_NODE == curoutNode->getNodeType()) {
        char* outName = XMLString::transcode(curoutNode->getNodeName());
        if (strcmp(outName,"vis")==0) {
          if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
            *voI_->os() << "Getting Output: Vis" << std::endl;
          }
          // get list of vis - this node MAY exist ONCE
          DOMNodeList* childList = curoutNode->getChildNodes();
          Teuchos::ParameterList visPL;
          for (int j=0; j<childList->getLength(); j++) {
            DOMNode* curKid = childList->item(j) ;
            textContent  = XMLString::transcode(curKid->getNodeName());
            if (strcmp(textContent,"base_filename")==0) {
	      textContent2 = XMLString::transcode(curKid->getTextContent());
              visPL.set<std::string>("File Name Base",trim_string(textContent2));
              XMLString::release(&textContent2);
	    }
            else if (strcmp(textContent,"num_digits")==0) {
	      textContent2 = XMLString::transcode(curKid->getTextContent());
              visPL.set<int>("File Name Digits",get_int_constant(textContent2,def_list));
              XMLString::release(&textContent2);
	    }
            else if (strcmp(textContent,"time_macros")==0 || strcmp(textContent,"time_macro")==0) {
	      textContent2 = XMLString::transcode(curKid->getTextContent());
	      Teuchos::Array<std::string> macro = make_regions_list(textContent2);
              visPL.set<Teuchos::Array<std::string> >("Time Macros",macro);
              XMLString::release(&textContent2);
	    }
            else if (strcmp(textContent,"cycle_macros")==0 || strcmp(textContent,"cycle_macro")==0) {
	      textContent2 = XMLString::transcode(curKid->getTextContent());
	      Teuchos::Array<std::string> macro = make_regions_list(textContent2);
              visPL.set<Teuchos::Array<std::string> >("Cycle Macros",macro);
              XMLString::release(&textContent2);
	    }
            else if (strcmp(textContent,"write_regions")==0) {

	      DOMNodeList* wrChildList = curKid->getChildNodes();
	      int numChild = wrChildList->getLength();
	      if (numChild > 0) {
		Teuchos::ParameterList wrPL;
		for (int L=0; L<numChild; L++) {
		  DOMNode* wrChildNode = wrChildList->item(L);
		  char* name = XMLString::transcode(wrChildNode->getNodeName());
		  if (strcmp(name,"field")==0) {
		    DOMElement* wrElem = static_cast<DOMElement*>(wrChildNode);
		    char *nameStr, *regStr;
		    if (wrElem->hasAttribute(XMLString::transcode("name"))) {
		      nameStr = XMLString::transcode(wrElem->getAttribute(XMLString::transcode("name")));
		    } else {
		      Exceptions::amanzi_throw("write_regions elements each require a \"name\" attribute");
		    }
		    if (wrElem->hasAttribute(XMLString::transcode("regions"))) {
		      regStr = XMLString::transcode(wrElem->getAttribute(XMLString::transcode("regions")));
		    } else {
		      Exceptions::amanzi_throw("write_regions elements each require a \"regions\" attribute");
		    }

		    Teuchos::Array<std::string> region_list = make_regions_list(regStr);
		    if (!compare_region_names(region_list, def_list)) {
		      Errors::Message msg;
		      msg << "Amanzi::InputTranslator: ERROR - invalid region in write_regions Section - " ;
		      msg << "valid regions are: \n" ;
		      for (int r=0; r<regionNames_string_.size(); r++) {
			msg << "    " << regionNames_string_[r] << "\n";
		      }
		      msg << "specified regions: \n" ;
		      for (int r=0; r<region_list.size(); r++) {
			msg << "    " << region_list[r] << "\n";
		      }
		      msg << "  Please correct and try again \n" ;
		      Exceptions::amanzi_throw(msg);
		    }
		    wrPL.set<Teuchos::Array<std::string> >(nameStr,region_list);
		    XMLString::release(&nameStr);
		    XMLString::release(&regStr);
		  }
		}
		visPL.sublist("Write Regions") = wrPL;
	      }
	    }
            XMLString::release(&textContent);
          }
          list.sublist("Visualization Data") = visPL;

	}
        else if (strcmp(outName,"checkpoint")==0) {
          if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
            *voI_->os() << "Getting Output: Checkpoint" << std::endl;
          }
          // get list of checkpoint - this node MAY exist ONCE
          Teuchos::ParameterList chkPL;
          DOMNodeList* childList = curoutNode->getChildNodes();
          for (int j=0; j<childList->getLength(); j++) {
            DOMNode* curKid = childList->item(j) ;
            textContent  = XMLString::transcode(curKid->getNodeName());
            if (strcmp(textContent,"base_filename")==0) {
	      textContent2 = XMLString::transcode(curKid->getTextContent());
              chkPL.set<std::string>("File Name Base",trim_string(textContent2));
              XMLString::release(&textContent2);
	    }
            else if (strcmp(textContent,"num_digits")==0) {
	      textContent2 = XMLString::transcode(curKid->getTextContent());
              chkPL.set<int>("File Name Digits",get_int_constant(textContent2,def_list));
              XMLString::release(&textContent2);
	    }
            else if (strcmp(textContent,"cycle_macros")==0 || strcmp(textContent,"cycle_macro")==0) {
              textContent2 = XMLString::transcode(curKid->getTextContent());
              Teuchos::Array<std::string> macro = make_regions_list(textContent2);
              chkPL.set<Teuchos::Array<std::string> >("Cycle Macros",macro);
              XMLString::release(&textContent2);
	    }
            XMLString::release(&textContent);
          }
          list.sublist("Checkpoint Data") = chkPL;
	}
        else if (strcmp(outName,"walkabout")==0) {
          if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
            *voI_->os() << "Getting Output: Walkabout" << std::endl;
          }
          // get list of walkabout - this node MAY exist ONCE
          Teuchos::ParameterList chkPL;
          DOMNodeList* childList = curoutNode->getChildNodes();
          for (int j=0; j<childList->getLength(); j++) {
            DOMNode* curKid = childList->item(j) ;
            textContent  = XMLString::transcode(curKid->getNodeName());
            if (strcmp(textContent,"base_filename")==0) {
	      textContent2 = XMLString::transcode(curKid->getTextContent());
              chkPL.set<std::string>("File Name Base",trim_string(textContent2));
              XMLString::release(&textContent2);
	    }
            else if (strcmp(textContent,"num_digits")==0) {
	      textContent2 = XMLString::transcode(curKid->getTextContent());
              chkPL.set<int>("File Name Digits",get_int_constant(textContent2,def_list));
              XMLString::release(&textContent2);
	    }
            else if (strcmp(textContent,"cycle_macros")==0 || strcmp(textContent,"cycle_macro")==0) {
              textContent2 = XMLString::transcode(curKid->getTextContent());
              Teuchos::Array<std::string> macro = make_regions_list(textContent2);
              chkPL.set<Teuchos::Array<std::string> >("Cycle Macros",macro);
              XMLString::release(&textContent2);
	    }
            XMLString::release(&textContent);
          }
          list.sublist("Walkabout Data") = chkPL;
	}
        else if (strcmp(outName,"observations")==0) {
          if (voI_->getVerbLevel() >= Teuchos::VERB_HIGH) {
            *voI_->os() << "Getting Output: Obs" << std::endl;
          }
          Teuchos::ParameterList obsPL;
          DOMNodeList* OBList = curoutNode->getChildNodes();
          for (int i=0; i<OBList->getLength(); i++) {
            DOMNode* curNode = OBList->item(i) ;
            if (DOMNode::ELEMENT_NODE == curNode->getNodeType()) {
              textContent  = XMLString::transcode(curNode->getNodeName());
              if (strcmp(textContent,"filename")==0) {
	        textContent2 = XMLString::transcode(curNode->getTextContent());
                obsPL.set<std::string>("Observation Output Filename",trim_string(textContent2));
	        XMLString::release(&textContent2);
              }
              else if (strcmp(textContent,"liquid_phase")==0) {
                DOMNamedNodeMap* attrMap = curNode->getAttributes();
                DOMNode* nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
                std::string phaseName;
		if (nodeAttr) {
                  textContent2 = XMLString::transcode(nodeAttr->getNodeValue());
                  phaseName = std::string(textContent2);
                  if (phaseName=="water") {
                    phaseName="Water";
                  }
                  XMLString::release(&textContent2);
	        }
                else {
                  throw_error_missattr("observations", "attribute", "name", "liquid_phase");
	        }

	        // loop over observations
	        DOMNodeList* childList = curNode->getChildNodes();
                for (int j=0; j<childList->getLength(); j++) {
	          Teuchos::ParameterList obPL;
                  DOMNode* curObs = childList->item(j) ;
                  if (DOMNode::ELEMENT_NODE == curObs->getNodeType()) {
                    char* obsType = XMLString::transcode(curObs->getNodeName());
                    if (strcmp(obsType,"aqueous_pressure")==0) {
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Aqueous pressure");
                      }
                      else {
                        obPL.set<std::string>("Variable","Aqueous_Pressure");
                      }
	            }
                    else if (strcmp(obsType,"integrated_mass")==0) {
	              // TODO: EIB can't find matching version
	              //obPL.set<std::string>("Variable","Aqueous pressure");
                    }
                    else if (strcmp(obsType,"volumetric_water_content")==0) {
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Volumetric water content");
                      }
                      else {
                        obPL.set<std::string>("Variable","Volumetric_Water_Content");
                      }
                    }
                    else if (strcmp(obsType,"gravimetric_water_content")==0) {
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Gravimetric water content");
                      }
                      else {
                        //TODO: EIB - don't think this is in structured
                      }
                    }
                    else if (strcmp(obsType,"x_aqueous_volumetric_flux")==0) {
                      // TODO: EIB needs double checking
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","X-Aqueous volumetric flux");
                      }
                      else {
                        obPL.set<std::string>("Variable","Aqueous_Volumetric_Flux_X");
                      }
                    }
                    else if (strcmp(obsType,"y_aqueous_volumetric_flux")==0) {
                      // TODO: EIB needs double checking
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Y-Aqueous volumetric flux");
                      }
                      else {
                        obPL.set<std::string>("Variable","Aqueous_Volumetric_Flux_Y");
                      }
                    }
                    else if (strcmp(obsType,"z_aqueous_volumetric_flux")==0) {
                      // TODO: EIB needs double checking
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Z-Aqueous volumetric flux");
                      }
                      else {
                        obPL.set<std::string>("Variable","Aqueous_Volumetric_Flux_Z");
                      }
                    }
                    else if (strcmp(obsType,"material_id")==0) {
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","MaterialID");
                      }
                      else {
                        obPL.set<std::string>("Variable","Material_ID");
                      }
                    }
                    else if (strcmp(obsType,"hydraulic_head")==0) {
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Hydraulic Head");
                      }
                      else {
                        obPL.set<std::string>("Variable","Hydraulic_Head");
                      }
                    }
                    else if (strcmp(obsType,"aqueous_mass_flow_rate")==0) {
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Aqueous mass flow rate");
                      }
                      else {
                        //TODO: EIB - don't think this is in structured
                      }
                    }
                    else if (strcmp(obsType,"aqueous_volumetric_flow_rate")==0) {
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Aqueous volumetric flow rate");
                      }
                      else {
                        //TODO: EIB - don't think this is in structured
                      }
		    }
                    else if (strcmp(obsType,"aqueous_saturation")==0) {
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Aqueous saturation");
                      }
                      else {
                        obPL.set<std::string>("Variable","Aqueous Saturation");
                      }
                    }
                    else if (strcmp(obsType,"aqueous_conc")==0) {
	              // get solute name
                      DOMNamedNodeMap* attrMap = curObs->getAttributes();
                      DOMNode* nodeAttr = attrMap->getNamedItem(XMLString::transcode("solute"));
	              char* soluteName ;
		      if (nodeAttr) {
                        soluteName = XMLString::transcode(nodeAttr->getNodeValue());
	              }
                      else {
                        throw_error_missattr("observations", "attribute", "solute", "aqueous_conc");
	              }

	              std::stringstream name;
                      if (isUnstr_) {
                        name<< soluteName << " Aqueous concentration";
                      }
                      else {
                        name<< soluteName << "_Aqueous_Concentration";
                      }
	              obPL.set<std::string>("Variable",name.str());
                    }
                    else if (strcmp(obsType,"drawdown")==0) {
                      if (isUnstr_) {
                        obPL.set<std::string>("Variable","Drawdown");
                      }
                      else {
                        //TODO: EIB - don't think this is in structured
                      }
	            }
                    else if (strcmp(obsType,"solute_volumetric_flow_rate")==0) {
                      // get solute name
                      DOMNamedNodeMap* attrMap = curObs->getAttributes();
                      DOMNode* nodeAttr = attrMap->getNamedItem(XMLString::transcode("solute"));
                      char* soluteName ;
                      if (nodeAttr) {
                        soluteName = XMLString::transcode(nodeAttr->getNodeValue());
                      }
                      else {
                        throw_error_missattr("observations", "attribute", "solute", "solute_volumetric_flow_rate");
                      }
                      
                      std::stringstream name;
                      if (isUnstr_) {
                        name<< soluteName << " volumetric flow rate";
                      }
                      else {
                        // TODO: EIB not sure this is in structured yet
                        name<< soluteName << "_volumetric_flow_rate";
                      }
                      obPL.set<std::string>("Variable",name.str());
                    }
	            DOMNodeList* kidList = curObs->getChildNodes();
                    for (int k=0; k<kidList->getLength(); k++) {
                      DOMNode* curElem = kidList->item(k) ;
                      if (DOMNode::ELEMENT_NODE == curElem->getNodeType()) {
                        char* Elem =  XMLString::transcode(curElem->getNodeName());
                        char* Value =  XMLString::transcode(curElem->getTextContent());
		        if (strcmp(Elem,"assigned_regions")==0) {
		          //TODO: EIB - really a note, REGION != ASSIGNED REGIONS, this isn't consistent!!!
		          /*
	                  Teuchos::Array<std::string> regs;
	                  char* char_array;
	                  char_array = strtok(Value,",");
	                  while(char_array!=NULL){
	                    regs.append(char_array);
	                    char_array = strtok(NULL,",");
	                  }
                          obPL.set<Teuchos::Array<std::string> >("Region",regs);
		          */
                          obPL.set<std::string>("Region",trim_string(Value));
		        }
                        else if (strcmp(Elem,"functional")==0) {
	                  if (strcmp(Value,"point")==0) {
	                    obPL.set<std::string>("Functional","Observation Data: Point");
	                  }
                          else if (strcmp(Value,"integral")==0) {
	                    obPL.set<std::string>("Functional","Observation Data: Integral");
	                  }
                          else if (strcmp(Value,"mean")==0) {
			    obPL.set<std::string>("Functional","Observation Data: Mean");
	                  }
		        }
                        // Keeping singular macro around to help users.  This will go away
                        else if (strcmp(Elem,"time_macros")==0 || strcmp(Elem,"time_macro")==0) {
                          Teuchos::Array<std::string> macros = make_regions_list(Value);
	                  obPL.set<Teuchos::Array<std::string> >("Time Macros",macros);
                        }
                        else if (strcmp(Elem,"cycle_macros")==0 || strcmp(Elem,"cycle_macro")==0) {
                          Teuchos::Array<std::string> macros = make_regions_list(Value);
                          obPL.set<Teuchos::Array<std::string> >("Cycle Macros",macros);
                        }
                        XMLString::release(&Elem);
                        XMLString::release(&Value);
	              }
	            }
	            std::stringstream listName;
	            listName << "observation-"<<j+1<<":"<<phaseName;
	            obsPL.sublist(listName.str()) = obPL;
	          }
	        }
              }
              XMLString::release(&textContent);
              list.sublist("Observation Data") = obsPL;
            }
          }
	}
      }
    }
  }

  return list;
  
}

/*
 ******************************************************************
 * Empty
 ******************************************************************
 */

  //MSD - get gslib information.
  void get_gslib_info(Teuchos::ParameterList& propertyList,const Teuchos::ParameterList& defList, const xercesc::DOMElement* propElement, const std::string& propName, const std::string& sectionName, const std::string& dfile_DEF)
  {
    const std::string pfile_label("parameter_file");
    const std::string dfile_label("data_file");
    const std::string value_label("value");
    const std::string gslib_pfile_label("GSLib Parameter File");
    const std::string gslib_dfile_label("GSLib Data File");
    const std::string gslib_value_label("Value");

    std::string pfile, dfile = dfile_DEF;;
    double value;

    if (propElement->hasAttribute(XMLString::transcode(pfile_label.c_str()))) {
      pfile = std::string(XMLString::transcode(propElement->getAttribute(XMLString::transcode(pfile_label.c_str()))));
    }
    else {
      throw_error_missattr(sectionName, "attribute", pfile_label, propName);
    }

    if (propElement->hasAttribute(XMLString::transcode(value_label.c_str()))) {
      std::string vstr = std::string(XMLString::transcode(propElement->getAttribute(XMLString::transcode(value_label.c_str()))));
      value = get_double_constant(vstr,defList);
    }
    else {
      throw_error_missattr(sectionName, "attribute", value_label, propName);
    }

    if (propElement->hasAttribute(XMLString::transcode(dfile_label.c_str()))) {
      dfile = std::string(XMLString::transcode(propElement->getAttribute(XMLString::transcode(dfile_label.c_str()))));
    }

    propertyList.set<std::string>(gslib_pfile_label,pfile);
    propertyList.set<std::string>(gslib_dfile_label,dfile);
    propertyList.set<double>(gslib_value_label,value);
  }

/*
 ******************************************************************
 * Empty
 ******************************************************************
 */
  
  //TODO: EIB - get file information.
  
  Teuchos::ParameterList get_file_info(Teuchos::ParameterList propertyList, DOMElement* propElement, std::string propName, std::string sectionName)
  {
    char* textContent;
    std::string text;
    bool isExodus = false;
    
    if (propElement->hasAttribute(XMLString::transcode("type"))) {
      text = std::string(XMLString::transcode(propElement->getAttribute(XMLString::transcode("type"))));
      // TODO: move all to lower case
      if (text == "file" || text == "exodus ii") {
        propertyList.set<std::string>("Format","Exodus II");
        isExodus = true;
      }
      else if (text == "color") {
        propertyList.set<std::string>("Format","Color Function");
      }
      else {
        throw_error_illformed(sectionName, "attribute", "type", "'file', 'exodus ii' or 'color'");
      }
    }
    else {
      throw_error_missattr(sectionName, "attribute", "type", propName);
    }
    
    if (propElement->hasAttribute(XMLString::transcode("filename"))) {
      textContent = XMLString::transcode(propElement->getAttribute(XMLString::transcode("filename")));
      propertyList.set<std::string>("File",trim_string(textContent));
      XMLString::release(&textContent);
    }
    else {
      throw_error_missattr(sectionName, "attribute", "filename", propName);
    }
    
    if (isExodus) {
      if (propElement->hasAttribute(XMLString::transcode("attribute"))) {
        textContent = XMLString::transcode(propElement->getAttribute(XMLString::transcode("attribute")));
        propertyList.set<std::string>("Attribute",trim_string(textContent));
        XMLString::release(&textContent);
      }
      else {
        throw_error_missattr(sectionName, "attribute", "attribute", propName);
      }
    }
    
    return propertyList;
  }

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */

//TODO: EIB - get default time unit from units, convert plain time values if not seconds.

double get_time_value(std::string time_value, Teuchos::ParameterList def_list)
{

  double time;

  // Check if time_value is: listed in constants
  if (def_list.sublist("constants").sublist("constants").isSublist(time_value)) {
    time = def_list.sublist("constants").sublist("constants").sublist(time_value).get<double>("value");

  // Check if time_value is: listed in time_constants
  }
  else if (def_list.sublist("constants").sublist("time_constants").isSublist(time_value)) {
    time = def_list.sublist("constants").sublist("time_constants").sublist(time_value).get<double>("value");

  // Otherwise, it must be a string or typedef_labeled_time (0000.00,y)
  }
  else {
    char* tmp = strcpy(new char[time_value.size() + 1], time_value.c_str());
    time = convert_time_value(tmp);
    delete[] tmp;
  }

  return time;
}

/*
******************************************************************
* Empty
******************************************************************
*/
  
//TODO: EIB - get default time unit from units, convert plain time values if not seconds.
  
double convert_time_value(char* time_value)
{
  double time;
  char*  char_array;
  
  char_array = strtok(time_value,";, ");
  time = atof(char_array);
  char_array = strtok(NULL,";,");
  if (char_array!=NULL) {
    if (strcmp(char_array,"y")==0) { time = time*365.25*24.0*60.0*60.0; }
    else if (strcmp(char_array,"d")==0) { time = time*24.0*60.0*60.0; }
    else if (strcmp(char_array,"h")==0) { time = time*60.0*60.0; }
  }
  
  return time;
}
  
/*
 ******************************************************************
 * Empty
 ******************************************************************
 */
double get_double_constant(std::string pos_name, Teuchos::ParameterList def_list)
{

  double value;

  // Check if pos_name is: listed in constants
  if (def_list.sublist("constants").sublist("constants").isSublist(pos_name)) {
    value = def_list.sublist("constants").sublist("constants").sublist(pos_name).get<double>("value");

  // Check if pos_name is: listed in time_constants
  }
  else if (def_list.sublist("constants").sublist("numerical_constant").isSublist(pos_name)) {
    value = def_list.sublist("constants").sublist("numerical_constant").sublist(pos_name).get<double>("value");

  // Otherwise, we must assume it's already a value
  }
  else {
    value = atof(pos_name.c_str());
  }
  
  return value;
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
int get_int_constant(std::string pos_name, Teuchos::ParameterList def_list)
{

  int value;

  // Check if pos_name is: listed in constants
  if (def_list.sublist("constants").sublist("constants").isSublist(pos_name)) {
    value = def_list.sublist("constants").sublist("constants").sublist(pos_name).get<int>("value");

  // Check if pos_name is: listed in time_constants
  }
  else if (def_list.sublist("constants").sublist("numerical_constant").isSublist(pos_name)) {
    value = def_list.sublist("constants").sublist("numerical_constant").sublist(pos_name).get<int>("value");

  // Otherwise, we must assume it's already a value
  }
  else {
    value = atoi(pos_name.c_str());
  }
  
  return value;
}

/*
******************************************************************
* Empty
******************************************************************
*/
  Teuchos::Array<int> make_int_list(char* char_array)
  {
    Teuchos::Array<int> int_list;
    char* tmp;
    tmp = strtok(char_array,", ");
    while(tmp!=NULL){
      std::string str(tmp);
      boost::algorithm::trim(str);
      int_list.append(atoi(str.c_str()));
      tmp = strtok(NULL,", ");
    }
    
    return int_list;
  }
/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::Array<std::string> make_regions_list(char* char_array)
{
  Teuchos::Array<std::string> regs;
  char* tmp;
  tmp = strtok(char_array,",");
  while(tmp!=NULL){
    std::string str(tmp);
    boost::algorithm::trim(str);
    regs.append(str);
    tmp = strtok(NULL,",");
  }

  return regs;
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
std::string trim_string(char* tmp)
{
    std::string str(tmp);
    boost::algorithm::trim(str);
    return str;
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
bool compare_region_names(Teuchos::Array<std::string> regions, Teuchos::ParameterList def_list)
{
  int cnt;
  cnt = 0;
  bool status=false;
  for (int i = 0; i < regions.size(); i++) {
    for (int r = 0; r < regionNames_string_.size(); r++) {
      if (strcmp(regionNames_string_[r].c_str(),regions[i].c_str())==0) {
        status=true;
      }
    }
    if (!status) {
      std::cout << "Amanzi::InputTranslator: ERROR - region "<< regions[i] << " NOT in known regions!" << std::endl;
      return status;
    }
  }

  return status;
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::Array<double> make_coordinates(char* char_array, Teuchos::ParameterList def_list)
{
  Teuchos::Array<double> coords;
  char* tmp;
  tmp = strtok(char_array,"(, ");
  while(tmp!=NULL){
    std::string str(tmp);
    boost::algorithm::trim(str);
    coords.append(get_double_constant(str, def_list));
    tmp = strtok(NULL,",");
  }

  return coords;
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList make_chemistry(Teuchos::ParameterList def_list)
{
  Teuchos::ParameterList chemistryPL;
  Teuchos::ParameterList bgdPL;

  // Get common options
  Teuchos::Array<std::string> verb;
  if (def_list.sublist("simulation").isParameter("verbosity")) {
    if (voI_->getVerbLevel() == Teuchos::VERB_EXTREME) {
      verb.append("error");
      chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
    }
    else if (voI_->getVerbLevel() == Teuchos::VERB_HIGH) {
      verb.append("warning");
      chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
    }
    else if (voI_->getVerbLevel() == Teuchos::VERB_MEDIUM) {
      verb.append("verbose");
      chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
    }
    else if (voI_->getVerbLevel() == Teuchos::VERB_LOW) {
      verb.append("terse");
      chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
    }
    else {
      verb.append("silent");
      chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
    }
  }
  else {
    verb.append("silent");
    chemistryPL.set<Teuchos::Array<std::string> >("Verbosity",verb);
  }
  
  
  if (def_list.isParameter("chemistry_engine")) {
    // get Amanzi Native options
    if (def_list.get<std::string>("chemistry_engine") == "amanzi") {
      // go ahead and add bdg file to PL
      // build bgd filename
      std::string bgdfilename;
      if (def_list.isParameter("xmlfilename") ) {
        bgdfilename = def_list.get<std::string>("xmlfilename") ;
        std::string new_extension(".bgd");
        size_t pos = bgdfilename.find(".xml");
        bgdfilename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)4);
      }
      else {
        // defaulting to hardcoded name
        bgdfilename = "isotherms.bgd" ;
      }
      // add bgd file and parameters to list
      Teuchos::ParameterList bgdPL;
      bgdPL.set<std::string>("Format","simple");
      bgdPL.set<std::string>("File",bgdfilename);
      chemistryPL.sublist("Thermodynamic Database") = bgdPL;
      chemistryPL.set<std::string>("Activity Model","unit");
    }
    // get Alquimia options
    else if (def_list.get<std::string>("chemistry_engine") == "pflotran") {
      chemistryPL.set<std::string>("Engine","PFloTran");
      
      // check file *.in filename in def_list
      if (def_list.isSublist("chemistry_PL")) {
        if (def_list.sublist("chemistry_PL").isParameter("Engine Input File")) {
          chemistryPL.set<std::string>("Engine Input File",def_list.sublist("chemistry_PL").get<std::string>("Engine Input File"));
        }
      }
    }
  }

    // fill in default values
    //chemistryPL.set<double>("Tolerance",1e-12);
    //chemistryPL.set<int>("Maximum Newton Iterations",200);
    //chemistryPL.set<double>("Max Time Step (s)",9e9);
  
    return chemistryPL;
}


/*
 ******************************************************************
 * Empty
 ******************************************************************
 */
void write_BDG_file(Teuchos::ParameterList sorption_list, Teuchos::ParameterList def_list)
{

  std::ofstream bgd_file;
  std::stringstream species_string;
  std::stringstream isotherms_string;

  // build streams
  for (Teuchos::ParameterList::ConstIterator i = sorption_list.begin(); i != sorption_list.end(); i++) {
      Teuchos::ParameterList& tmpList = sorption_list.sublist(sorption_list.name(i)) ;
      species_string << sorption_list.name(i) << " ;   0.00 ;   0.00 ;   1.00 \n";
      if ( tmpList.isParameter("Langmuir b") ) {
          isotherms_string << sorption_list.name(i) << " ; langmuir ; " << tmpList.get<double>("Kd")<< " " <<tmpList.get<double>("Langmuir b") << std::endl;
      }
      else if ( tmpList.isParameter("Freundlich n") ) {
          isotherms_string << sorption_list.name(i) << " ; freundlich ; " << tmpList.get<double>("Kd")<< " " <<tmpList.get<double>("Freundlich n") << std::endl;
      }
      else {
          isotherms_string << sorption_list.name(i) << " ; linear ; " << tmpList.get<double>("Kd")<< std::endl;
      }
  }
  
  // build bgd filename
  std::string bgdfilename;
  if (def_list.isParameter("xmlfilename") ) {
      bgdfilename = def_list.get<std::string>("xmlfilename") ;
      std::string new_extension(".bgd");
      size_t pos = bgdfilename.find(".xml");
      bgdfilename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)4);
  }
  else {
      // defaulting to hardcoded name
      bgdfilename = "isotherms.bgd" ;
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

/*
**********************************************************************
* Generate unified warning message for unrecognized element (skipping)
**********************************************************************
*/
void throw_warning_skip(std::string element){
  
  *voI_->os() << "Amanzi::InputTranslator: WARNING - The following is an unrecognized element and is being skipped during translation." << std::endl << element << " was unrecognized option." << std::endl << " Please check input file schema for correct naming and syntax." << std::endl;
}
  
  
/*
 *******************************************************************
 * Generate unified error message for Structure/Unstructure conflict
 *******************************************************************
*/
  void throw_error_str_ustr(std::string section, std::string element_type, std::string sim_type){
    
    Errors::Message msg;
    msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing " << section << "- " ;
    msg << section << " type " << element_type << " is only available for " << sim_type << ". \n" ;
    msg << "  Please correct and try again \n" ;
    Exceptions::amanzi_throw(msg);
  }
  
/*
 *******************************************************************
 * Generate unified error message for ill-formed element
 *******************************************************************
*/
  void throw_error_illformed(std::string section, std::string element_type, std::string ill_formed){
    
    Errors::Message msg;
    msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing " << section << "- " ;
    msg << "  Missing or Ill-formed '" << element_type << "' for '" << ill_formed << "'. \n" ;
    msg << "  Please correct and try again \n" ;
    Exceptions::amanzi_throw(msg);
  }
  
/*
 *****************************************************************************
 * Generate unified error message for ill-formed element with options provided
 *****************************************************************************
*/
  void throw_error_illformed(std::string section, std::string element_type, std::string ill_formed, std::string options){
    
    Errors::Message msg;
    msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing " << section << " - " ;
    msg << "  Missing or Ill-formed '" << element_type << "' for '" << ill_formed << "'. Valid options are: " << options << "\n" ;
    msg << "  Please correct and try again \n" ;
    Exceptions::amanzi_throw(msg);
  }
  
/*
 *******************************************************************
 * Generate unified error message for missing item
 *******************************************************************
*/
  void throw_error_missattr(std::string section, std::string att_elem_type, std::string missing, std::string elem_name){
    
    Errors::Message msg;
    msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing " << section << " - \n" ;
    msg << "  No " << att_elem_type << " " << missing << " found for " << elem_name << ". \n" ;
    msg << "  Please correct and try again \n" ;
    Exceptions::amanzi_throw(msg);
  }

/*
 ***********************************************************************
 * Generate unified error message for missing item with options provided
 ***********************************************************************
*/
  void throw_error_missattr(std::string section, std::string att_elem_type, std::string missing, std::string elem_name, std::string options){
    
    Errors::Message msg;
    msg << "Amanzi::InputTranslator: ERROR - An error occurred during parsing " << section << " - " ;
    msg << "  No " << att_elem_type << " " << missing << " found for " << elem_name << ". Options are: " << options <<"\n" ;
    msg << "  Please correct and try again \n" ;
    Exceptions::amanzi_throw(msg);
  }
  
} // end namespace AmanziNewInput
} // end namespace Amanzi
