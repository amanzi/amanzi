#include "InputTranslator.hh"
#include "InputParserIS-defaults.hh"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <sstream>
#include <string>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>

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
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include "DOMTreeErrorReporter.hpp"
#include "ErrorHandler.hpp"

namespace Amanzi {
namespace AmanziNewInput {

/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList translate(const std::string& xmlfilename, const std::string& xmlSchemafile) {

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
  xercesc::XMLPlatformUtils::Initialize();
  xercesc::XercesDOMParser *parser = new xercesc::XercesDOMParser;
  if (parser->loadGrammar(xmlSchemafile.c_str(), 
			  xercesc::Grammar::SchemaGrammarType, true) == NULL) {
      amanzi_throw(Errors::Message("ERROR: didn't load grammar"));
  }
  parser->useCachedGrammarInParse( true );
  AmanziErrorHandler* errHandler = new AmanziErrorHandler();
  parser->setErrorHandler(errHandler);
  parser->setValidationScheme(xercesc::XercesDOMParser::Val_Always);
  parser->setDoNamespaces(true);
  parser->setDoSchema(true);
  bool errorsOccured = false;
  try{
      parser->parse(xmlfilename.c_str());
  }
  catch (const xercesc::OutOfMemoryException&)
  {
      std::cerr << "OutOfMemoryException" << std::endl;
      errorsOccured = true;
  }

  // check that it's validating here

  // grab the version number attribute
  new_list.set<std::string>("Amanzi Input Format Version", "1.1.0");

  // grab the mesh type
  //new_list.sublist(framework) = ...;

  // go through each section, if it exist in the file, translate it 
  // to the old format
  xercesc::DOMDocument *doc = parser->getDocument();
  
  def_list.sublist("constants") = get_constants(doc);

  new_list.sublist("General Description") = get_model_description(doc);
  new_list.sublist("Mesh") = get_Mesh(doc);
  new_list.sublist("Domain").set<int>("Spatial Dimension",dimension_);
  new_list.sublist("Execution Control") = get_execution_controls(doc,def_list);
  new_list.sublist("Phase Definitions") = get_phases(doc);
  new_list.sublist("Regions") = get_regions(doc);
  new_list.sublist("Material Properties") = get_materials(doc);
  new_list.sublist("Initial Conditions") = get_initial_conditions(doc);
  new_list.sublist("Boundary Conditions") = get_boundary_conditions(doc,def_list);
  new_list.sublist("Sources") = get_sources(doc,def_list);
  new_list.sublist("Output") = get_output(doc);
  
  
  xercesc::XMLPlatformUtils::Terminate();

  // return the completely translated input file as a parameter list
  return new_list;
}

/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList get_constants(xercesc::DOMDocument* xmlDoc) {

  Teuchos::ParameterList list;

  xercesc::DOMNamedNodeMap *attrMap;
  xercesc::DOMNode *namedNode;
  char* name;
  char* type;
  char* value;
  char* char_array;
  double time;

  // read in new stuff
  xercesc::DOMNodeList* nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("definitions"));

  if (nodeList->getLength() > 0) {
    xercesc::DOMNode* nodeD = nodeList->item(0);
    xercesc::DOMNodeList* childern = nodeD->getChildNodes();
    for (int i=0; i<childern->getLength(); i++) {
      xercesc::DOMNode* currentNode = childern->item(i) ;
      if (xercesc::DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
        char* tagname = xercesc::XMLString::transcode(currentNode->getNodeName());
	// deal with: constants, named_times, macros
        if (strcmp(tagname,"constants")==0) {
          xercesc::DOMNodeList* kids = currentNode->getChildNodes();
          for (int j=0; j<kids->getLength(); j++) {
            xercesc::DOMNode* currentKid = kids->item(j) ;
            if (xercesc::DOMNode::ELEMENT_NODE == currentKid->getNodeType()) {
              char* kidname = xercesc::XMLString::transcode(currentKid->getNodeName());
              // types: constant, time_constant, numerical_constant, area_mass_flux_constant
              if (strcmp(kidname,"constant")==0) {
	        attrMap = currentKid->getAttributes();
	        namedNode = attrMap->getNamedItem(XMLString::transcode("name"));
	        name = xercesc::XMLString::transcode(namedNode->getNodeValue());
	        namedNode = attrMap->getNamedItem(XMLString::transcode("type"));
	        type = xercesc::XMLString::transcode(namedNode->getNodeValue());
	        namedNode = attrMap->getNamedItem(XMLString::transcode("value"));
	        value = xercesc::XMLString::transcode(namedNode->getNodeValue());
		if (strcmp(type,"time")==0) {
		  // check if time and convert to seconds - year = 365.25
		  // TODO: EIB - verify this works with spaces
		  // TODO: EIB - expect Akuna to move to no deliminator, need to test for this
		  char_array = strtok(value,";, ");
		  time = atof(char_array);
		  char_array = strtok(NULL,";,");
		  if (strcmp(char_array,"y")==0) { time = time*31577600.0; }
		  else if (strcmp(char_array,"d")==0) { time = time*86100.0; }
		  else if (strcmp(char_array,"h")==0) { time = time*3600.0; }
		}
		// add to list
		Teuchos::ParameterList tmp;
		tmp.set<std::string>("type",type);
		tmp.set<double>("value",time);
		list.sublist("constants").sublist(name) = tmp;
                XMLString::release(&name);
                XMLString::release(&type);
                XMLString::release(&value);
	      } else if (strcmp(kidname,"time_constant")==0) {
	        attrMap = currentKid->getAttributes();
	        namedNode = attrMap->getNamedItem(XMLString::transcode("name"));
	        name = xercesc::XMLString::transcode(namedNode->getNodeValue());
	        namedNode = attrMap->getNamedItem(XMLString::transcode("value"));
	        value = xercesc::XMLString::transcode(namedNode->getNodeValue());
		// check if time and convert to seconds - year = 365.25
		// TODO: EIB - verify this works with spaces
		// TODO: EIB - expect Akuna to move to no deliminator, need to test for this
		char_array = strtok(value,";, ");
		time = atof(char_array);
		char_array = strtok(NULL,";,");
		if (strcmp(char_array,"y")==0) { time = time*31577600.0; }
		else if (strcmp(char_array,"d")==0) { time = time*86100.0; }
		else if (strcmp(char_array,"h")==0) { time = time*3600.0; }
		// add to list
		Teuchos::ParameterList tmp;
		tmp.set<double>("value",time);
		list.sublist("time_constants").sublist(name) = tmp;
                XMLString::release(&name);
                XMLString::release(&value);
	      } else if (strcmp(kidname,"numerical_constant")==0) {
	        attrMap = currentKid->getAttributes();
	        namedNode = attrMap->getNamedItem(XMLString::transcode("name"));
	        name = xercesc::XMLString::transcode(namedNode->getNodeValue());
	        namedNode = attrMap->getNamedItem(XMLString::transcode("value"));
	        value = xercesc::XMLString::transcode(namedNode->getNodeValue());
		// add to list
		Teuchos::ParameterList tmp;
		tmp.set<double>("value",atof(value));
		list.sublist("numerical_constant").sublist(name) = tmp;
                XMLString::release(&name);
                XMLString::release(&value);
	      } else if (strcmp(kidname,"area_mass_flux_constant")==0) {
	        attrMap = currentKid->getAttributes();
	        namedNode = attrMap->getNamedItem(XMLString::transcode("name"));
	        name = xercesc::XMLString::transcode(namedNode->getNodeValue());
	        namedNode = attrMap->getNamedItem(XMLString::transcode("value"));
	        value = xercesc::XMLString::transcode(namedNode->getNodeValue());
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
	} else if (strcmp(tagname,"named_times")==0) {
	  //TODO: EIB - deal with named times
          xercesc::DOMNodeList* kids = currentNode->getChildNodes();
          for (int j=0; j<kids->getLength(); j++) {
            xercesc::DOMNode* currentKid = kids->item(j) ;
            if (xercesc::DOMNode::ELEMENT_NODE == currentKid->getNodeType()) {
              char* kidname = xercesc::XMLString::transcode(currentKid->getNodeName());
              // types: time
              if (strcmp(kidname,"time")==0) {
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

/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList get_model_description(xercesc::DOMDocument* xmlDoc) {

  Teuchos::ParameterList list;

  // read in new stuff
  XMLCh* tag = XMLString::transcode("model_description");
  xercesc::DOMNodeList* nodeList = xmlDoc->getElementsByTagName(tag);
  XMLString::release(&tag);

  // write to old format, mostly won't be read so just write out as strings
  const XMLSize_t nodeCount = nodeList->getLength() ;
  if (nodeList->getLength() > 0) {
    xercesc::DOMNode* nodeGD = nodeList->item(0);
    xercesc::DOMElement* elementGD = static_cast<xercesc::DOMElement*>(nodeGD);
    char* model_name = xercesc::XMLString::transcode(elementGD->getAttribute(
		       xercesc::XMLString::transcode("name")));
    list.set<std::string>("model_name",model_name);

    xercesc::DOMNodeList* childern = nodeGD->getChildNodes();
    for (int i=0; i<childern->getLength(); i++) {
      xercesc::DOMNode* currentNode = childern->item(i) ;
      if (xercesc::DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
        char* tagname = xercesc::XMLString::transcode(currentNode->getNodeName());
        xercesc::DOMNode::NodeType type = currentNode->getNodeType();
        if (strcmp(tagname,"units")!=0) {
          char* textContent = XMLString::transcode(currentNode->getTextContent());
          list.set<std::string>(tagname,textContent);
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

/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList get_Mesh(xercesc::DOMDocument* xmlDoc ) {

  //TODO: EIB - see if mesh list ending up in right order/structure

  Teuchos::ParameterList list;

  bool unstructured = true;
  bool generate = true;
  bool file = false;
  char *framework;
  Teuchos::ParameterList mesh_list;

  // read in new stuff
  XMLCh* tag = XMLString::transcode("mesh");
  xercesc::DOMNodeList* nodeList = xmlDoc->getElementsByTagName(tag);
  XMLString::release(&tag);

  // read the attribute to set the framework sublist
  const XMLSize_t nodeCount = nodeList->getLength() ;
  if (nodeList->getLength() > 0) {
    xercesc::DOMNode* nodeMesh = nodeList->item(0);
    xercesc::DOMElement* elementMesh = static_cast<xercesc::DOMElement*>(nodeMesh);
    if(nodeMesh->hasAttributes()) {
      // unstructured
      framework = xercesc::XMLString::transcode(elementMesh->getAttribute(
		  xercesc::XMLString::transcode("framework")));
      unstructured = true;
    }
    else { 
      // structured
      unstructured = false;
    }

    // loop over child nodes
    xercesc::DOMNodeList* children = nodeMesh->getChildNodes();
    for (int i=0; i<children->getLength(); i++) {
      xercesc::DOMNode* currentNode = children->item(i) ;
      if (xercesc::DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
	char* tagname = xercesc::XMLString::transcode(currentNode->getNodeName());
	if (strcmp(tagname,"dimension")==0) {
	  char* temp = xercesc::XMLString::transcode(currentNode->getTextContent());
	  dimension_ = atoi(temp);
	} else if (strcmp(tagname,"generate")==0) {
	  generate = true;
	  file = false;
	  xercesc::DOMElement* elementGen = static_cast<xercesc::DOMElement*>(currentNode);

	  // get Number of Cells
	  Teuchos::Array<int> ncells; 
	  xercesc::DOMNodeList* nodeList = elementGen->getElementsByTagName(
			                   XMLString::transcode("number_of_cells"));

	  xercesc::DOMNode* node = nodeList->item(0);
	  xercesc::DOMElement* elementNode = static_cast<xercesc::DOMElement*>(node);
	  xercesc::DOMNamedNodeMap *attrMap = node->getAttributes();
	  xercesc::DOMNode* namednode = attrMap->getNamedItem(XMLString::transcode("nx"));
	  char* temp = xercesc::XMLString::transcode(namednode->getNodeValue());
	  ncells.append(atoi(temp));
	  XMLString::release(&temp);
	  namednode = attrMap->getNamedItem(XMLString::transcode("ny"));
	  temp = xercesc::XMLString::transcode(namednode->getNodeValue());
	  ncells.append(atoi(temp));
	  XMLString::release(&temp);
	  namednode = attrMap->getNamedItem(XMLString::transcode("nz"));
	  temp = xercesc::XMLString::transcode(namednode->getNodeValue());
	  ncells.append(atoi(temp));
	  XMLString::release(&temp);
          mesh_list.set<Teuchos::Array<int> >("Number of Cells",ncells);

	  // get Box - generalize
	  Teuchos::Array<double> low;
	  Teuchos::Array<double> high;
	  char* char_array;
	  nodeList = elementGen->getElementsByTagName( XMLString::transcode("box"));
	  node = nodeList->item(0);
	  elementNode = static_cast<xercesc::DOMElement*>(node);
	  temp = xercesc::XMLString::transcode(elementNode->getAttribute(
		  xercesc::XMLString::transcode("low_coordinates")));
	  // translate to array
	  char_array = strtok(temp,",");
	  low.append(atof(char_array));
	  char_array = strtok(NULL,",");
	  low.append(atof(char_array));
	  char_array = strtok(NULL,",");
	  low.append(atof(char_array));
	  XMLString::release(&temp);
	  temp = xercesc::XMLString::transcode(elementNode->getAttribute(
		  xercesc::XMLString::transcode("high_coordinates")));
	  // translate to array
	  char_array = strtok(temp,",");
	  high.append(atof(char_array));
	  char_array = strtok(NULL,",");
	  high.append(atof(char_array));
	  char_array = strtok(NULL,",");
	  high.append(atof(char_array));
	  XMLString::release(&temp);
          mesh_list.set<Teuchos::Array<double> >("Domain Low Corner",low);
          mesh_list.set<Teuchos::Array<double> >("Domain High Corner",high);
	}
	if (strcmp(tagname,"file")==0) {
	  file = true;
	  generate = false;
	  char* filename = XMLString::transcode(currentNode->getTextContent());
          mesh_list.set<std::string>("File",filename);
	}
      }
    }

    if (generate) {
      if (unstructured) {
        list.sublist("Unstructured").sublist("Generate Mesh").sublist("Uniform Structured") = mesh_list;
        if (strcmp(framework,"mstk")==0) {
          list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MSTK");
        } else if (strcmp(framework,"moab")==0) {
          list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MOAB");
        } else if (strcmp(framework,"simple")==0) {
          list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","Simple");
        } else if (strcmp(framework,"stk::mesh")==0) {
          list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","stk::mesh");
        }
      } else {
        list.sublist("Structured") = mesh_list;
      }
    }
    else if (file) {
      mesh_list.set<std::string>("Format",framework);
      list.sublist("Unstructured").sublist("Read Mesh File") = mesh_list;
    }
    else {
      // bad mesh, again if validated shouldn't need this
    }
  }
  else {
    // framework didn't exist, report an error -> already validated, shouldn't need this
  }

  XMLString::release(&framework);

  isUnstr_ = unstructured;
  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_execution_controls(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list ) {

  Teuchos::ParameterList list;

  xercesc::DOMNodeList* nodeList;
  xercesc::DOMNode* nodeTmp;
  xercesc::DOMNode* nodeAttr;
  xercesc::DOMNamedNodeMap* attrMap;
  char* tagName;
  char* attrName;
  char* textContent;
  std::string meshbase;

  // This actually includes: process kernels, execution controls, and numerical controls
  // all three map back to the old exection controls
  
  if (isUnstr_) {
    meshbase = std::string("Unstructured Algorithm");
  } else {
    meshbase = std::string("Structured Algorithm");
  }

  // get execution contorls node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("execution_controls"));
  Teuchos::ParameterList ecsPL;
  Teuchos::ParameterList defPL;
  bool hasSteady = false;
  bool hasTrans = false;

  for (int i=0; i<nodeList->getLength(); i++) {
    xercesc::DOMNode* ecNode = nodeList->item(i);
    if (xercesc::DOMNode::ELEMENT_NODE == ecNode->getNodeType()) {
      //loop over children
      xercesc::DOMNodeList* children = ecNode->getChildNodes();
      for (int j=0; j<children->getLength(); j++) {
        xercesc::DOMNode* currentNode = children->item(j) ;
        if (xercesc::DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
    	  char* tagname = xercesc::XMLString::transcode(currentNode->getNodeName());

          if (strcmp(tagname,"verbosity")==0) {
              attrMap = currentNode->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("level"));
              textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
              list.set<std::string>("Verbosity",textContent);
              XMLString::release(&textContent);

	  } else if (strcmp(tagname,"execution_control_defaults")==0) {
              attrMap = currentNode->getAttributes();
              for (int k=0; k<attrMap->getLength(); k++) {
		nodeAttr = attrMap->item(k);
		attrName =xercesc::XMLString::transcode(nodeAttr->getNodeName());
		textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		defPL.set<std::string>(attrName,textContent);
		if (strcmp(attrName,"mode")==0) {
		  if (strcmp(textContent,"steady")==0) {
		    hasSteady = true;
		  } else {
		    hasTrans = true;
		  }
		}
	      }
	  } else if (strcmp(tagname,"execution_control")==0) {
              Teuchos::ParameterList ecPL;
              attrMap = currentNode->getAttributes();
	      char* name;
              for (int k=0; k<attrMap->getLength(); k++) {
		nodeAttr = attrMap->item(k);
		attrName =xercesc::XMLString::transcode(nodeAttr->getNodeName());
		textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		ecPL.set<std::string>(attrName,textContent);
		if (strcmp(attrName,"start")==0) name=textContent;
		if (strcmp(attrName,"mode")==0) {
		  if (strcmp(textContent,"steady")==0) {
		    hasSteady = true;
		  } else {
		    hasTrans = true;
		  }
		}
	      }
	      ecsPL.sublist(name) = ecPL;
	  }
	}
      }
    }
  }

  // Now, go back and sort things out
  bool haveSSF = false; // have steady-steady incr/red factors for later
  bool haveTF = false;  // have transient incr/red factors for later
  
  // Steady case
  if (hasSteady && !hasTrans) {
    Teuchos::ParameterList steadyPL;
    // look for values from default list
    // if not there, grab from ec list
    std::string value;
    std::string method;
    bool gotValue;
    if (defPL.isParameter("start")) {
      value = defPL.get<std::string>("start");
      gotValue = true;
    } else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
        if (ecsPL.sublist(it->first).isParameter("start")) {
          value = ecsPL.sublist(it->first).get<std::string>("start");
          gotValue = true;
	}
      }
    }
    if (gotValue) {
      double time = get_time_value(value, def_list);
      steadyPL.set<double>("Start",time);
      gotValue = false;
    } else {
      // ERROR - for unstructured, optional for structured;
    }
    if (defPL.isParameter("end")) {
      value = defPL.get<std::string>("end");
      gotValue = true;
    } else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
        if (ecsPL.sublist(it->first).isParameter("end")) {
          value = ecsPL.sublist(it->first).get<std::string>("end");
          gotValue = true;
	}
      }
    }
    if (gotValue) {
      double time = get_time_value(value, def_list);
      steadyPL.set<double>("End",time);
      gotValue = false;
    } else {
      // ERROR ;
    }
    if (defPL.isParameter("init_dt")) {
      value = defPL.get<std::string>("init_dt");
      gotValue = true;
    } else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
        if (ecsPL.sublist(it->first).isParameter("init_dt")) {
          value = ecsPL.sublist(it->first).get<std::string>("init_dt");
          gotValue = true;
	}
      }
    }
    if (gotValue) {
      steadyPL.set<double>("Initial Time Step",atof(value.c_str()));
      gotValue = false;
    } else {
      // ERROR ;
    }
    if (defPL.isParameter("method")) {
      method = defPL.get<std::string>("method");
      gotValue = true;
    } else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
        if (ecsPL.sublist(it->first).isParameter("method")) {
          method = ecsPL.sublist(it->first).get<std::string>("method");
          gotValue = true;
	}
      }
    }
    if (gotValue && strcmp(method.c_str(),"picard")==0) {
      steadyPL.set<bool>("Use Picard","true");
      gotValue = false;
    } else {
      // ERROR ;
    }
    if (defPL.isParameter("reduction_factor")) {
      haveSSF = true;
    } else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
        if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
          value = ecsPL.sublist(it->first).get<std::string>("reduction_factor");
          haveSSF = true;
	}
      }
    }
    if (defPL.isParameter("increase_factor")) {
      haveSSF = true;
    } else {
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
        if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
          haveSSF = true;
	}
      }
    }
    list.sublist("Time Integration Mode").sublist("Steady") = steadyPL;

  } else {
    if (!hasSteady) {
    // Transient case
      Teuchos::ParameterList transPL;
      // loop over ecs to set up, TPC lists
      Teuchos::Array<double> start_times;
      Teuchos::Array<double> init_steps;
      Teuchos::Array<double> max_steps;
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	bool gotValue;
	std::string Value;
	if (ecsPL.sublist(it->first).isParameter("start")) {
            Value = ecsPL.sublist(it->first).get<std::string>("start");
            gotValue = true;
	} else {
            Value = defPL.get<std::string>("start");
            gotValue = true;
	}
        if (gotValue) {
	    double time = get_time_value(Value,def_list);
	    start_times.append(time);
	    gotValue = false;
	}
	if (ecsPL.sublist(it->first).isParameter("end")) {
            Value = ecsPL.sublist(it->first).get<std::string>("end");
            gotValue = true;
	} else {
            Value = defPL.get<std::string>("end");
            gotValue = true;
	}
        if (gotValue) {
	    double time = get_time_value(Value,def_list);
	    transPL.set<double>("End",time);
	    gotValue = false;
	}
	if (ecsPL.sublist(it->first).isParameter("init_dt")) {
            Value = ecsPL.sublist(it->first).get<std::string>("init_dt");
            gotValue = true;
	} else {
            Value = defPL.get<std::string>("init_dt");
            gotValue = true;
	}
        if (gotValue) {
	    init_steps.append(atof(Value.c_str()));
	    gotValue = false;
	}
	if (ecsPL.sublist(it->first).isParameter("max_dt")) {
            Value = ecsPL.sublist(it->first).get<std::string>("max_dt");
            gotValue = true;
	} else {
	  if (defPL.isParameter("max_dt")) {
            Value = defPL.get<std::string>("max_dt");
            gotValue = true;
	  }
	}
        if (gotValue) {
	    max_steps.append(atof(Value.c_str()));
	    gotValue = false;
	}
	if (ecsPL.sublist(it->first).isParameter("reduction_factor") || defPL.isParameter("reduction_factor")) {
            haveTF = true;
        }
	if (ecsPL.sublist(it->first).isParameter("increase_factor") || defPL.isParameter("increase_factor")) {
            haveTF = true;
        }
      }
      transPL.set<double>("Start",start_times[0]);
      transPL.set<double>("Initial Time Step",init_steps[0]);
      if ( max_steps.length()>0) transPL.set<double>("Maximum Time Step Size",max_steps[0]);
      list.sublist("Time Integration Mode").sublist("Transient") = transPL;
      if (start_times.length() > 1) {
	// to include "Time Period Control" list
        Teuchos::ParameterList tpcPL;
	tpcPL.set<Teuchos::Array<double> >("Start Times",start_times);
	tpcPL.set<Teuchos::Array<double> >("Initial Time Step",init_steps);
	if ( max_steps.length()>0) tpcPL.set<Teuchos::Array<double> >("Maximum Time Step",max_steps);
	list.sublist("Time Period Control") = tpcPL;
      }
    } else {
    // Initialize to Steady case
      Teuchos::Array<double> start_times;
      Teuchos::Array<double> init_steps;
      Teuchos::Array<double> max_steps;
      Teuchos::ParameterList timesPL;
      Teuchos::ParameterList initPL;
      std::string value;
      bool gotValue = false;
      for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
	std::string mode("none");
	if (defPL.isParameter("mode")) mode = defPL.get<std::string>("mode");
        if (ecsPL.sublist(it->first).isParameter("mode")) mode = ecsPL.sublist(it->first).get<std::string>("mode");
	if (strcmp(mode.c_str(),"steady")==0) {
	  if (ecsPL.sublist(it->first).isParameter("start")) {
            value = ecsPL.sublist(it->first).get<std::string>("start");
	    double time = get_time_value(value,def_list);
	    initPL.set<double>("Start",time);
	  }
	  if (ecsPL.sublist(it->first).isParameter("end")) {
            value = ecsPL.sublist(it->first).get<std::string>("end");
	    double time = get_time_value(value,def_list);
	    initPL.set<double>("Switch",time);
	  }
	  if (ecsPL.sublist(it->first).isParameter("init_dt")) {
            value = ecsPL.sublist(it->first).get<std::string>("init_dt");
	    initPL.set<double>("Steady Initial Time Step",atof(value.c_str()));
	  } 
	  if (ecsPL.sublist(it->first).isParameter("method")) {
            value = ecsPL.sublist(it->first).get<std::string>("method");
	    if (strcmp(value.c_str(),"true")==0) 
	      initPL.set<bool>("Use Picard",true);
	  } 
	  if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
            haveSSF = true;
	  } 
	  if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
            haveSSF = true;
	  } 
	} else {
	  if (ecsPL.sublist(it->first).isParameter("start")) {
            value = ecsPL.sublist(it->first).get<std::string>("start");
	    double time = get_time_value(value,def_list);
	    start_times.append(time);
	  }
	  if (ecsPL.sublist(it->first).isParameter("end")) {
            value = ecsPL.sublist(it->first).get<std::string>("end");
	    double time = get_time_value(value,def_list);
	    initPL.set<double>("End",time);
	  }
	  if (ecsPL.sublist(it->first).isParameter("init_dt")) {
            value = ecsPL.sublist(it->first).get<std::string>("init_dt");
	    gotValue = true;
	  } else {
	    if (defPL.isParameter("init_dt")) {
              value = defPL.get<std::string>("init_dt");
	      gotValue = true;
	    }
	  }
	  if (gotValue) {
	    init_steps.append(atof(value.c_str()));
	    gotValue = false;
	  }
	}
      }
      initPL.set<double>("Transient Initial Time Step",init_steps[0]);
      list.sublist("Time Integration Mode").sublist("Initialize To Steady") = initPL;
      timesPL.set<Teuchos::Array<double> >("Start Times",start_times);
      timesPL.set<Teuchos::Array<double> >("Initial Time Step",init_steps);
      list.sublist("Time Period Control") = timesPL;
    }
  }
  /*
    for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
      if (ecsPL.sublist(it->first).isParameter("mode")) {
        std::string comp_mode = ecsPL.sublist(it->first).get<std::string>("mode");
	if (strcmp(comp_mode.c_str(),mode.c_str())==0) {
          // Transient case
	 } else {
          // Initialize to Steady case
          Teuchos::ParameterList itsPL;
	  Teuchos::ParameterList tpcPL;
	  // look for init_dt and max_dt from defaults, then execution list
	  std::string value;
          double def_init_dt;
          bool gotValue;
          if (defPL.isParameter("init_dt")) {
            value = defPL.get<std::string>("init_dt");
            gotValue = true;
          } else if (ecsPL.sublist(it->first).isParameter("init_dt")) {
            value = ecsPL.sublist(it->first).get<std::string>("init_dt");
            gotValue = true;
          }
          if (gotValue) {
            itsPL.set<double>("Transient Initial Time Step",atof(value.c_str()));
            gotValue = false;
	    def_init_dt = atof(value.c_str());
          } else {
            // ERROR ;
          }
	  if (defPL.isParameter("reduction_factor")) {
            haveTF = true;
          }
	  if (defPL.isParameter("increase_factor")) {
            haveTF = true;
          }
          // sitting on steady - grab start, end(switch) 
	  if (ecsPL.sublist(it->first).isParameter("start")) {
            value = ecsPL.sublist(it->first).get<std::string>("start");
	    double time = get_time_value(value,def_list);
	    itsPL.set<double>("Start",time);
	  } else {
	    //ERROR
	  }
	  if (ecsPL.sublist(it->first).isParameter("end")) {
            value = ecsPL.sublist(it->first).get<std::string>("end");
	    double time = get_time_value(value,def_list);
	    itsPL.set<double>("Switch",time);
	  } else {
	    //ERROR
	  }
	  if (ecsPL.sublist(it->first).isParameter("init_dt")) {
            value = ecsPL.sublist(it->first).get<std::string>("init_dt");
	    itsPL.set<double>("Steady Initial Time Step",atof(value.c_str()));
	  } 
	  if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
            haveSSF = true;
	  } 
	  if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
            haveSSF = true;
	  } 
	  // loop over execution list for:
	  //     - build list of start times, initial step times (use default value if none given)
	  Teuchos::Array<double> start_times;
	  Teuchos::Array<double> init_steps;
          for (Teuchos::ParameterList::ConstIterator jt = ecsPL.begin(); jt != ecsPL.end(); ++jt) {
            if (ecsPL.sublist(jt->first).isParameter("mode")) {
	    } else {
	      if (ecsPL.sublist(jt->first).isParameter("start")) {
                value = ecsPL.sublist(jt->first).get<std::string>("start");
	        double time = get_time_value(value,def_list);
		start_times.append(time);
		if (ecsPL.sublist(jt->first).isParameter("init_dt")){
		  value = ecsPL.sublist(jt->first).get<std::string>("init_dt");
		  init_steps.append(atof(value.c_str()));
		} else {
		  init_steps.append(def_init_dt);
		}
	      }
	      if (ecsPL.sublist(jt->first).isParameter("end")) {
                value = ecsPL.sublist(jt->first).get<std::string>("end");
	        double time = get_time_value(value,def_list);
		itsPL.set<double>("End",time);
	      }
	    } 
	  }
	  tpcPL.set<Teuchos::Array<double> >("Start Times",start_times);
	  tpcPL.set<Teuchos::Array<double> >("Initial Time Step",init_steps);
	  list.sublist("Time Integration Mode").sublist("Initialize To Steady") = itsPL;
	  list.sublist("Time Period Control") = tpcPL;
	}
      }
    }
  }
  */
  // set so we don't have to reread transport and chemisty list later
  Teuchos::ParameterList tpkPL, cpkPL;
  bool transportON=false;
  bool chemistryON=false;
  // get process kernels node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("process_kernels"));
  for (int i=0; i<nodeList->getLength(); i++) {
    xercesc::DOMNode* pkNode = nodeList->item(i);
    if (xercesc::DOMNode::ELEMENT_NODE == pkNode->getNodeType()) {
      xercesc::DOMElement* pkElement = static_cast<xercesc::DOMElement*>(pkNode);
      // get flow
      xercesc::DOMNodeList* tmpList = pkElement->getElementsByTagName(XMLString::transcode("flow"));
      attrMap = tmpList->item(0)->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("state"));
      textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
      if (strcmp(textContent,"off")==0){
        list.set<std::string>("Flow Model","Off");
      } else {
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("model"));
        XMLString::release(&textContent);
        textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	if (strcmp(textContent,"saturated")==0) {
          list.set<std::string>("Flow Model","Single Phase");
	} else if (strcmp(textContent,"richards")==0) {
          list.set<std::string>("Flow Model","Richards");
	} else {
          list.set<std::string>("Flow Model",textContent);
	}
      }
      XMLString::release(&textContent);
      // get transport
      tmpList = pkElement->getElementsByTagName(XMLString::transcode("transport"));
      attrMap = tmpList->item(0)->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("state"));
      textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
      if (strcmp(textContent,"off")==0) {
        list.set<std::string>("Transport Model","Off");
      } else {
        list.set<std::string>("Transport Model","On");
	transportON=true;
	//TODO: EIB - now get algorithm option
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("algorithm"));
        textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	if (strcmp(textContent,"explicit first-order")) {
	  tpkPL.set<std::string>("Transport Integration Algorithm","Explicit First-Order");
	} else if (strcmp(textContent,"explicit second-order")) {
	  tpkPL.set<std::string>("Transport Integration Algorithm","Explicit Second-Order");
	}
      }
      XMLString::release(&textContent);
      // get chemisty - TODO: EIB - assuming this will be set to OFF!!!!!
      // NOTE: EIB - old spec options seem to be ON/OFF, algorithm option goes under Numerical Control Parameters
      tmpList = pkElement->getElementsByTagName(XMLString::transcode("chemistry"));
      attrMap = tmpList->item(0)->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("state"));
      textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
      if (strcmp(textContent,"off")==0) {
        list.set<std::string>("Chemistry Model","Off");
      } else {
        list.set<std::string>("Chemistry Model","On");
	chemistryON=true;
	//TODO: EIB - now get process model option
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("process_model"));
        textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	if (strcmp(textContent,"implicit operator split")) {
	  cpkPL.set<double>("max chemistry to transport timestep ratio",atof(textContent));
	}
      }
      XMLString::release(&textContent);
    }
  }
  // get numerical controls node
  Teuchos::ParameterList ncPL;
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("numerical_controls"));
  for (int i=0; i<nodeList->getLength(); i++) {
    xercesc::DOMNode* ncNode = nodeList->item(i);
    if (xercesc::DOMNode::ELEMENT_NODE == ncNode->getNodeType()) {
      xercesc::DOMNodeList* childList = ncNode->getChildNodes();
      for(int j=0; j<childList->getLength(); j++) {
        xercesc::DOMNode* tmpNode = childList->item(j) ;
        if (xercesc::DOMNode::ELEMENT_NODE == tmpNode->getNodeType()) {
          char* nodeName = xercesc::XMLString::transcode(tmpNode->getNodeName());
          if (strcmp(nodeName,"steady-state_controls")==0) {
	    Teuchos::ParameterList ssPL;
	    // check for incr/red factors from execution_controls first
	    if (haveSSF) {
	      // look in defPL first
	      if (defPL.isParameter("mode")) {
                std::string mode = defPL.get<std::string>("mode");
		if (strcmp(mode.c_str(),"steady")==0) {
		    if (defPL.isParameter("reduction_factor")) {
		      std::string value = defPL.get<std::string>("reduction_factor");
                      ssPL.set<double>("steady time step reduction factor",atof(value.c_str()));
		    }
		    if (defPL.isParameter("increase_factor")) {
                      std::string value = defPL.get<std::string>("increase_factor");
                      ssPL.set<double>("steady time step increase factor",atof(value.c_str()));
		    }
		} else {
		  // look in ecsPL
                  for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
                    if (ecsPL.sublist(it->first).isParameter("mode")) {
                      std::string mode = ecsPL.sublist(it->first).get<std::string>("mode");
	              if (strcmp(mode.c_str(),"steady")==0) {
		        if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
                          std::string value = ecsPL.sublist(it->first).get<std::string>("reduction_factor");
                          ssPL.set<double>("steady time step reduction factor",atof(value.c_str()));
		        }
		        if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
                          std::string value = ecsPL.sublist(it->first).get<std::string>("increase_factor");
                          ssPL.set<double>("steady time step increase factor",atof(value.c_str()));
		        }
		      }
		    }
		  }
		}
	      }
	    }
            // loop through children and deal with them
            xercesc::DOMNodeList* children = tmpNode->getChildNodes();
            for (int k=0; k<children->getLength(); k++) {
              xercesc::DOMNode* currentNode = children->item(k) ;
              if (xercesc::DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
    	        char* tagname = xercesc::XMLString::transcode(currentNode->getNodeName());
                if (strcmp(tagname,"min_iterations")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    ssPL.set<int>("steady min iterations",atoi(textContent));
                    XMLString::release(&textContent);
                } else if (strcmp(tagname,"max_iterations")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    ssPL.set<int>("steady max iterations",atoi(textContent));
                    XMLString::release(&textContent);
                } else if (strcmp(tagname,"max_preconditioner_lag_iterations")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    ssPL.set<int>("steady max preconditioner lag iterations",atoi(textContent));
                    XMLString::release(&textContent);
                } else if (strcmp(tagname,"nonlinear_tolerance")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    ssPL.set<double>("steady nonlinear tolerance",atoi(textContent));
                    XMLString::release(&textContent);
                } else if (strcmp(tagname,"error_rel_tol")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    ssPL.set<double>("steady error rel tol",atoi(textContent));
                    XMLString::release(&textContent);
                } else if (strcmp(tagname,"error_abs_tol")==0) {
                    textContent = XMLString::transcode(currentNode->getTextContent());
                    ssPL.set<double>("steady error abs tol",atoi(textContent));
                    XMLString::release(&textContent);
                } else if (strcmp(tagname,"pseudo_time_integrator")==0) {
                    Teuchos::ParameterList ptiPL;
                    xercesc::DOMNodeList* kids = currentNode->getChildNodes();
                    for (int l=0; l<kids->getLength(); l++) {
                        xercesc::DOMNode* curNode = kids->item(l) ;
                        if (xercesc::DOMNode::ELEMENT_NODE == curNode->getNodeType()) {
                            char* tag = xercesc::XMLString::transcode(curNode->getNodeName());
                            if (strcmp(tag,"method")==0) {
                                textContent = XMLString::transcode(curNode->getTextContent());
                                ptiPL.set<std::string>("pseudo time integrator time integration method",textContent);
                                XMLString::release(&textContent);
                            } else if (strcmp(tag,"preconditioner")==0) {
                                textContent = XMLString::transcode(curNode->getTextContent());
                                ptiPL.set<std::string>("pseudo time integrator preconditioner",textContent);
                                XMLString::release(&textContent);
                            } else if (strcmp(tag,"linear_solver")==0) {
                                textContent = XMLString::transcode(curNode->getTextContent());
                                ptiPL.set<std::string>("pseudo time integrator linear solver",textContent);
                                XMLString::release(&textContent);
                            } else if (strcmp(tag,"control_options")==0) {
                                textContent = XMLString::transcode(curNode->getTextContent());
                                ptiPL.set<std::string>("pseudo time integrator error control options",textContent);
                                XMLString::release(&textContent);
                            } else if (strcmp(tag,"divergent_max_iterations")==0) {
                                textContent = XMLString::transcode(curNode->getTextContent());
                                ptiPL.set<double>("pseudo time integrator picard maximum number of iterations",atof(textContent));
                                XMLString::release(&textContent);
                            } else if (strcmp(tag,"clipping_saturation")==0) {
                                textContent = XMLString::transcode(curNode->getTextContent());
                                ptiPL.set<double>("pseudo time integrator clipping saturation value",atof(textContent));
                                XMLString::release(&textContent);
                            } else if (strcmp(tag,"convergence_tolerance")==0) {
                                textContent = XMLString::transcode(curNode->getTextContent());
                                ptiPL.set<double>("pseudo time integrator picard convergence tolerance",atof(textContent));
                                XMLString::release(&textContent);
                            } else if (strcmp(tag,"initialize_with_darcy")==0) {
                                textContent = XMLString::transcode(curNode->getTextContent());
                                ptiPL.set<std::string>("pseudo time integrator initialize with darcy",textContent);
                                XMLString::release(&textContent);
                            }
                        }
                    }
                    ssPL.sublist("Steady-State Psuedo-Time Implicit Solver") = ptiPL;
                }
              }
	    }
            list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Steady-State Implicit Time Integration") = ssPL;
	  }
          else if (strcmp(nodeName,"transient_controls")==0) {
	    Teuchos::ParameterList tcPL;
	    // check for incr/red factors from execution_controls first
	    if (haveTF) {
	      // look in defPL first
	      if (defPL.isParameter("mode")) {
                std::string mode = defPL.get<std::string>("mode");
		if (strcmp(mode.c_str(),"transient")==0) {
		    if (defPL.isParameter("reduction_factor")) {
		      std::string value = defPL.get<std::string>("reduction_factor");
                      tcPL.set<double>("transient time step reduction factor",atof(value.c_str()));
		    }
		    if (defPL.isParameter("increase_factor")) {
                      std::string value = defPL.get<std::string>("increase_factor");
                      tcPL.set<double>("transient time step increase factor",atof(value.c_str()));
		    }
		} else {
		  // look in ecsPL
                  for (Teuchos::ParameterList::ConstIterator it = ecsPL.begin(); it != ecsPL.end(); ++it) {
                    if (ecsPL.sublist(it->first).isParameter("mode")) {
                      std::string mode = ecsPL.sublist(it->first).get<std::string>("mode");
	              if (strcmp(mode.c_str(),"transient")==0) {
		        if (ecsPL.sublist(it->first).isParameter("reduction_factor")) {
                          std::string value = ecsPL.sublist(it->first).get<std::string>("reduction_factor");
                          tcPL.set<double>("transient time step reduction factor",atof(value.c_str()));
		        }
		        if (ecsPL.sublist(it->first).isParameter("increase_factor")) {
                          std::string value = ecsPL.sublist(it->first).get<std::string>("increase_factor");
                          tcPL.set<double>("transient time step increase factor",atof(value.c_str()));
		        }
		      }
		    }
		  }
		}
	      }
	    }
	    // grab integration method, then loop through it's children
            xercesc::DOMElement* tcElement = static_cast<xercesc::DOMElement*>(tmpNode);
            xercesc::DOMNodeList* tmpList = tcElement->getElementsByTagName(XMLString::transcode("integration_method"));
            if (tmpList->getLength()>0) {
              xercesc::DOMNode* tmpNode = tmpList->item(0);
              xercesc::DOMNodeList* children = tmpNode->getChildNodes();
              for (int k=0; k<children->getLength(); k++) {
              xercesc::DOMNode* currentNode = children->item(k) ;
                if (xercesc::DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
    	          char* tagname = xercesc::XMLString::transcode(currentNode->getNodeName());
	          if (strcmp(tagname,"convergence_criteria")==0) {
		    // do nothing, don't see where these map!!!
		    /*
                    xercesc::DOMNamedNodeMap* attrMap = currentNode->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("error_rel_tol"));
                    textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
                    tcPL.set<>(,textContent);
                    XMLString::release(&textContent);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("error_abs_tol"));
                    textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
                    tcPL.set<>(,textContent);
                    XMLString::release(&textContent);
		    */
                  } else if (strcmp(tagname,"nonlinear_solver_parameters")==0) {
                    xercesc::DOMElement* curElem = static_cast<xercesc::DOMElement*>(currentNode);
		    if (curElem->hasAttribute(xercesc::XMLString::transcode("min_iterations"))){
		      textContent = xercesc::XMLString::transcode(
				      curElem->getAttribute(xercesc::XMLString::transcode("min_iterations")));
                      tcPL.set<int>("transient min iterations",atoi(textContent));
                      XMLString::release(&textContent);
		    }
		    if (curElem->hasAttribute(xercesc::XMLString::transcode("max_iterations"))){
		      textContent = xercesc::XMLString::transcode(
				      curElem->getAttribute(xercesc::XMLString::transcode("max_iterations")));
                      tcPL.set<int>("transient max iterations",atoi(textContent));
                      XMLString::release(&textContent);
		    }
		    if (curElem->hasAttribute(xercesc::XMLString::transcode("limit_iterations"))){
		      textContent = xercesc::XMLString::transcode(
				      curElem->getAttribute(xercesc::XMLString::transcode("limit_iterations")));
                      tcPL.set<int>("transient nonlinear tolerance",atoi(textContent));
                      XMLString::release(&textContent);
		    }
		    if (curElem->hasAttribute(xercesc::XMLString::transcode("nonlinear_tolerance"))){
		      textContent = xercesc::XMLString::transcode(
				      curElem->getAttribute(xercesc::XMLString::transcode("nonlinear_tolerance")));
                      tcPL.set<int>("transient nonlinear tolerance",atoi(textContent));
                      XMLString::release(&textContent);
		    }
		    if (curElem->hasAttribute(xercesc::XMLString::transcode("max_divergent_iterations"))){
		      textContent = xercesc::XMLString::transcode(
				      curElem->getAttribute(xercesc::XMLString::transcode("max_divergent_iterations")));
                      tcPL.set<int>("transient max divergent iterations",atoi(textContent));
                      XMLString::release(&textContent);
		    }
		    if (curElem->hasAttribute(xercesc::XMLString::transcode("max_preconditioner_lag"))){
		      textContent = xercesc::XMLString::transcode(
				      curElem->getAttribute(xercesc::XMLString::transcode("max_preconditioner_lag")));
                      tcPL.set<int>("transient max preconditioner lag iterations",atoi(textContent));
                      XMLString::release(&textContent);
		    }
                  }
                }
              }
            }
            list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Transient Implicit Time Integration") = tcPL;
          }
          else if (strcmp(nodeName,"linear_solver")==0) {
	    Teuchos::ParameterList lsPL;
	    Teuchos::ParameterList pcPL;
	    bool usePCPL=false;
            // loop through children and deal with them
            xercesc::DOMNodeList* children = tmpNode->getChildNodes();
            for (int k=0; k<children->getLength(); k++) {
              xercesc::DOMNode* currentNode = children->item(k) ;
              if (xercesc::DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
    	        char* tagname = xercesc::XMLString::transcode(currentNode->getNodeName());
                if (strcmp(tagname,"method")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        lsPL.set<std::string>("linear solver preconditioner",textContent);
                        XMLString::release(&textContent);
                } else if (strcmp(tagname,"max_iterations")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        lsPL.set<int>("linear solver maximum iterations",atoi(textContent));
                        XMLString::release(&textContent);
                } else if (strcmp(tagname,"tolerance")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        lsPL.set<double>("linear solver tolerance",atof(textContent));
                        XMLString::release(&textContent);
                } else if (strcmp(tagname,"ml_cycle_applications")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Trilinos ML").set<int>("ML cycle applications",atoi(textContent));
                        XMLString::release(&textContent);
			usePCPL = true;
                } else if (strcmp(tagname,"use_hypre_amg")==0) {
                        //TODO: EIB - not sure what to do with this, sublists get created later
                } else if (strcmp(tagname,"use_block_ilu")==0) {
                        //TODO: EIB - not sure what to do with this, sublists get created later
                } else if (strcmp(tagname,"hypre_amg_cycle_applications")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Hypre AMG").set<int>("Hypre AMG cycle applications",atoi(textContent));
                        XMLString::release(&textContent);
			usePCPL = true;
                } else if (strcmp(tagname,"hypre_amg_smoother_sweeps")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Hypre AMG").set<int>("Hypre AMG smoother sweeps",atoi(textContent));
                        XMLString::release(&textContent);
			usePCPL = true;
                } else if (strcmp(tagname,"hypre_amg_tolerance")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Hypre AMG").set<double>("Hypre AMG tolerance",atof(textContent));
                        XMLString::release(&textContent);
			usePCPL = true;
                } else if (strcmp(tagname,"hypre_amg_threshold")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Hypre AMG").set<double>("Hypre AMG strong threshold",atof(textContent));
                        XMLString::release(&textContent);
			usePCPL = true;
                } else if (strcmp(tagname,"ml_smoother_type")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Trilinos ML").set<std::string>("ML smoother type",textContent);
                        XMLString::release(&textContent);
			usePCPL = true;
                } else if (strcmp(tagname,"sub_cycling")==0) {
                        //TODO: EIB - don't know where this goes
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        //pcPL.set<>("",textContent);
                        XMLString::release(&textContent);
                } else if (strcmp(tagname,"transport_sub_cycling")==0 && transportON) {
                        // add this to tpkPL
                        textContent = XMLString::transcode(currentNode->getTextContent());
			if(strcmp(textContent,"true")) {
                            tpkPL.set<bool>("transport subcycling",true);
			} else {
                            tpkPL.set<bool>("transport subcycling",false);
			}
                        XMLString::release(&textContent);
                }
              }
            }
            list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Linear Solver") = lsPL;
	    if (usePCPL)
              list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Preconditioners") = pcPL;
	    if (transportON) list.sublist("Numerical Control Parameters").sublist(meshbase) = tpkPL;
	    if (chemistryON) list.sublist("Numerical Control Parameters").sublist(meshbase) = cpkPL;
          }
          else if (strcmp(nodeName,"nonlinear_solver")==0) {
            // EIB: creating sub for section that doesn't actually exist yet in the New Schema, but does in the Input Spec
	    Teuchos::ParameterList nlsPL;
            xercesc::DOMNodeList* children = tmpNode->getChildNodes();
            for (int k=0; k<children->getLength(); k++) {
              xercesc::DOMNode* currentNode = children->item(k) ;
              if (xercesc::DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
    	        char* tagname = xercesc::XMLString::transcode(currentNode->getNodeName());
                if (strcmp(tagname,"nonlinear_solver_type")==0) {
                   textContent = XMLString::transcode(currentNode->getTextContent());
                   nlsPL.set<std::string>("Nonlinear Solver Type",textContent);
                   XMLString::release(&textContent);
		}
	      }
	    }
            list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Nonlinear Solver") = nlsPL;
          }
        }
      }
    }      
  }
  
  // TODO: EIB - got back and get the transport algorithm and chemisty model names
  
  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_phases(xercesc::DOMDocument* xmlDoc) {

  Teuchos::ParameterList list;

  xercesc::DOMNodeList* nodeList;
  xercesc::DOMNodeList* nodeList2;
  xercesc::DOMNode* nodeTmp;
  xercesc::DOMNode* nodeTmp2;
  xercesc::DOMNode* nodeAttr;
  xercesc::DOMNamedNodeMap* attrMap;
  char* tagName;
  char* textContent;
  char* textContent2;

  // get phases node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("phases"));
  xercesc::DOMNode* nodeEC = nodeList->item(0);
  xercesc::DOMElement* elementEC = static_cast<xercesc::DOMElement*>(nodeEC);

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
    char* phaseName = xercesc::XMLString::transcode(nodeTmp->getNodeName());
    xercesc::DOMNodeList* childern = nodeTmp->getChildNodes();
    for (int i=0; i<childern->getLength(); i++) {
      xercesc::DOMNode* cur = childern->item(i) ;
      if (xercesc::DOMNode::ELEMENT_NODE == cur->getNodeType()) {
        tagName  = xercesc::XMLString::transcode(cur->getNodeName());
        textContent = XMLString::transcode(cur->getTextContent());
	//TODO: NOTE: EIB - skipping EOS, not currently supported
	if (strcmp(tagName,"viscosity")==0){
          list.sublist("Aqueous").sublist("Phase Properties").sublist("Viscosity: Uniform").set<double>("Viscosity",atof(textContent));
	}
	else if (strcmp(tagName,"density")==0) {
          list.sublist("Aqueous").sublist("Phase Properties").sublist("Density: Uniform").set<double>("Density",atof(textContent));
	}
	else if (strcmp(tagName,"dissolved_components")==0) {
	  Teuchos::ParameterList dcPL;
	  Teuchos::Array<double> diffusion;
	  Teuchos::Array<std::string> solutes;
          xercesc::DOMElement* discompElem = static_cast<xercesc::DOMElement*>(cur);
          nodeList2 = discompElem->getElementsByTagName(XMLString::transcode("solutes"));
          if (nodeList2->getLength() > 0) {
            nodeTmp2 = nodeList2->item(0);
            xercesc::DOMNodeList* kids = nodeTmp2->getChildNodes();
            for (int j=0; j<kids->getLength(); j++) {
              xercesc::DOMNode* curKid = kids->item(j) ;
              if (xercesc::DOMNode::ELEMENT_NODE == curKid->getNodeType()) {
                tagName  = xercesc::XMLString::transcode(curKid->getNodeName());
	        if (strcmp(tagName,"solute")==0){
		  // put value in solutes array
                  textContent2 = XMLString::transcode(curKid->getTextContent());
		  solutes.append(textContent2);
                  XMLString::release(&textContent2);
		  // put attribute - coefficient_of_diffusion in diffusion array
	          attrMap = curKid->getAttributes();
                  nodeAttr = attrMap->getNamedItem(XMLString::transcode("coefficient_of_diffusion"));
                  textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		  diffusion.append(atof(textContent2));
                  XMLString::release(&textContent2);
		}
	      }
	    }
	  }
	  dcPL.set<Teuchos::Array<std::string> >("Component Solutes",solutes);
	  list.sublist("Aqueous").sublist("Phase Components").sublist(phaseName) = dcPL;
	}
        XMLString::release(&textContent);
      }
    }
  }

  // get any solid_phase node
  nodeList = elementEC->getElementsByTagName(XMLString::transcode("solid_phase"));
  if (nodeList->getLength() > 0) {
    Teuchos::ParameterList spPL;
    Teuchos::Array<std::string> minerals;
    nodeTmp = nodeList->item(0);
    xercesc::DOMElement* solidElem = static_cast<xercesc::DOMElement*>(nodeTmp);
    nodeList2 = solidElem->getElementsByTagName(XMLString::transcode("minerals"));
    if (nodeList2->getLength() > 0) {
      nodeTmp2 = nodeList2->item(0);
      xercesc::DOMNodeList* kids = nodeTmp2->getChildNodes();
      for (int i=0; i<kids->getLength(); i++) {
        xercesc::DOMNode* curKid = kids->item(i) ;
        if (xercesc::DOMNode::ELEMENT_NODE == curKid->getNodeType()) {
          tagName  = xercesc::XMLString::transcode(curKid->getNodeName());
	  if (strcmp(tagName,"mineral")==0){
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
Teuchos::ParameterList get_regions(xercesc::DOMDocument* xmlDoc) {

  Teuchos::ParameterList list;

  xercesc::DOMNodeList* nodeList;
  xercesc::DOMNode* nodeTmp;
  xercesc::DOMNode* nodeAttr;
  xercesc::DOMNamedNodeMap* attrMap;
  char* tagName;
  char* nodeName;
  char* textContent;
  char* textContent2;
  char* char_array;

  // get regions node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("regions"));
  xercesc::DOMNode* nodeRgn = nodeList->item(0);
  xercesc::DOMElement* elementRgn = static_cast<xercesc::DOMElement*>(nodeRgn);

  // just loop over the children and deal with them as they come
  // new options: comment, region, box, point
  xercesc::DOMNodeList* childern = nodeRgn->getChildNodes();
  for (int i=0; i<childern->getLength(); i++) {
    xercesc::DOMNode* cur = childern->item(i) ;
    if (xercesc::DOMNode::ELEMENT_NODE == cur->getNodeType()) {
      tagName  = xercesc::XMLString::transcode(cur->getNodeName());
      /* NOTE: EIB - geometry doesn't deal with extra comment node
      if (strcmp(tagName,"comments") == 0){
        textContent = XMLString::transcode(cur->getTextContent());
        list.set<std::string>("comments",textContent);
        XMLString::release(&textContent);
      } */
      if  (strcmp(tagName,"region") == 0){
	attrMap = cur->getAttributes();
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
        textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	// deal with children: comments, box/file
        xercesc::DOMNodeList* kids = cur->getChildNodes();
        for (int j=0; j<kids->getLength(); j++) {
          xercesc::DOMNode* curKid = kids->item(j) ;
          if (xercesc::DOMNode::ELEMENT_NODE == curKid->getNodeType()) {
            nodeName  = xercesc::XMLString::transcode(curKid->getNodeName());
	    /*             
	    if (strcmp(tagName,"comments") == 0){
              textContent = XMLString::transcode(curKid->getTextContent());
              list.set<std::string>("comments",textContent);
              XMLString::release(&textContent);
            }
	    */
            if  (strcmp(nodeName,"box") == 0){
	      attrMap = curKid->getAttributes();
	      // get low coord
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("low_coordinates"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      char_array = strtok(textContent2,"(,");
              Teuchos::Array<double> low;
              Teuchos::Array<double> high;
	      low.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      low.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      low.append(atof(char_array));
              list.sublist(textContent).sublist("Region: Box").set<Teuchos::Array<double> >("Low Coordinate",low);
	      XMLString::release(&textContent2);
	      // get high coord
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("high_coordinates"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      char_array = strtok(textContent2,"(,");
	      high.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      high.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      high.append(atof(char_array));
              list.sublist(textContent).sublist("Region: Box").set<Teuchos::Array<double> >("High Coordinate",high);
	      XMLString::release(&textContent2);
	      XMLString::release(&textContent);
	    }
            else if  (strcmp(nodeName,"plane") == 0){
              Teuchos::Array<double> loc;
              Teuchos::Array<double> dir;
	      attrMap = curKid->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("location"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      char_array = strtok(textContent2,"(,");
	      loc.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      loc.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      loc.append(atof(char_array));
              list.sublist(textContent).sublist("Region: Plane").set<Teuchos::Array<double> >("Location",loc);
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("normal"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      char_array = strtok(textContent2,"(,");
	      dir.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      dir.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      dir.append(atof(char_array));
              list.sublist(textContent).sublist("Region: Plane").set<Teuchos::Array<double> >("Direction",dir);
	      XMLString::release(&textContent2);
	    } else if  (strcmp(nodeName,"file") == 0){
	      //TODO: EIB - add file
	    }
	    XMLString::release(&nodeName);
	  }
	}
      }
      else if  (strcmp(tagName,"box") == 0){
	attrMap = cur->getAttributes();
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
        textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	// get low coord
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("low_coordinates"));
        textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	char_array = strtok(textContent2,"(,");
        Teuchos::Array<double> low;
        Teuchos::Array<double> high;
	low.append(atof(char_array));
	char_array = strtok(NULL,",");
	low.append(atof(char_array));
	char_array = strtok(NULL,",");
	low.append(atof(char_array));
        list.sublist(textContent).sublist("Region: Box").set<Teuchos::Array<double> >("Low Coordinate",low);
	XMLString::release(&textContent2);
	// get high coord
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("high_coordinates"));
        textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	char_array = strtok(textContent2,"(,");
	high.append(atof(char_array));
	char_array = strtok(NULL,",");
	high.append(atof(char_array));
	char_array = strtok(NULL,",");
	high.append(atof(char_array));
        list.sublist(textContent).sublist("Region: Box").set<Teuchos::Array<double> >("High Coordinate",high);
	XMLString::release(&textContent2);
	XMLString::release(&textContent);
      }
      else if  (strcmp(tagName,"point") == 0){
	attrMap = cur->getAttributes();
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
        textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("coordinate"));
        textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	char_array = strtok(textContent2,"(,");
        Teuchos::Array<double> coord;
	coord.append(atof(char_array));
	char_array = strtok(NULL,",");
	coord.append(atof(char_array));
	char_array = strtok(NULL,",");
	coord.append(atof(char_array));
        list.sublist(textContent).sublist("Region: Point").set<Teuchos::Array<double> >("Coordinate",coord);
	XMLString::release(&textContent);
	XMLString::release(&textContent2);
      } else if  (strcmp(tagName,"plane") == 0){
	attrMap = cur->getAttributes();
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
        textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("location"));
        textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	char_array = strtok(textContent2,"(,");
        Teuchos::Array<double> loc;
	loc.append(atof(char_array));
	char_array = strtok(NULL,",");
	loc.append(atof(char_array));
	char_array = strtok(NULL,",");
	loc.append(atof(char_array));
        list.sublist(textContent).sublist("Region: Plane").set<Teuchos::Array<double> >("Location",loc);
        nodeAttr = attrMap->getNamedItem(XMLString::transcode("normal"));
        textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
        Teuchos::Array<double> dir;
	char_array = strtok(textContent2,"(,");
	dir.append(atof(char_array));
	char_array = strtok(NULL,",");
	dir.append(atof(char_array));
	char_array = strtok(NULL,",");
	dir.append(atof(char_array));
        list.sublist(textContent).sublist("Region: Plane").set<Teuchos::Array<double> >("Direction",dir);
	XMLString::release(&textContent);
	XMLString::release(&textContent2);
      }
      XMLString::release(&tagName);
    }
  }

  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_materials(xercesc::DOMDocument* xmlDoc) {

  Teuchos::ParameterList list;

  xercesc::DOMNodeList* nodeList;
  xercesc::DOMNode* nodeTmp;
  xercesc::DOMNode* nodeTmp2;
  xercesc::DOMNode* nodeAttr;
  xercesc::DOMNamedNodeMap* attrMap;
  char* tagName;
  char* propName;
  char* propValue;
  char* textContent;
  char* textContent2;
  char* char_array;
  char* attrName;
  char* attrValue;

  // get regions node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("materials"));
  xercesc::DOMNode* nodeMat = nodeList->item(0);
  xercesc::DOMElement* elementMat = static_cast<xercesc::DOMElement*>(nodeMat);

  // just loop over the children and deal with them as they come
  xercesc::DOMNodeList* childern = nodeMat->getChildNodes();
  for (int i=0; i<childern->getLength(); i++) {
    bool cappressON = false;
    Teuchos::ParameterList caplist;
    std::string capname;
    xercesc::DOMNode* cur = childern->item(i) ;
    if (xercesc::DOMNode::ELEMENT_NODE == cur->getNodeType()) {
      attrMap = cur->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
      textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
      Teuchos::ParameterList matlist(textContent);
      xercesc::DOMNodeList* kids = cur->getChildNodes();
      for (int j=0; j<kids->getLength(); j++) {
        xercesc::DOMNode* curkid = kids->item(j) ;
        if (xercesc::DOMNode::ELEMENT_NODE == curkid->getNodeType()) {
            tagName  = xercesc::XMLString::transcode(curkid->getNodeName());
	    if (strcmp("assigned_regions",tagName)==0){
	      //TODO: EIB - if this is more than 1 region -> assuming comma seperated list of strings????
              textContent2 = xercesc::XMLString::transcode(curkid->getTextContent());
	      Teuchos::Array<std::string> regs;
	      char* char_array;
	      char_array = strtok(textContent2,",");
	      while(char_array!=NULL){
		regs.append(char_array);
	        char_array = strtok(NULL,",");
	      }
	      matlist.set<Teuchos::Array<std::string> >("Assigned Regions",regs);
	      XMLString::release(&textContent2);
	    } else if  (strcmp("mechanical_properties",tagName)==0){
              xercesc::DOMNodeList* list = curkid->getChildNodes();
	      // loop over child: deal with porosity and density
              for (int k=0; k<list->getLength(); k++) {
                xercesc::DOMNode* curkiddy = list->item(k) ;
                if (xercesc::DOMNode::ELEMENT_NODE == curkiddy->getNodeType()) {
                  propName  = xercesc::XMLString::transcode(curkiddy->getNodeName());
	          if  (strcmp("porosity",propName)==0){
	            // TODO: EIB - assuming value, implement file later
                    attrMap = curkiddy->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            matlist.sublist("Porosity: Uniform").set<double>("Value",atof(textContent2));
	            XMLString::release(&textContent2);
	          } else if  (strcmp("particle_density",propName)==0){
	            // TODO: EIB - assuming value, implement file later
		    // TODO: EIB - should be check value >= 0.
                    attrMap = curkiddy->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            matlist.sublist("Particle Density: Uniform").set<double>("Value",atof(textContent2));
	            XMLString::release(&textContent2);
	          } else if  (strcmp("specific_storage",propName)==0){
		    // TODO: EIB - not handling file case
                    attrMap = curkiddy->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            matlist.sublist("Specific Storage: Uniform").set<double>("Value",atof(textContent2));
	            XMLString::release(&textContent2);
	          } else if  (strcmp("specific_yield",propName)==0){
		    // TODO: EIB - not handling file case
                    attrMap = curkiddy->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            matlist.sublist("Specific Yield: Uniform").set<double>("Value",atof(textContent2));
	            XMLString::release(&textContent2);
	          } else if  (strcmp("dispersion_tensor",propName)==0){
		    // TODO: EIB - not handling file case
                    attrMap = curkiddy->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("alphaL"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            matlist.sublist("Dispersion Tensor: Uniform Isotropic").set<double>("alphaL",atof(textContent2));
	            XMLString::release(&textContent2);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("alphaT"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            matlist.sublist("Dispersion Tensor: Uniform Isotropic").set<double>("alphaT",atof(textContent2));
	            XMLString::release(&textContent2);
	          } else if  (strcmp("molecular_diffusion",propName)==0){
		    // TODO: EIB - not handling file case
                    attrMap = curkiddy->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            matlist.sublist("Molecular Diffusion: Uniform").set<double>("Value",atof(textContent2));
	            XMLString::release(&textContent2);
	          } else if  (strcmp("tortuosity",propName)==0){
		    // TODO: EIB - not handling file case
                    attrMap = curkiddy->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            matlist.sublist("Tortuosity: Uniform").set<double>("Value",atof(textContent2));
	            XMLString::release(&textContent2);
	          }
		}
	      }
	    } else if  (strcmp("permeability",tagName)==0){
	      // loop over attributes to get x,y,z
	      char *x,*y,*z;
              attrMap = curkid->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("x"));
              x = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("y"));
              y = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("z"));
              z = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      Teuchos::ParameterList perm;
	      if (strcmp(x,y)==0 && strcmp(y,z)==0) {
		perm.set<double>("Value",atof(x));
	        matlist.sublist("Intrinsic Permeability: Uniform") = perm;
	      } else {
		perm.set<double>("x",atof(x));
	        perm.set<double>("y",atof(y));
	        perm.set<double>("z",atof(z));
	        matlist.sublist("Intrinsic Permeability: Anisotropic Uniform") = perm;
	      }
	      XMLString::release(&x);
	      XMLString::release(&y);
	      XMLString::release(&z);
	    } else if  (strcmp("cap_pressure",tagName)==0){
              attrMap = curkid->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("model"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      if  (strcmp("van_genuchten",textContent2)==0){
		  cappressON = true;
                  xercesc::DOMNodeList* paramList= curkid->getChildNodes();
                  for (int k=0; k<paramList->getLength(); k++) {
                    xercesc::DOMNode* paramNode = paramList->item(k) ;
                    if (xercesc::DOMNode::ELEMENT_NODE == paramNode->getNodeType()) {
                      propName  = xercesc::XMLString::transcode(paramNode->getNodeName());
	              if  (strcmp("parameters",propName)==0){
		        attrMap = paramNode->getAttributes();
		        for (int l=0; l<attrMap->getLength(); l++) {
                            xercesc::DOMNode* attr = attrMap->item(l) ;
                            attrName  = xercesc::XMLString::transcode(attr->getNodeName());
                            attrValue  = xercesc::XMLString::transcode(attr->getNodeValue());
			    if (strcmp(attrName,"sr")==0) {
		              caplist.set<double>("Sr",atof(attrValue));
			    } else {
		              caplist.set<double>(attrName,atof(attrValue));
			    }
	                    XMLString::release(&attrName);
	                    XMLString::release(&attrValue);
			}
		      }
		    }
		  }
		  capname = "Capillary Pressure: van Genuchten";
	      }else if  (strcmp("brooks_corey",textContent2)==0){
		  cappressON = true;
                  xercesc::DOMNodeList* paramList= curkid->getChildNodes();
                  for (int k=0; k<paramList->getLength(); k++) {
                    xercesc::DOMNode* paramNode = paramList->item(k) ;
                    if (xercesc::DOMNode::ELEMENT_NODE == paramNode->getNodeType()) {
                      propName  = xercesc::XMLString::transcode(paramNode->getNodeName());
	              if  (strcmp("parameters",propName)==0){
		        attrMap = paramNode->getAttributes();
		        for (int l=0; l<attrMap->getLength(); l++) {
                            xercesc::DOMNode* attr = attrMap->item(l) ;
                            attrName  = xercesc::XMLString::transcode(attr->getNodeName());
                            attrValue  = xercesc::XMLString::transcode(attr->getNodeValue());
			    if (strcmp(attrName,"sr")==0) {
		              caplist.set<double>("Sr",atof(attrValue));
			    } else {
		              caplist.set<double>(attrName,atof(attrValue));
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
	    } else if  (strcmp("rel_perm",tagName)==0){
	      // TODO: EIB - how to handle if cappress=false? ie, caplist not setup yet?
              attrMap = curkid->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("model"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      if (strcmp(textContent2,"burdine")==0) {
	        caplist.set<std::string>("Relative Permeability","Burdine");
              } else if (strcmp(textContent2,"mualem")==0) {
	        caplist.set<std::string>("Relative Permeability","Mualem");
	      }
	      if (strcmp(textContent2,"none")!=0) {
                  xercesc::DOMNodeList* paramList= curkid->getChildNodes();
                  for (int k=0; k<paramList->getLength(); k++) {
                    xercesc::DOMNode* paramNode = paramList->item(k) ;
                    if (xercesc::DOMNode::ELEMENT_NODE == paramNode->getNodeType()) {
                      propName  = xercesc::XMLString::transcode(paramNode->getNodeName());
	              if  (strcmp("optional_krel_smoothing_interval",propName)==0){
                        propValue  = xercesc::XMLString::transcode(paramNode->getTextContent());
		        caplist.set<double>("krel smoothing interval",atof(propValue));
		      }
		    }
		  }
	      }

	    }
	    XMLString::release(&tagName);
	}
      }
      if(cappressON) matlist.sublist(capname) = caplist;
      list.sublist(textContent) = matlist;
      XMLString::release(&textContent2);
    }

  }

  return list;
  
}

/* 
 ******************************************************************
 * Empty
 ******************************************************************
 */
Teuchos::ParameterList get_initial_conditions(xercesc::DOMDocument* xmlDoc) {

  Teuchos::ParameterList list;

  xercesc::DOMNodeList* nodeList;
  xercesc::DOMNode* nodeTmp;
  xercesc::DOMNode* nodeAttr;
  xercesc::DOMNamedNodeMap* attrMap;
  char* tagName;
  char* propName;
  char* textContent;
  char* textContent2;
  char* char_array;
  char* attrName;
  char* attrValue;
  char* phaseName;

  // get regions node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("initial_conditions"));
  xercesc::DOMNode* nodeIC = nodeList->item(0);
  xercesc::DOMElement* elementIC = static_cast<xercesc::DOMElement*>(nodeIC);

  // just loop over the children and deal with them as they come
  xercesc::DOMNodeList* childern = nodeIC->getChildNodes();
  for (int i=0; i<childern->getLength(); i++) {
    xercesc::DOMNode* cur = childern->item(i) ;
    if (xercesc::DOMNode::ELEMENT_NODE == cur->getNodeType()) {
      // get name of IC, then loop over it's children to fill it in
      attrMap = cur->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
      textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
      Teuchos::ParameterList iclist(textContent);
      xercesc::DOMNodeList* IC = cur->getChildNodes();
      for (int j=0; j<IC->getLength(); j++) {
        xercesc::DOMNode* ICNode = IC->item(j) ;
        tagName  = xercesc::XMLString::transcode(ICNode->getNodeName());
        //NOTE: EIB - ignoring comments for now
        if (strcmp(tagName,"assigned_regions")==0) {
	  //TODO: EIB - if this is more than 1 region -> assuming comma seperated list of strings????
          textContent2 = xercesc::XMLString::transcode(ICNode->getTextContent());
	  Teuchos::Array<std::string> regs;
	  char* char_array;
	  char_array = strtok(textContent2,",");
	  while(char_array!=NULL){
	    regs.append(char_array);
	    char_array = strtok(NULL,",");
	  }
	  iclist.set<Teuchos::Array<std::string> >("Assigned Regions",regs);
	  XMLString::release(&textContent2);
        }
        else if (strcmp(tagName,"liquid_phase")==0) {
          //TODO: EIB - deal with liquid phase
          attrMap = ICNode->getAttributes();
          nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
          phaseName = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	  // loop over children, deal with liquid_component, solute_component, geomchemistry
          xercesc::DOMNodeList* compList = ICNode->getChildNodes();
          for (int k=0; k<compList->getLength(); k++) {
            xercesc::DOMNode* compNode = compList->item(k) ;
            char* compName  = xercesc::XMLString::transcode(compNode->getNodeName());
            if (strcmp(compName,"liquid_component")==0) {
	      // loop over children to find pressure
              xercesc::DOMNodeList* childList = compNode->getChildNodes();
              for (int l=0; l<childList->getLength(); l++) {
                xercesc::DOMNode* pressure = childList->item(l) ;
                char* pressName  = xercesc::XMLString::transcode(pressure->getNodeName());
		//TODO: EIB - saturation
                if (strcmp(pressName,"pressure")==0 || strcmp(pressName,"saturation")==0) {
	          // loop over attributes to get info
	          attrMap = pressure->getAttributes();
                  nodeAttr = attrMap->getNamedItem(XMLString::transcode("function"));
                  textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	          Teuchos::ParameterList pressureList;
	          pressureList.set<std::string>("Phase","Aqueous");
	          // build tags based on function value
	          if (strcmp(textContent2,"uniform")==0) {
		    //value
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    attrValue = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		    pressureList.set<double>("Value",atof(attrValue));
	            XMLString::release(&attrValue);
                    if (strcmp(pressName,"pressure")==0 ) {
		      iclist.sublist("IC: Uniform Pressure") = pressureList;
		    } else {
		      iclist.sublist("IC: Uniform Saturation") = pressureList;
		    }
	          } else {
	            Teuchos::Array<double> coord;
	            Teuchos::Array<double> grad;
	            char* char_array;
		    //value
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    attrValue = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		    pressureList.set<double>("Reference Value",atof(attrValue));
	            XMLString::release(&attrValue);
		    //reference_coord
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("reference_coord"));
                    attrValue = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            char_array = strtok(attrValue,"(,");
	            coord.append(atof(char_array));
	            char_array = strtok(NULL,",");
	            coord.append(atof(char_array));
	            char_array = strtok(NULL,",");
	            coord.append(atof(char_array));
		    pressureList.set<Teuchos::Array<double> >("Reference Coordinate",coord);
	            XMLString::release(&attrValue);
		    //gradient
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("gradient"));
                    attrValue = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            char_array = strtok(attrValue,"(,");
	            grad.append(atof(char_array));
	            char_array = strtok(NULL,",");
	            grad.append(atof(char_array));
	            char_array = strtok(NULL,",");
	            grad.append(atof(char_array));
		    pressureList.set<Teuchos::Array<double> >("Gradient Value",grad);
	            XMLString::release(&attrValue);
                    if (strcmp(pressName,"pressure")==0 ) {
		      iclist.sublist("IC: Linear Pressure") = pressureList;
		    } else {
		      iclist.sublist("IC: Linear Saturation") = pressureList;
		    }
	          }
                  XMLString::release(&textContent2);
		}
	      }
			      
	    }
	    else if (strcmp(compName,"solute_component")==0) {
	      char* solName;
	      char* funcType;
	      Teuchos::ParameterList sclist;
	      attrMap = compNode->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
              solName = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      attrMap = compNode->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("function"));
              funcType = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      if (strcmp(funcType,"uniform")==0){
                nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		sclist.sublist("IC: Uniform Concentration").set<double>("Value",atof(textContent2));
	        XMLString::release(&textContent2);
	      }else if (strcmp(funcType,"lineaer")==0){
		// TODO: EIB - currently can't handle this
	      }
	      //TODO: EIB - not added concerntation units, confused by what to add. grab from units?
	      iclist.sublist("Solute IC").sublist("Aqueous").sublist(phaseName).sublist(solName) = sclist;
	      XMLString::release(&solName);
	      XMLString::release(&funcType);
	    }
	    else if (strcmp(compName,"geochemistry")==0) {
              //TODO: EIB - deal with geochemisty later
	    }
	  }
	  XMLString::release(&phaseName);
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
Teuchos::ParameterList get_boundary_conditions(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  Teuchos::ParameterList list;

  xercesc::DOMNodeList* nodeList;
  xercesc::DOMNodeList* BCList;
  xercesc::DOMNode* nodeTmp;
  xercesc::DOMNode* nodeAttr;
  xercesc::DOMNamedNodeMap* attrMap;
  char* tagName;
  char* propName;
  char* phaseName;
  char* textContent;
  char* textContent2;
  char* char_array;
  char* attrName;
  char* attrValue;


  // get BCs node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("boundary_conditions"));
  xercesc::DOMNode* nodeBC = nodeList->item(0);
  xercesc::DOMElement* elementBC = static_cast<xercesc::DOMElement*>(nodeBC);

  // get list of BCs
  BCList = elementBC->getElementsByTagName(XMLString::transcode("boundary_condition"));
  for (int i=0; i<BCList->getLength(); i++) {
    xercesc::DOMNode* cur = BCList->item(i) ;
    if (xercesc::DOMNode::ELEMENT_NODE == cur->getNodeType()) {
      // get name of BC, then loop over it's children to fill it in
      attrMap = cur->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
      textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
      Teuchos::ParameterList bclist(textContent);
      xercesc::DOMNodeList* BC = cur->getChildNodes();
      for (int j=0; j<BC->getLength(); j++) {
        xercesc::DOMNode* BCNode = BC->item(j) ;
        tagName  = xercesc::XMLString::transcode(BCNode->getNodeName());
        //NOTE: EIB - ignoring comments for now
        if (strcmp(tagName,"assigned_regions")==0) {
	  //TODO: EIB - if this is more than 1 region -> assuming comma seperated list of strings????
          textContent2 = xercesc::XMLString::transcode(BCNode->getTextContent());
	  Teuchos::Array<std::string> regs;
	  char* char_array;
	  char_array = strtok(textContent2,",");
	  while(char_array!=NULL){
	    regs.append(char_array);
	    char_array = strtok(NULL,",");
	  }
	  bclist.set<Teuchos::Array<std::string> >("Assigned Regions",regs);
	  XMLString::release(&textContent2);
        }
        else if (strcmp(tagName,"liquid_phase")==0) {
          //TODO: EIB - deal with liquid phase
          attrMap = BCNode->getAttributes();
          nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
          phaseName = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	  // loop over children, deal with liquid_component, solute_component, geomchemistry
          xercesc::DOMNodeList* compList = BCNode->getChildNodes();
          for (int k=0; k<compList->getLength(); k++) {
            xercesc::DOMNode* compNode = compList->item(k) ;
            char* compName  = xercesc::XMLString::transcode(compNode->getNodeName());
            if (strcmp(compName,"liquid_component")==0) {
	      Teuchos::Array<double> vals;
	      Teuchos::Array<double> times;
	      Teuchos::Array<std::string> funcs;
	      std::string bcname;
	      std::string valname;
	      bool hasCoordsys = false;
	      std::string coordsys = "Absolute";
	      // loop over children and deal with bc's
              xercesc::DOMNodeList* bcChildList = compNode->getChildNodes();
              for (int l=0; l<bcChildList->getLength(); l++) {
                xercesc::DOMNode* bcChildNode = bcChildList->item(l) ;
                if (xercesc::DOMNode::ELEMENT_NODE == bcChildNode->getNodeType()) {
                 char* bcChildName  = xercesc::XMLString::transcode(bcChildNode->getNodeName());
		 if (strcmp(bcChildName,"seepage_face")==0){
		  bcname = "BC: Seepage";
		  valname = "Inward Mass Flux";
                  xercesc::DOMElement* bcElem = static_cast<xercesc::DOMElement*>(bcChildNode);
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("function")));
		  if (strcmp(textContent2,"linear")==0) {
		    funcs.append("Linear");
		  } else if (strcmp(textContent2,"constant")==0) {
		    funcs.append("Constant");
		  } else if (strcmp(textContent2,"uniform")==0) {
		    funcs.append("Uniform");
		  }
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("inward_mass_flux")));
                  double time = get_time_value(textContent2, def_list);
		  vals.append(time);
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("start")));
                  time = get_time_value(textContent2, def_list);
		  times.append(time);
		 } else if (strcmp(bcChildName,"no_flow")==0) {
		  bcname = "BC: Zero Flow";
                  xercesc::DOMElement* bcElem = static_cast<xercesc::DOMElement*>(bcChildNode);
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("function")));
		  if (strcmp(textContent2,"linear")==0) {
		    funcs.append("Linear");
		  } else if (strcmp(textContent2,"constant")==0) {
		    funcs.append("Constant");
		  } else if (strcmp(textContent2,"uniform")==0) {
		    funcs.append("Uniform");
		  }
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("start")));
                  double time = get_time_value(textContent2, def_list);
		  times.append(time);
		 } else {
                  xercesc::DOMElement* bcElem = static_cast<xercesc::DOMElement*>(bcChildNode);
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("function")));
		  if (strcmp(textContent2,"linear")==0) {
		    funcs.append("Linear");
		  } else if (strcmp(textContent2,"constant")==0) {
		    funcs.append("Constant");
		  } else if (strcmp(textContent2,"uniform")==0) {
		    funcs.append("Uniform");
		  }
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("value")));
                  double time = get_time_value(textContent2, def_list);
		  vals.append(time);
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("start")));
                  time = get_time_value(textContent2, def_list);
		  times.append(time);

		  // translate boundary condition name here
		  if (strcmp(bcChildName,"inward_mass_flux")==0) {
		    bcname = "BC: Flux";
		    valname = "Inward Mass Flux";
		  } else if (strcmp(bcChildName,"inward_volumetric_flux")==0) {
		    bcname = "BC: Flux";
		    valname = "Outward Volumetric Flux";
		    } else if (strcmp(bcChildName,"outward_mass_flux")==0) {
		    bcname = "BC: Flux";
		    valname = "Inward Mass Flux";
		  } else if (strcmp(bcChildName,"outward_volumetric_flux")==0) {
		    bcname = "BC: Flux";
		    valname = "Outward Volumetric Flux";
		  } else if (strcmp(bcChildName,"uniform_pressure")==0) {
		    bcname = "BC: Uniform Pressure";
		    valname = "Values";
		  } else if (strcmp(bcChildName,"linear_pressure")==0) {
		    bcname = "BC: Linear Pressure";
		    valname = "Values";
		  } else if (strcmp(bcChildName,"hydrostatic")==0) {
		    bcname = "BC: Hydrostatic";
		    valname = "Water Table Height";
		    //TODO: EIB - update if necessary
		    if (bcElem->hasAttribute(xercesc::XMLString::transcode("coordinate_system"))) {
                      textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		                     xercesc::XMLString::transcode("coordinate_system")));
		      hasCoordsys = true;
		      if (strcmp(textContent2,"relative")==0){
			coordsys = "Relative";
		      }
		    }
		  }
		 }
		}
	      }
	      // if len array == 1: add dummy vals
	      if (times.length()==1 && bcname != "BC: Zero Flow" ){
		times.append(times[0]+1.);
		vals.append(vals[0]);
	      }
	      if (times.length()==funcs.length()) funcs.remove(0); //EIB - this is iffy!!!
	        // create a new list here
	        Teuchos::ParameterList newbclist;
	        newbclist.set<Teuchos::Array<double> >("Times",times);
	        newbclist.set<Teuchos::Array<std::string> >("Time Functions",funcs);
		if (bcname != "BC: Zero Flow") newbclist.set<Teuchos::Array<double> >(valname,vals);
		if (bcname == "BC: Hydrostatic" && hasCoordsys) newbclist.set<std::string>("Coordinate System",coordsys);
	        bclist.sublist(bcname) = newbclist;
	      }
            if (strcmp(compName,"solute_component")==0) {
	      char* solName;
	      Teuchos::ParameterList sclist;
	      // loop over elements to build time series, add to list
	      Teuchos::Array<double> vals;
	      Teuchos::Array<double> times;
	      Teuchos::Array<std::string> funcs;
              xercesc::DOMNodeList* acList = compNode->getChildNodes();
              for (int l=0; l<acList->getLength(); l++) {
		//TODO: EIB - not bother to look at nodeName, only aqueous_conc available for now
                xercesc::DOMNode* cur = acList->item(l) ;
                if (xercesc::DOMNode::ELEMENT_NODE == cur->getNodeType()) {
                  xercesc::DOMElement* bcElem = static_cast<xercesc::DOMElement*>(cur);
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("function")));
		  if (strcmp(textContent2,"linear")==0) {
		    funcs.append("Linear");
		  } else if (strcmp(textContent2,"constant")==0) {
		    funcs.append("Constant");
		  } else if (strcmp(textContent2,"uniform")==0) {
		    funcs.append("Uniform");
		  }
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("value")));
                  double time = get_time_value(textContent2, def_list);
		  vals.append(time);
                  textContent2 = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("start")));
                  time = get_time_value(textContent2, def_list);
		  times.append(time);
                  solName = xercesc::XMLString::transcode(bcElem->getAttribute(
		              xercesc::XMLString::transcode("name")));
		}
	      }
	      //TODO: EIB - not added concerntation units, need to grab from units
	      sclist.sublist("BC: Uniform Concentration").set<Teuchos::Array<double> >("Times",times);
	      sclist.sublist("BC: Uniform Concentration").set<Teuchos::Array<std::string> >("Time Functions",funcs);
	      sclist.sublist("BC: Uniform Concentration").set<Teuchos::Array<double> >("Values",vals);
	      bclist.sublist("Solute BC").sublist("Aqueous").sublist(phaseName).sublist(solName) = sclist;
	      XMLString::release(&solName);
	    }
            if (strcmp(compName,"geochemistry")==0) {
              //TODO: EIB - deal with geochemisty later
	    }
	  }
          XMLString::release(&phaseName);
        }
        else if (strcmp(tagName,"solid_phase")==0) {
          //TODO: EIB - deal with solid phase -> mineral, geochemisty
        }
        XMLString::release(&tagName);
      }
      list.sublist(textContent) = bclist;
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
Teuchos::ParameterList get_sources(xercesc::DOMDocument* xmlDoc, Teuchos::ParameterList def_list) {

  Teuchos::ParameterList list;

  // get Sources node
  xercesc::DOMNodeList* nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("sources"));
  if (nodeList->getLength() > 0) {
  xercesc::DOMNode* nodeSC = nodeList->item(0);
  xercesc::DOMElement* elementSC = static_cast<xercesc::DOMElement*>(nodeSC);

  // get list of Source
  xercesc::DOMNodeList* SCList = elementSC->getElementsByTagName(XMLString::transcode("source"));
  for (int i=0; i<SCList->getLength(); i++) {
    xercesc::DOMNode* cur = SCList->item(i) ;
    if (xercesc::DOMNode::ELEMENT_NODE == cur->getNodeType()) {
      xercesc::DOMNamedNodeMap* attrMap = cur->getAttributes();
      xercesc::DOMNode* nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
      char* textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
      Teuchos::ParameterList sclist;
      xercesc::DOMNodeList* SC = cur->getChildNodes();
      for (int j=0; j<SC->getLength(); j++) {
        xercesc::DOMNode* SCNode = SC->item(j) ;
        if (xercesc::DOMNode::ELEMENT_NODE == SCNode->getNodeType()) {
          char* tagName  = xercesc::XMLString::transcode(SCNode->getNodeName());
          if (strcmp(tagName,"assigned_regions")==0) {
            char* textContent2 = xercesc::XMLString::transcode(SCNode->getTextContent());
	    Teuchos::Array<std::string> regs;
	    char* char_array;
	    char_array = strtok(textContent2,",");
	    while(char_array!=NULL){
	      regs.append(char_array);
	      char_array = strtok(NULL,",");
	    }
	    sclist.set<Teuchos::Array<std::string> >("Assigned Regions",regs);
	    XMLString::release(&textContent2);
          } else if (strcmp(tagName,"liquid_phase")==0) {
            attrMap = SCNode->getAttributes();
            xercesc::DOMNodeList* compList = SCNode->getChildNodes();
            for (int k=0; k<compList->getLength(); k++) {
              xercesc::DOMNode* compNode = compList->item(k) ;
              char* compName  = xercesc::XMLString::transcode(compNode->getNodeName());
              if (strcmp(compName,"liquid_component")==0) {
	        Teuchos::Array<double> vals;
	        Teuchos::Array<double> times;
	        Teuchos::Array<std::string> funcs;
	        std::string scname;
                xercesc::DOMNodeList* scChildList = compNode->getChildNodes();
                for (int l=0; l<scChildList->getLength(); l++) {
                  xercesc::DOMNode* scChildNode = scChildList->item(l) ;
                  if (xercesc::DOMNode::ELEMENT_NODE == scChildNode->getNodeType()) {
                    char* scChildName  = xercesc::XMLString::transcode(scChildNode->getNodeName());
		    if (strcmp(scChildName,"volume_weighted")==0){
		     scname = "Source: Volume Weighted";
		    } else if (strcmp(scChildName,"permeability_weighted")==0){
		     scname = "Source: Permeability Weighted";
		    }
                    xercesc::DOMElement* scElem = static_cast<xercesc::DOMElement*>(scChildNode);
                    char* textContent2 = xercesc::XMLString::transcode(scElem->getAttribute(
		                    xercesc::XMLString::transcode("function")));
		    if (strcmp(textContent2,"linear")==0) {
		      funcs.append("Linear");
		    } else if (strcmp(textContent2,"constant")==0) {
		      funcs.append("Constant");
		    } else if (strcmp(textContent2,"uniform")==0) {
		      funcs.append("Uniform");
		    }
                    textContent2 = xercesc::XMLString::transcode(scElem->getAttribute(
		                   xercesc::XMLString::transcode("start")));
                    double time = get_time_value(textContent2, def_list);
		    times.append(time);
                    textContent2 = xercesc::XMLString::transcode(scElem->getAttribute(
		                   xercesc::XMLString::transcode("value")));
                    time = get_time_value(textContent2, def_list);
		    vals.append(time);
		  }
		}
	        if (times.length()==1 ){
	  	  times.append(times[0]+1.);
		  vals.append(vals[0]);
	        }
	        Teuchos::ParameterList newsclist;
	        newsclist.set<Teuchos::Array<double> >("Times",times);
	        newsclist.set<Teuchos::Array<std::string> >("Time Functions",funcs);
	        newsclist.set<Teuchos::Array<double> >("Values",vals);
	        sclist.sublist(scname) = newsclist;
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
Teuchos::ParameterList get_output(xercesc::DOMDocument* xmlDoc) {

  Teuchos::ParameterList list;

  xercesc::DOMNodeList* nodeList;
  xercesc::DOMNodeList* visList;
  xercesc::DOMNodeList* chkList;
  xercesc::DOMNodeList* obsList;
  xercesc::DOMNodeList* tmpList;
  xercesc::DOMNamedNodeMap* attrMap;
  xercesc::DOMNode* tmpNode;
  xercesc::DOMNode* nodeAttr;
  char* textContent;
  char* textContent2;


  // get definitions node - this node MAY exist ONCE
  // this contains any time macros and cycle macros
  // they are stored in the outputs of the old format
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("definitions"));
  if (nodeList->getLength() > 0) {
    xercesc::DOMNode* defNode = nodeList->item(0);
    xercesc::DOMElement* defElement = static_cast<xercesc::DOMElement*>(defNode);
    xercesc::DOMNodeList* macroList = xmlDoc->getElementsByTagName(XMLString::transcode("macros"));
    Teuchos::ParameterList tmPL;
    Teuchos::ParameterList cmPL;
    //loop over children
    xercesc::DOMNodeList* children = macroList->item(0)->getChildNodes();
    for (int i=0; i<children->getLength(); i++) {
      xercesc::DOMNode* currentNode = children->item(i) ;
      if (xercesc::DOMNode::ELEMENT_NODE == currentNode->getNodeType()) {
	char* tagname = xercesc::XMLString::transcode(currentNode->getNodeName());
	if (strcmp(tagname,"time_macro")==0) {
          Teuchos::ParameterList tm_parameter;
          attrMap = currentNode->getAttributes();
          nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
          textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	  // deal differently if "times" or "start-inter-stop"
          xercesc::DOMNodeList* childList = currentNode->getChildNodes();
	  bool isTime = false;
          for (int j=0; j<childList->getLength(); j++) {
            xercesc::DOMNode* timeNode = childList->item(j) ;
            if (xercesc::DOMNode::ELEMENT_NODE == timeNode->getNodeType()) {
	      if (strcmp(xercesc::XMLString::transcode(timeNode->getNodeName()),"time")==0)
		      isTime = true;
	    }
	  }
	  if ( isTime ) {
            Teuchos::Array<double> times;
            for (int j=0; j<childList->getLength(); j++) {
              xercesc::DOMNode* timeNode = childList->item(j) ;
              if (xercesc::DOMNode::ELEMENT_NODE == timeNode->getNodeType()) {
	        char* nodeTxt = xercesc::XMLString::transcode(timeNode->getTextContent());
	        times.append(atof(nodeTxt));
	        XMLString::release(&nodeTxt);
	      }
	    }
	    tm_parameter.set<Teuchos::Array<double> >("Values", times);
	  } else {
            xercesc::DOMElement* curElement = static_cast<xercesc::DOMElement*>(currentNode);
            xercesc::DOMNodeList* curList = curElement->getElementsByTagName(XMLString::transcode("start"));
            tmpNode = curList->item(0);
	    char* nodeTxt = xercesc::XMLString::transcode(tmpNode->getTextContent());
            Teuchos::Array<double> sps;
	    sps.append(atof(nodeTxt));
	    XMLString::release(&nodeTxt);
            curList = curElement->getElementsByTagName(XMLString::transcode("timestep_interval"));
	    if (curList->getLength() >0) {
              tmpNode = curList->item(0);
	      nodeTxt = xercesc::XMLString::transcode(tmpNode->getTextContent());
	      sps.append(atof(nodeTxt));
	      XMLString::release(&nodeTxt);
              curList = curElement->getElementsByTagName(XMLString::transcode("stop"));
	      if (curList->getLength() >0) {
                tmpNode = curList->item(0);
	        nodeTxt = xercesc::XMLString::transcode(tmpNode->getTextContent());
	        sps.append(atof(nodeTxt));
	        XMLString::release(&nodeTxt);
	      } else {
	        sps.append(-1.0);
	      }
	      tm_parameter.set<Teuchos::Array<double> >("Start_Period_Stop", sps);
	    } else {
	      tm_parameter.set<Teuchos::Array<double> >("Values", sps);
	    }
	  }
	  tmPL.sublist(textContent) = tm_parameter;
	} else if (strcmp(tagname,"cycle_macro")==0) {
          Teuchos::ParameterList cm_parameter;
          attrMap = currentNode->getAttributes();
          nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
          textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
          xercesc::DOMElement* curElement = static_cast<xercesc::DOMElement*>(currentNode);
          xercesc::DOMNodeList* curList = curElement->getElementsByTagName(XMLString::transcode("start"));
          tmpNode = curList->item(0);
	  char* nodeTxt = xercesc::XMLString::transcode(tmpNode->getTextContent());
          Teuchos::Array<int> sps;
	  sps.append(atoi(nodeTxt));
	  XMLString::release(&nodeTxt);
          curList = curElement->getElementsByTagName(XMLString::transcode("timestep_interval"));
	  if (curList->getLength() >0) {
            tmpNode = curList->item(0);
	    nodeTxt = xercesc::XMLString::transcode(tmpNode->getTextContent());
	    sps.append(atoi(nodeTxt));
	    XMLString::release(&nodeTxt);
            curList = curElement->getElementsByTagName(XMLString::transcode("stop"));
	    if (curList->getLength() >0) {
              tmpNode = curList->item(0);
	      nodeTxt = xercesc::XMLString::transcode(tmpNode->getTextContent());
	      sps.append(atoi(nodeTxt));
	      XMLString::release(&nodeTxt);
	    } else {
	      sps.append(-1.0);
	    }
	    cm_parameter.set<Teuchos::Array<int> >("Start_Period_Stop", sps);
	  } else {
	    cm_parameter.set<Teuchos::Array<int> >("Values", sps);
	  }
	  cmPL.sublist(textContent) = cm_parameter;
	}
      }
    }
    list.sublist("Time Macros") = tmPL;
    list.sublist("Cycle Macros") = cmPL;
  }

  // get output node - this node must exist ONCE
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("output"));
  xercesc::DOMNode* outNode = nodeList->item(0);
  xercesc::DOMElement* outElement = static_cast<xercesc::DOMElement*>(outNode);

  // get list of vis - this node MAY exist ONCE
  visList = outElement->getElementsByTagName(XMLString::transcode("vis"));
  for (int i=0; i<visList->getLength(); i++) {
    xercesc::DOMNode* curNode = visList->item(i) ;
    if (xercesc::DOMNode::ELEMENT_NODE == curNode->getNodeType()) {
      Teuchos::ParameterList visPL;
      xercesc::DOMElement* curElement = static_cast<xercesc::DOMElement*>(curNode);
      // get base_filename element
      tmpList = curElement->getElementsByTagName(XMLString::transcode("base_filename"));
      tmpNode = tmpList->item(0);
      textContent = XMLString::transcode(tmpNode->getTextContent());
      visPL.set<std::string>("File Name Base",textContent);
      XMLString::release(&textContent);
      // get num_digits element
      tmpList = curElement->getElementsByTagName(XMLString::transcode("num_digits"));
      tmpNode = tmpList->item(0);
      textContent = XMLString::transcode(tmpNode->getTextContent());
      visPL.set<int>("File Name Digits",atoi(textContent));
      XMLString::release(&textContent);
      // get time_macro element
      tmpList = curElement->getElementsByTagName(XMLString::transcode("time_macro"));
      tmpNode = tmpList->item(0);
      textContent = XMLString::transcode(tmpNode->getTextContent());
      Teuchos::Array<std::string> macro;
      macro.append(textContent);
      visPL.set<Teuchos::Array<std::string> >("Time Macro",macro);
      XMLString::release(&textContent);
      list.sublist("Visualization Data") = visPL;
    }
  }

  // get list of checkpoint - this node MAY exist ONCE
  chkList = outElement->getElementsByTagName(XMLString::transcode("checkpoint"));
  for (int i=0; i<chkList->getLength(); i++) {
    xercesc::DOMNode* curNode = chkList->item(i) ;
    if (xercesc::DOMNode::ELEMENT_NODE == curNode->getNodeType()) {
      Teuchos::ParameterList chkPL;
      xercesc::DOMElement* curElement = static_cast<xercesc::DOMElement*>(curNode);
      // get base_filename element
      tmpList = curElement->getElementsByTagName(XMLString::transcode("base_filename"));
      tmpNode = tmpList->item(0);
      textContent = XMLString::transcode(tmpNode->getTextContent());
      chkPL.set<std::string>("File Name Base",textContent);
      XMLString::release(&textContent);
      // get num_digits element
      tmpList = curElement->getElementsByTagName(XMLString::transcode("num_digits"));
      tmpNode = tmpList->item(0);
      textContent = XMLString::transcode(tmpNode->getTextContent());
      chkPL.set<int>("File Name Digits",atoi(textContent));
      XMLString::release(&textContent);
      // get time_macro element
      tmpList = curElement->getElementsByTagName(XMLString::transcode("time_macro"));
      tmpNode = tmpList->item(0);
      textContent = XMLString::transcode(tmpNode->getTextContent());
      Teuchos::Array<std::string> macro;
      macro.append(textContent);
      chkPL.set<Teuchos::Array<std::string> >("Time Macro",macro);
      XMLString::release(&textContent);
      list.sublist("Checkpoint Data") = chkPL;
    }
  }

  // get list of observations - this node MAY exist ONCE
  obsList = outElement->getElementsByTagName(XMLString::transcode("observations"));
  if (obsList->getLength() > 0) {
  xercesc::DOMNode* nodeObs = obsList->item(0);

  xercesc::DOMNodeList* OBList = nodeObs->getChildNodes();
  Teuchos::ParameterList obsPL;
  for (int i=0; i<OBList->getLength(); i++) {
    xercesc::DOMNode* curNode = OBList->item(i) ;
    if (xercesc::DOMNode::ELEMENT_NODE == curNode->getNodeType()) {
      textContent  = xercesc::XMLString::transcode(curNode->getNodeName());
      if (strcmp(textContent,"filename")==0) {
	textContent2 = xercesc::XMLString::transcode(curNode->getTextContent());
        obsPL.set<std::string>("Observation Output Filename",textContent2);
	XMLString::release(&textContent2);
      } else if (strcmp(textContent,"liquid_phase")==0) {
        xercesc::DOMNamedNodeMap* attrMap = curNode->getAttributes();
        xercesc::DOMNode* nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
	char* phaseName = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	// loop over observations
	xercesc::DOMNodeList* childList = curNode->getChildNodes();
        for (int j=0; j<childList->getLength(); j++) {
	  Teuchos::ParameterList obPL;
          xercesc::DOMNode* curObs = childList->item(j) ;
          if (xercesc::DOMNode::ELEMENT_NODE == curObs->getNodeType()) {
            char* obsType = xercesc::XMLString::transcode(curObs->getNodeName());
            if (strcmp(obsType,"aqueous_pressure")==0) {
	      obPL.set<std::string>("Variable","Aqueous pressure");
	    } else if (strcmp(obsType,"integrated_mass")==0) {
	      // TODO: EIB can't find matching version
	      //obPL.set<std::string>("Variable","Aqueous pressure");
            } else if (strcmp(obsType,"volumetric_water_content")==0) {
	      obPL.set<std::string>("Variable","Volumetric water content");
            } else if (strcmp(obsType,"gravimetric_water_content")==0) {
	      obPL.set<std::string>("Variable","Gravimetric water content");
            } else if (strcmp(obsType,"x_aqueous_volumetric_flux")==0) {
	      // TODO: EIB needs double checking
	      obPL.set<std::string>("Variable","X-Aqueous volumetric flux");
            } else if (strcmp(obsType,"y_aqueous_volumetric_flux")==0) {
	      // TODO: EIB needs double checking
	      obPL.set<std::string>("Variable","Y-Aqueous volumetric flux");
            } else if (strcmp(obsType,"z_aqueous_volumetric_flux")==0) {
	      // TODO: EIB needs double checking
	      obPL.set<std::string>("Variable","Z-Aqueous volumetric flux");
            } else if (strcmp(obsType,"material_id")==0) {
	      obPL.set<std::string>("Variable","MaterialID");
            } else if (strcmp(obsType,"hydraulic_head")==0) {
	      obPL.set<std::string>("Variable","Hydraulic Head");
            } else if (strcmp(obsType,"aqueous_mass_flow_rate")==0) {
	      obPL.set<std::string>("Variable","Aqueous mass flow rate");
            } else if (strcmp(obsType,"aqueous_volumetric_flow_rate")==0) {
	      obPL.set<std::string>("Variable","Aqueous volumetric flow rate");
            } else if (strcmp(obsType,"aqueous_conc")==0) {
	      // get solute name
              xercesc::DOMNamedNodeMap* attrMap = curObs->getAttributes();
              xercesc::DOMNode* nodeAttr = attrMap->getNamedItem(XMLString::transcode("name"));
	      char* soluteName = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      std::stringstream name;
	      name<< soluteName << " Aqueous concentration";
	      obPL.set<std::string>("Variable",name.str());
	    }
	    xercesc::DOMNodeList* kidList = curObs->getChildNodes();
            for (int k=0; k<kidList->getLength(); k++) {
              xercesc::DOMNode* curElem = kidList->item(k) ;
              if (xercesc::DOMNode::ELEMENT_NODE == curElem->getNodeType()) {
                char* Elem =  xercesc::XMLString::transcode(curElem->getNodeName());
                char* Value =  xercesc::XMLString::transcode(curElem->getTextContent());
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
                  obPL.set<std::string>("Region",Value);
		} else if (strcmp(Elem,"functional")==0) {
	          if (strcmp(Value,"point")==0) {
	            obPL.set<std::string>("Functional","Observation Data: Point");
	          } else if (strcmp(Value,"integral")==0) {
	            obPL.set<std::string>("Functional","Observation Data: Integral");
	          } else if (strcmp(Value,"mean")==0) {
	            obPL.set<std::string>("Functional","Observation Data: Mean");
	          }
		} else if (strcmp(Elem,"time_macro")==0) {
	          obPL.set<std::string>("Time Macro",Value);
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
	XMLString::release(&phaseName);
      }
      XMLString::release(&textContent);
      list.sublist("Observation Data") = obsPL;
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

//TODO: EIB - get default time unit from units, convert plain time values if not seconds.

double get_time_value(std::string time_value, Teuchos::ParameterList def_list)
{

  double time;

  // Check if time_value is: listed in constants
  if (def_list.sublist("constants").sublist("constants").isSublist(time_value)) {
    time = def_list.sublist("constants").sublist("constants").sublist(time_value).get<double>("value");

  // Check if time_value is: listed in time_constants
  } else if (def_list.sublist("constants").sublist("time_constants").isSublist(time_value)) {
    time = def_list.sublist("constants").sublist("time_constants").sublist(time_value).get<double>("value");

  // Otherwise, it must be a string or typedef_labeled_time (0000.00,y)
  } else {
    char* tmp = strcpy(new char[time_value.size() + 1], time_value.c_str());
    char* char_array = strtok(tmp,";, ");
    time = atof(char_array);
    char_array = strtok(NULL,";, ");
    if (char_array!=NULL) {
      if (strcmp(char_array,"y")==0) { time = time*31577600.0; }
      else if (strcmp(char_array,"d")==0) { time = time*86100.0; }
      else if (strcmp(char_array,"h")==0) { time = time*3600.0; }
    }
    delete[] tmp;
  }

  return time;
}

}
}
