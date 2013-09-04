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
  
  new_list.sublist("General Description") = get_model_description(doc);
  new_list.sublist("Mesh") = get_Mesh(doc);
  new_list.sublist("Domain").set<int>("Spatial Dimension",dimension_);
  new_list.sublist("Execution Control") = get_execution_controls(doc );
  new_list.sublist("Phase Definitions") = get_phases(doc);
  new_list.sublist("Regions") = get_regions(doc);
  new_list.sublist("Material Properties") = get_materials(doc);
  new_list.sublist("Initial Conditions") = get_initial_conditions(doc);
  new_list.sublist("Boundary Conditions") = get_boundary_conditions(doc);
  new_list.sublist("Output") = get_output(doc);
  
  
  xercesc::XMLPlatformUtils::Terminate();

  // return the completely translated input file as a parameter list
  return new_list;
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
Teuchos::ParameterList get_execution_controls(xercesc::DOMDocument* xmlDoc) {

  Teuchos::ParameterList list;

  xercesc::DOMNodeList* nodeList;
  xercesc::DOMNode* nodeTmp;
  xercesc::DOMNode* nodeAttr;
  xercesc::DOMNamedNodeMap* attrMap;
  char* tagName;
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
	      //TODO: EIB - don't understand what to do with these compared with the next list??????
	  } else if (strcmp(tagname,"execution_control")==0) {
              attrMap = currentNode->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("mode"));
              textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      Teuchos::ParameterList modePL;
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("start"));
              char* txtVal = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      modePL.set<double>("Start",atof(txtVal));
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("end"));
              txtVal = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      modePL.set<double>("End",atof(txtVal));
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("init_dt"));
	      if (nodeAttr > 0) {
                  txtVal = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	          modePL.set<double>("Initial Time Step",atof(txtVal));
                  XMLString::release(&txtVal);
	      }
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("method"));
              txtVal = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      if (strcmp(txtVal,"picard")==0){
	          modePL.set<bool>("Use Picard",true);
	      }
              XMLString::release(&txtVal);
	      if (strcmp(textContent,"steady")==0){
		  list.sublist("Time Integration Mode").sublist("Steady") = modePL;
	      } else {
		  list.sublist("Time Integration Mode").sublist("Transient") = modePL;
	      }
	      //NOTE - EIB option "Initialize To Steady" doesn't seen to be implement in new spec yet, or am I wrong???
	  }
	}
      }
    }
  }

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
      }
      XMLString::release(&textContent);
      // get chemisty - TODO: EIB - assuming this will be set to OFF!!!!!
      // NOTE: EIB - old spec options seem to be ON/OFF, algorithm option goes under Numerical Control Parameters
      tmpList = pkElement->getElementsByTagName(XMLString::transcode("transport"));
      attrMap = tmpList->item(0)->getAttributes();
      nodeAttr = attrMap->getNamedItem(XMLString::transcode("state"));
      textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
      if (strcmp(textContent,"off")==0) {
        list.set<std::string>("Chemistry Model","Off");
      } else {
        list.set<std::string>("Chemistry Model","On");
	chemistryON=true;
	//TODO: EIB - now get process model option
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
	        Teuchos::ParameterList ssPL("Steady-State Implicit Time Integration");
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
            list.sublist("Numerical Control Parameters").sublist(meshbase) = ssPL;
	  }
          else if (strcmp(nodeName,"transient_controls")==0) {
	    Teuchos::ParameterList tcPL("Transient Implicit Time Integration");
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
	                xercesc::DOMNamedNodeMap* attrMap = currentNode->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("min_iterations"));
                    textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
                    tcPL.set<int>("transient min iterations",atoi(textContent));
                    XMLString::release(&textContent);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("max_iterations"));
                    textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
                    tcPL.set<int>("transient max iterations",atoi(textContent));
                    XMLString::release(&textContent);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("limit_iterations"));
                    textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
                    tcPL.set<int>("transient limit iterations",atoi(textContent));
                    XMLString::release(&textContent);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("nonlinear_tolerance"));
                    textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
                    tcPL.set<double>("transient nonlinear tolerance",atof(textContent));
                    XMLString::release(&textContent);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("max_divergent_iterations"));
                    textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
                    tcPL.set<int>("transient max divergent iterations",atoi(textContent));
                    XMLString::release(&textContent);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("max_preconditioner_lag"));
                    textContent = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
                    tcPL.set<int>("transient max preconditioner lag iterations",atoi(textContent));
                    XMLString::release(&textContent);
                  }
                }
              }
            }
            list.sublist("Numerical Control Parameters").sublist(meshbase) = tcPL;
          }
          else if (strcmp(nodeName,"linear_solver")==0) {
	    Teuchos::ParameterList lsPL;
	    Teuchos::ParameterList pcPL;
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
                } else if (strcmp(tagname,"use_hypre_amg")==0) {
                        //TODO: EIB - not sure what to do with this, sublists get created later
                } else if (strcmp(tagname,"use_block_ilu")==0) {
                        //TODO: EIB - not sure what to do with this, sublists get created later
                } else if (strcmp(tagname,"hypre_amg_cycle_applications")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Hypre AMG").set<int>("Hypre AMG cycle applications",atoi(textContent));
                        XMLString::release(&textContent);
                } else if (strcmp(tagname,"hypre_amg_smoother_sweeps")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Hypre AMG").set<int>("Hypre AMG smoother sweeps",atoi(textContent));
                        XMLString::release(&textContent);
                } else if (strcmp(tagname,"hypre_amg_tolerance")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Hypre AMG").set<double>("Hypre AMG tolerance",atof(textContent));
                        XMLString::release(&textContent);
                } else if (strcmp(tagname,"hypre_amg_threshold")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Hypre AMG").set<double>("Hypre AMG strong threshold",atof(textContent));
                        XMLString::release(&textContent);
                } else if (strcmp(tagname,"ml_smoother_type")==0) {
                        textContent = XMLString::transcode(currentNode->getTextContent());
                        pcPL.sublist("Trilinos ML").set<std::string>("ML smoother type",textContent);
                        XMLString::release(&textContent);
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
            list.sublist("Numerical Control Parameters").sublist(meshbase).sublist("Preconditioners") = pcPL;
	    if (transportON) list.sublist("Numerical Control Parameters").sublist(meshbase) = tpkPL;
	    if (chemistryON) list.sublist("Numerical Control Parameters").sublist(meshbase) = cpkPL;
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
  xercesc::DOMNode* nodeTmp;
  char* tagName;
  char* textContent;

  // get phases node
  nodeList = xmlDoc->getElementsByTagName(XMLString::transcode("phases"));
  xercesc::DOMNode* nodeEC = nodeList->item(0);
  xercesc::DOMElement* elementEC = static_cast<xercesc::DOMElement*>(nodeEC);

  // get comments node
  nodeList = elementEC->getElementsByTagName(XMLString::transcode("comments"));
  if (nodeList->getLength() > 0) {
    nodeTmp = nodeList->item(0);
    textContent = XMLString::transcode(nodeTmp->getTextContent());
    list.set<std::string>("comments",textContent);
    XMLString::release(&textContent);
  }

  // get liquid_phase node
  nodeList = elementEC->getElementsByTagName(XMLString::transcode("liquid_phase"));
  if (nodeList->getLength() > 0) {
    nodeTmp = nodeList->item(0);
    xercesc::DOMNodeList* childern = nodeTmp->getChildNodes();
    for (int i=0; i<childern->getLength(); i++) {
      xercesc::DOMNode* cur = childern->item(i) ;
      if (xercesc::DOMNode::ELEMENT_NODE == cur->getNodeType()) {
        tagName  = xercesc::XMLString::transcode(cur->getNodeName());
        textContent = XMLString::transcode(cur->getTextContent());
	if (strcmp(tagName,"viscosity")==0){
          list.sublist("Aqueous").sublist("Phase Properties").sublist("Viscosity: Uniform").set<double>("Viscosity",atof(textContent));}
	else if (strcmp(tagName,"density")==0) {
          list.sublist("Aqueous").sublist("Phase Properties").sublist("Density: Uniform").set<double>("Density",atof(textContent));
	}
      }
    }
    //TODO: EIB - get eos, dissolved components
  }

  // get any solid_phase node
  nodeList = elementEC->getElementsByTagName(XMLString::transcode("solid_phase"));
  if (nodeList->getLength() > 0) {
    // do something with them
    //TODO: EIB - solid phases
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
            tagName  = xercesc::XMLString::transcode(curKid->getNodeName());
	    /*             
	    if (strcmp(tagName,"comments") == 0){
              textContent = XMLString::transcode(curKid->getTextContent());
              list.set<std::string>("comments",textContent);
              XMLString::release(&textContent);
            }
	    */
            if  (strcmp(tagName,"box") == 0){
	      attrMap = curKid->getAttributes();
	      // get low coord
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("low_coordinates"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      char_array = strtok(textContent2,",");
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
	      char_array = strtok(textContent2,",");
	      high.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      high.append(atof(char_array));
	      char_array = strtok(NULL,",");
	      high.append(atof(char_array));
              list.sublist(textContent).sublist("Region: Box").set<Teuchos::Array<double> >("High Coordinate",high);
	      XMLString::release(&textContent2);
	      XMLString::release(&textContent);
	    }
            else if  (strcmp(tagName,"file") == 0){
	      //TODO: EIB - add file
	    }
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
	char_array = strtok(textContent2,",");
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
	char_array = strtok(textContent2,",");
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
                    attrMap = curkiddy->getAttributes();
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	            matlist.sublist("Mass Density: Uniform").set<double>("Value",atof(textContent2));
	            XMLString::release(&textContent2);
	          }
		}
	      }
	    } else if  (strcmp("permeability",tagName)==0){
	      // loop over attributes to get x,y,z
	      // TODO: x,y,z translates to Horizontal and Vertical; how do people want this done?????
              attrMap = curkid->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("x"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      matlist.sublist("Intrinsic Permeability: Anisotropic Uniform").set<double>("x",atof(textContent2));
	      XMLString::release(&textContent2);
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("y"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      matlist.sublist("Intrinsic Permeability: Anisotropic Uniform").set<double>("y",atof(textContent2));
	      XMLString::release(&textContent2);
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("z"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      matlist.sublist("Intrinsic Permeability: Anisotropic Uniform").set<double>("z",atof(textContent2));
	      XMLString::release(&textContent2);
	    } else if  (strcmp("cap_pressure",tagName)==0){
 	      //TODO: EIB - for now, van Genuchten only
	      // read model - van Genuchten or none
              attrMap = curkid->getAttributes();
              nodeAttr = attrMap->getNamedItem(XMLString::transcode("x"));
              textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	      if  (strcmp("van_genuchten",propName)==0){
                  xercesc::DOMNode* paramNode= curkid->getFirstChild();
		  attrMap = paramNode->getAttributes();
		  Teuchos::ParameterList caplist;
		  for (int k=0; k<attrMap->getLength(); k++) {
                    xercesc::DOMNode* attr = attrMap->item(i) ;
                    attrName  = xercesc::XMLString::transcode(attr->getNodeName());
                    attrValue  = xercesc::XMLString::transcode(attr->getNodeValue());
		    caplist.set<double>(attrName,atof(attrValue));
	            XMLString::release(&attrName);
	            XMLString::release(&attrValue);
		  }
		  matlist.sublist("Capillary Pressure: van Genuchten") = caplist;
	      }
	      XMLString::release(&textContent2);
	    } else if  (strcmp("rel_perm",tagName)==0){
		    //TODO: EIB - do something 
		    // this adds to cap_pressure tag
	    }
	}
      }
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
                if (strcmp(pressName,"pressure")==0) {
	          // loop over attributes to get info
	          attrMap = pressure->getAttributes();
                  nodeAttr = attrMap->getNamedItem(XMLString::transcode("function"));
                  textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
	          Teuchos::ParameterList pressureList;
	          pressureList.set<std::string>("Phase","Aqueous");
		  //value
                  nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                  attrValue = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		  pressureList.set<double>("Value",atof(attrValue));
	          XMLString::release(&attrValue);
	          // build tags based on function value
	          if (strcmp(textContent2,"uniform")==0) {
		    iclist.sublist("IC: Uniform Pressure") = pressureList;
	          } else {
		  //TODO: EIB - check for reference_coord
		  //TODO: EIB - check for gradient
		    iclist.sublist("IC: Linear Pressure") = pressureList;
	          }
                  XMLString::release(&textContent2);
		}
	      }
			      
	    }
            if (strcmp(compName,"solute_component")==0) {
              //TODO: EIB - deal with solute_component later
	    }
            if (strcmp(compName,"geochemistry")==0) {
              //TODO: EIB - deal with geochemisty later
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
Teuchos::ParameterList get_boundary_conditions(xercesc::DOMDocument* xmlDoc) {

  Teuchos::ParameterList list;

  xercesc::DOMNodeList* nodeList;
  xercesc::DOMNodeList* BCList;
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
	  // loop over children, deal with liquid_component, solute_component, geomchemistry
          xercesc::DOMNodeList* compList = BCNode->getChildNodes();
          for (int k=0; k<compList->getLength(); k++) {
            xercesc::DOMNode* compNode = compList->item(k) ;
            char* compName  = xercesc::XMLString::transcode(compNode->getNodeName());
            if (strcmp(compName,"liquid_component")==0) {
	      // loop over children and deal with bc's
              xercesc::DOMNodeList* bcChildList = compNode->getChildNodes();
              for (int l=0; l<bcChildList->getLength(); l++) {
                xercesc::DOMNode* bcChildNode = bcChildList->item(l) ;
                char* bcChildName  = xercesc::XMLString::transcode(bcChildNode->getNodeName());
		if (strcmp(bcChildName,"hydrostatic")==0) {
		  Teuchos::ParameterList nlist;
	          attrMap = bcChildNode->getAttributes();
		  //TODO: EIB - need to add testing to verify the attributes are actually there!!!
                  nodeAttr = attrMap->getNamedItem(XMLString::transcode("function"));
                  textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		  if (strcmp(textContent2,"constant")==0) {
		    Teuchos::Array<std::string> tmpStr;
		    tmpStr.append("Constant");
		    nlist.set<Teuchos::Array<std::string> >("Time Functions",tmpStr);
	            XMLString::release(&textContent2);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("start"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		    Teuchos::Array<double> tmp;
		    tmp.append(atof(textContent2));
		    tmp.append(atof(textContent2)+1.0);
		    nlist.set<Teuchos::Array<double> >("Times",tmp);
	            XMLString::release(&textContent2);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		    Teuchos::Array<double> tmp2;
		    tmp2.append(atof(textContent2));
		    tmp2.append(atof(textContent2));
	            XMLString::release(&textContent2);
		    nlist.set<Teuchos::Array<double> >("Water Table Height",tmp2);
		  }
		  // TODO: EIB - old has linear, new has uniform ???
	 	  bclist.sublist("BC: Hydrostatic") = nlist;
		} else if (strcmp(bcChildName,"inward_mass_flux")==0) {
		  Teuchos::ParameterList nlist;
	          attrMap = bcChildNode->getAttributes();
                  nodeAttr = attrMap->getNamedItem(XMLString::transcode("function"));
                  textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		  // TODO: EIB - deal with linear, uniform
		  if (strcmp(textContent2,"constant")==0) {
		    Teuchos::Array<std::string> tmpStr;
		    tmpStr.append("Constant");
		    nlist.set<Teuchos::Array<std::string> >("Time Functions",tmpStr);
	            XMLString::release(&textContent2);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("start"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		    Teuchos::Array<double> tmp;
		    tmp.append(atof(textContent2));
		    tmp.append(atof(textContent2)+1.0);
		    nlist.set<Teuchos::Array<double> >("Times",tmp);
	            XMLString::release(&textContent2);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		    Teuchos::Array<double> tmp2;
		    tmp2.append(atof(textContent2));
		    tmp2.append(atof(textContent2));
	            XMLString::release(&textContent2);
		    nlist.set<Teuchos::Array<double> >("Inward Mass Flux",tmp2);
		  }
		  bclist.sublist("BC: Flux") = nlist;
		} else if (strcmp(bcChildName,"inward_volumetric_flux")==0) {
		  Teuchos::ParameterList nlist;
	          attrMap = bcChildNode->getAttributes();
                  nodeAttr = attrMap->getNamedItem(XMLString::transcode("function"));
                  textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		  // TODO: EIB - deal with linear, uniform
		  if (strcmp(textContent2,"constant")==0) {
		    Teuchos::Array<std::string> tmpStr;
		    tmpStr.append("Constant");
		    nlist.set<Teuchos::Array<std::string> >("Time Functions",tmpStr);
	            XMLString::release(&textContent2);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("start"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		    Teuchos::Array<double> tmp;
		    tmp.append(atof(textContent2));
		    tmp.append(atof(textContent2)+1.0);
		    nlist.set<Teuchos::Array<double> >("Times",tmp);
	            XMLString::release(&textContent2);
                    nodeAttr = attrMap->getNamedItem(XMLString::transcode("value"));
                    textContent2 = xercesc::XMLString::transcode(nodeAttr->getNodeValue());
		    Teuchos::Array<double> tmp2;
		    tmp2.append(atof(textContent2));
		    tmp2.append(atof(textContent2));
	            XMLString::release(&textContent2);
		    nlist.set<Teuchos::Array<double> >("Inward Volumetric Flux",tmp2);
		  }
		  bclist.sublist("BC: Flux") = nlist;
		} else if (strcmp(bcChildName,"uniform_pressure")==0) {
	          //TODO: EIB - do something here
		}
                //TODO: EIB - "Outward Volumetric Flux"
	        //TODO: EIB - "Outward Mass Flux"
	      }
            }
            if (strcmp(compName,"solute_component")==0) {
              //TODO: EIB - deal with solute_component later
	    }
            if (strcmp(compName,"geochemistry")==0) {
              //TODO: EIB - deal with geochemisty later
	    }
	  }
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
      visPL.set<std::string>("File Name Digits",textContent);
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
      chkPL.set<std::string>("File Name Digits",textContent);
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
  for (int i=0; i<obsList->getLength(); i++) {
    xercesc::DOMNode* curNode = obsList->item(i) ;
    if (xercesc::DOMNode::ELEMENT_NODE == curNode->getNodeType()) {
      Teuchos::ParameterList obsPL;
      xercesc::DOMElement* curElement = static_cast<xercesc::DOMElement*>(curNode);
      // get filename element 
      tmpList = curElement->getElementsByTagName(XMLString::transcode("filename"));
      tmpNode = tmpList->item(0);
      textContent = XMLString::transcode(tmpNode->getTextContent());
      obsPL.set<std::string>("Observation Output Filename",textContent);
      XMLString::release(&textContent);
      // get liquid phase
      tmpList = curElement->getElementsByTagName(XMLString::transcode("liquid_phase"));
      // get list of observation(s)
      xercesc::DOMElement* tmpElement = static_cast<xercesc::DOMElement*>(tmpList->item(0));
      xercesc::DOMNodeList* childList = tmpElement->getElementsByTagName(XMLString::transcode("observation"));
      for (int j=0; j<childList->getLength(); j++) {
        // loop over observation(s)
        xercesc::DOMNode* childNode = childList->item(j) ;
	// get variable name (attribute)
	xercesc::DOMNamedNodeMap *attrMap = childNode->getAttributes();
	xercesc::DOMNode* namednode = attrMap->getNamedItem(XMLString::transcode("variable"));
	textContent = xercesc::XMLString::transcode(namednode->getNodeValue());
	Teuchos::ParameterList obPL(textContent);
	obPL.set<std::string>("Variable",textContent);
	// get assigned regions (child element)
        xercesc::DOMElement* childElement = static_cast<xercesc::DOMElement*>(childNode);
	xercesc::DOMNodeList* kidList = childElement->getElementsByTagName(XMLString::transcode("assigned_regions"));
	if (kidList->getLength() > 0) {
          textContent2 = xercesc::XMLString::transcode(kidList->item(0)->getTextContent());
	  obPL.set<std::string>("Region",textContent2);
	  XMLString::release(&textContent2);
	}
	// get functional (child element)
	kidList = childElement->getElementsByTagName(XMLString::transcode("functional"));
	if (kidList->getLength() > 0) {
          textContent2 = xercesc::XMLString::transcode(kidList->item(0)->getTextContent());
	  if (strcmp(textContent2,"point")==0) {
	    obPL.set<std::string>("Functional","Observation Data: Point");
	  } else if (strcmp(textContent2,"integral")==0) {
	    obPL.set<std::string>("Functional","Observation Data: Integral");
	  } else if (strcmp(textContent2,"mean")==0) {
	    obPL.set<std::string>("Functional","Observation Data: Mean");
	  }
	  XMLString::release(&textContent2);
	}
	// get time_macro (child element)
	kidList = childElement->getElementsByTagName(XMLString::transcode("time_macro"));
	if (kidList->getLength() > 0) {
          textContent2 = xercesc::XMLString::transcode(kidList->item(0)->getTextContent());
	  obPL.set<std::string>("Time Macro",textContent2);
	  XMLString::release(&textContent2);
	}
	std::stringstream listName;
	listName << "observation-"<<j+1<<":"<<textContent;
	obsPL.sublist(listName.str()) = obPL;
	XMLString::release(&textContent);
      }
      list.sublist("Observation Data") = obsPL;
    }
  }

  return list;
  
}


}
}
