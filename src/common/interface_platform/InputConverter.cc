/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>

// TPLs
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

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/XMLString.hpp"
#include "xercesc/util/PlatformUtils.hpp"
#include "xercesc/parsers/DOMLSParserImpl.hpp"
#include "xercesc/framework/StdOutFormatTarget.hpp"
#include "xercesc/util/OutOfMemoryException.hpp"

// Amanzi's
#include "ErrorHandler.hpp"
#include "InputConverter.hh"


namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Initialization of xercecs document.
****************************************************************** */
void InputConverter::Init(const std::string& xmlfilename)
{
  Teuchos::ParameterList out_list;
  
  XercesDOMParser *parser = new XercesDOMParser();
  parser->setExitOnFirstFatalError(true);
  parser->setValidationConstraintFatal(true);
  parser->setValidationScheme(XercesDOMParser::Val_Never);
  parser->setDoNamespaces(true);

  AmanziErrorHandler* errorHandler = new AmanziErrorHandler();
  parser->setErrorHandler(errorHandler);
  parser->useCachedGrammarInParse(true);
 
  bool errorsOccured = false;

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

  doc_ = parser->getDocument();

  delete errorHandler;
}


/* ******************************************************************
* Returns node tag1->tag2 where both tag1 and tag2 are unique leaves
* of the tree.
****************************************************************** */
DOMNode* InputConverter::getUniqueElementByTagNames_(
    const std::string& tag1, const std::string& tag2, bool& flag)
{
  flag = false;
  DOMNode* node;
  DOMNodeList* node_list = doc_->getElementsByTagName(XMLString::transcode(tag1.c_str()));
  if (node_list->getLength() != 1) return node;

  int ntag2(0);
  DOMNodeList* children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      char* tagname = XMLString::transcode(inode->getNodeName());   
      if (strcmp(tagname, tag2.c_str()) == 0) {
        node = inode;
        ntag2++;
      }
      XMLString::release(&tagname);
    }
  }

  if (ntag2 == 1) flag = true;
  return node;
}


/* ******************************************************************
* Returns node tag1->tag2->tag3 where tag1, tag2 iand tag3 are unique
* leaves of the tree.
****************************************************************** */
DOMNode* InputConverter::getUniqueElementByTagNames_(
    const std::string& tag1, const std::string& tag2, const std::string& tag3, bool& flag)
{
  int ntag2(0), ntag3(0);
  DOMNode* node;

  flag = false;
  DOMNodeList* node_list = doc_->getElementsByTagName(XMLString::transcode(tag1.c_str()));
  if (node_list->getLength() != 1) return node;

  // first leaf
  DOMNodeList* children = node_list->item(0)->getChildNodes();
  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      char* tagname = XMLString::transcode(inode->getNodeName());   
      if (strcmp(tagname, tag2.c_str()) == 0) {
        node = inode;
        ntag2++;
      }
      XMLString::release(&tagname);
    }
  }
  if (ntag2 != 1) return node;

  // second leaf
  children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      char* tagname = XMLString::transcode(inode->getNodeName());   
      if (strcmp(tagname, tag3.c_str()) == 0) {
        node = inode;
        ntag3++;
      }
      XMLString::release(&tagname);
    }
  }
  if (ntag3 == 1) flag = true;

  return node;
}


/* ******************************************************************
* Return node described by the list of consequtive names tags 
* separated by commas. It 
****************************************************************** */
DOMNode* InputConverter::getUniqueElementByTagsString_(
    const std::string& tags, bool& flag)
{
  DOMNode* node;

  flag = false;
  std::vector<std::string> tag_names = CharToStrings_(tags.c_str());
  if (tag_names.size() == 0) return node;

  // get the first node
  DOMNodeList* node_list = doc_->getElementsByTagName(XMLString::transcode(tag_names[0].c_str()));
  if (node_list->getLength() != 1) return node;

  for (int n = 1; n < tag_names.size(); ++n) {
    DOMNodeList* children = node->getChildNodes();
    int ntag(0);
    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
        char* tagname = XMLString::transcode(inode->getNodeName());   
        if (strcmp(tagname, tag_names[n].c_str()) == 0) {
          node = inode;
          ntag++;
        }
        XMLString::release(&tagname);
      }
    }
    if (ntag != 1) return node;
  }

  flag = true;
  return node;
}


/* ******************************************************************
* Returns node tag1->tag2 where both tag1 and tag2 are unique leaves
* of the tree.
****************************************************************** */
DOMNode* InputConverter::getUniqueElementByTagNames_(
    const DOMNode* node1, const std::string& tag2, bool& flag)
{
  flag = false;
  int ntag2(0);
  DOMNode* node;
  DOMNodeList* children = node1->getChildNodes();

  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      char* tagname = XMLString::transcode(inode->getNodeName());   
      if (strcmp(tagname, tag2.c_str()) == 0) {
        node = inode;
        ntag2++;
      }
      XMLString::release(&tagname);
    }
  }

  if (ntag2 == 1) flag = true;
  return node;
}


/* ******************************************************************
* Returns node tag1->tag2 where both tag1 and tag2 are unique leaves
* of the tree.
****************************************************************** */
DOMNode* InputConverter::getUniqueElementByTagNames_(
    const DOMNode* node1, const std::string& tag2, const std::string& tag3, bool& flag)
{
  flag = false;
  int ntag2(0), ntag3(0);
  DOMNode* node;
  DOMNodeList* children = node1->getChildNodes();

  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      char* tagname = XMLString::transcode(inode->getNodeName());   
      if (strcmp(tagname, tag2.c_str()) == 0) {
        node = inode;
        ntag2++;
      }
      XMLString::release(&tagname);
    }
  }
  if (ntag2 != 1) return node;

  // second leaf
  children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      char* tagname = XMLString::transcode(inode->getNodeName());   
      if (strcmp(tagname, tag3.c_str()) == 0) {
        node = inode;
        ntag3++;
      }
      XMLString::release(&tagname);
    }
  }
  if (ntag3 == 1) flag = true;

  return node;
}


/* ******************************************************************
* Return node described by the list of consequtive names tags 
* separated by commas.
****************************************************************** */
DOMNode* InputConverter::getUniqueElementByTagsString_(
    const DOMNode* node1, const std::string& tags, bool& flag)
{
  DOMNode* node;

  flag = false;
  std::vector<std::string> tag_names = CharToStrings_(tags.c_str());
  if (tag_names.size() == 0) return node;

  // get the first node
  node = const_cast<DOMNode*>(node1);

  for (int n = 0; n < tag_names.size(); ++n) {
    DOMNodeList* children = node->getChildNodes();
    int ntag(0);
    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
        char* tagname = XMLString::transcode(inode->getNodeName());   
        if (strcmp(tagname, tag_names[n].c_str()) == 0) {
          node = inode;
          ntag++;
        }
        XMLString::release(&tagname);
      }
    }
    if (ntag != 1) return node;
  }

  flag = true;
  return node;
}


/* ******************************************************************
* Extract atribute of type double.
****************************************************************** */
double InputConverter::GetAttributeValueD_(
    DOMElement* elem, const char* attr_name, bool exception, double default_val)
{
  double val;
  if (elem->hasAttribute(XMLString::transcode(attr_name))) {
    char* text_content = XMLString::transcode(elem->getAttribute(XMLString::transcode(attr_name)));
    val = std::strtod(text_content, NULL);
    XMLString::release(&text_content);
  } else if (!exception) {
    val = default_val;
  } else {
    char* tagname = XMLString::transcode(elem->getNodeName());
    ThrowErrorMissattr_(tagname, "attribute", attr_name, tagname);
  }
  return val;
}


/* ******************************************************************
* Extract atribute of type int.
****************************************************************** */
int InputConverter::GetAttributeValueL_(
    DOMElement* elem, const char* attr_name, bool exception, int default_val)
{
  int val;
  if (elem->hasAttribute(XMLString::transcode(attr_name))) {
    char* text_content = XMLString::transcode(elem->getAttribute(XMLString::transcode(attr_name)));
    val = std::strtol(text_content, NULL, 10);
    XMLString::release(&text_content);
  } else if (! exception) {
    val = default_val;
  } else {
    char* tagname = XMLString::transcode(elem->getNodeName());
    ThrowErrorMissattr_(tagname, "attribute", attr_name, tagname);
  }
  return val;
}


/* ******************************************************************
* Extract atribute of type char*.
* The returned memory should be freeded using XMLString::release.
****************************************************************** */
char* InputConverter::GetAttributeValueC_(
    DOMElement* elem, const char* attr_name, bool exception, char* default_val)
{
  char* val;
  if (elem->hasAttribute(XMLString::transcode(attr_name))) {
    val = XMLString::transcode(elem->getAttribute(XMLString::transcode(attr_name)));
  } else if (!exception) {
    val = default_val;
  } else {
    char* tagname = XMLString::transcode(elem->getNodeName());
    ThrowErrorMissattr_(tagname, "attribute", attr_name, tagname);
  }
  return val;
}


/* ******************************************************************
* Extract attribute of type vector.
****************************************************************** */
std::vector<double> InputConverter::GetAttributeVector_(DOMElement* elem, const char* attr_name)
{
  std::vector<double> val;
  if (elem->hasAttribute(XMLString::transcode(attr_name))) {
    char* text_content = XMLString::transcode(elem->getAttribute(XMLString::transcode(attr_name)));
    val = MakeCoordinates_(text_content);
    XMLString::release(&text_content);
  } else {
    char* tagname = XMLString::transcode(elem->getNodeName());
    ThrowErrorMissattr_(tagname, "attribute", attr_name, tagname);
  }
  return val;
}


/* ******************************************************************
* Find positing in the array.
****************************************************************** */
int InputConverter::GetPosition_(const std::vector<std::string>& names, char* name)
{
  for (int i = 0; i < names.size(); ++i) {
    if (strcmp(names[i].c_str(), name) == 0) return i;
  }

  Errors::Message msg;
  msg << "Amanzi::InputConverter: vector of names (e.g. solutes) has no \"" << name << "\".\n";
  msg << "  Please correct and try again.\n";
  Exceptions::amanzi_throw(msg);

  return -1;
}

/* ******************************************************************
* Converts string of names separated by comma to array of strings.
****************************************************************** */
std::vector<std::string> InputConverter::CharToStrings_(const char* namelist)
{
  char* tmp1 = new char[strlen(namelist) + 1];
  strcpy(tmp1, namelist);

  char* tmp2;
  tmp2 = strtok(tmp1, ",");

  std::vector<std::string> regs;
  while (tmp2 != NULL) {
    std::string str(tmp2);
    boost::algorithm::trim(str);
    regs.push_back(str);
    tmp2 = strtok(NULL, ",");
  }

  delete[] tmp1;
  return regs;
}


/* ******************************************************************
* Empty
****************************************************************** */
double InputConverter::TimeStringToValue_(std::string& time_value)
{
  double time;
  char* tmp = strcpy(new char[time_value.size() + 1], time_value.c_str());
  time = TimeCharToValue_(tmp);
  delete[] tmp;

  return time;
}


/* ******************************************************************
* Get default time unit from units, convert plain time values if not seconds.
****************************************************************** */
double InputConverter::TimeCharToValue_(char* time_value)
{
  double time;
  char* char_array;
  
  char_array = strtok(time_value, ";, ");
  time = std::strtod(char_array, NULL);
  char_array = strtok(NULL, ";,");

  if (char_array != NULL) {
    if (strcmp(char_array, "y") == 0) { 
      time *= 365.25 * 24.0 * 60.0 * 60.0;
    } else if (strcmp(char_array, "d") == 0) {
      time *= 24.0*60.0*60.0;
    } else if (strcmp(char_array, "h") == 0) {
      time *= 60.0*60.0;
    }
  }
  
  return time;
}


/* ******************************************************************
* Converts coordinate string to an array of doubles.
****************************************************************** */
std::vector<double> InputConverter::MakeCoordinates_(char* char_array)
{
  std::vector<double> coords;
  char* tmp;
  tmp = strtok(char_array, "(, ");

  while (tmp != NULL) {
    std::string str(tmp);
    boost::algorithm::trim(str);
    coords.push_back(std::strtod(str.c_str(), NULL));
    tmp = strtok(NULL, ",");
  }

  return coords;
}


/* ******************************************************************
* Empty
****************************************************************** */
std::string InputConverter::TrimString_(char* tmp)
{
  std::string str(tmp);
  boost::algorithm::trim(str);
  return str;
}


/* *******************************************************************
* Generate unified error message for ill-formed element
******************************************************************* */
void InputConverter::ThrowErrorIllformed_(
    const std::string& section, const std::string& type, const std::string& ill_formed)
{
  Errors::Message msg;
  msg << "Amanzi::InputConverter: an error occurred during parsing node \"" << section << "\"\n";
  msg << "  Missing or ill-formed " << type << " for \"" << ill_formed << "\".\n";
  msg << "  Please correct and try again.\n";
  Exceptions::amanzi_throw(msg);
}


/* *****************************************************************************
* Generate unified error message for ill-formed element with options provided
***************************************************************************** */
void InputConverter::ThrowErrorIllformed_(
    const std::string& section, const std::string& type, const std::string& ill_formed, const std::string& options)
{
  Errors::Message msg;
  msg << "Amanzi::InputTranslator: an error occurred during parsing node \"" << section << "\"\n";
  msg << "  Missing or ill-formed " << type << " for \"" << ill_formed << "\"\n";
  msg << "  Valid options are: " << options << "\n";
  msg << "  Please correct and try again.\n" ;
  Exceptions::amanzi_throw(msg);
}


/* *******************************************************************
* Generate unified error message for missing item
******************************************************************* */
void InputConverter::ThrowErrorMissattr_(
    const std::string& section, const std::string& type, const std::string& missing, const std::string& name)
{
  Errors::Message msg;
  msg << "Amanzi::InputConverter: an error occurred during parsing node \"" << section << "\"\n";
  msg << "  No " << type << " \"" << missing << "\" found for \"" << name << "\".\n";
  msg << "  Please correct and try again \n";
  Exceptions::amanzi_throw(msg);
}

}  // namespace AmanziInput
}  // namespace Amanzi
