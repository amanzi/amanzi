/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
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
void InputConverter::Init(const std::string& xmlfilename, bool& valid)
{
  xmlfilename_ = xmlfilename;
  
  parser = new XercesDOMParser();
  parser->setExitOnFirstFatalError(true);
  parser->setValidationConstraintFatal(true);
  parser->setValidationScheme(XercesDOMParser::Val_Never);
  parser->setDoNamespaces(true);
  parser->setCreateCommentNodes(false);

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
  FilterNodes(doc_, "comments");

  // check that XML has version 2.x or version 1.x spec
  MemoryManager mm;
  DOMElement* element = doc_->getDocumentElement();
  valid = (strcmp(mm.transcode(element->getTagName()), "amanzi_input") == 0);

  delete errorHandler;
}


/* ******************************************************************
* Populates protected std::map constants_.
****************************************************************** */
void InputConverter::ParseConstants_()
{
  MemoryManager mm;

  char *tagname, *text;
  DOMNodeList *node_list, *children;

  // process constants: we ignore type of generic constants.
  node_list = doc_->getElementsByTagName(mm.transcode("constants"));
  if (node_list->getLength() == 0) return;

  children = node_list->item(0)->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    
    std::string name, type, value;
    char* tagname = mm.transcode(inode->getNodeName());   
    if (strcmp(tagname, "constant") == 0 ||
        strcmp(tagname, "time_constant") == 0 ||
        strcmp(tagname, "numerical_constant") == 0 ||
        strcmp(tagname, "area_mass_flux_constant") == 0) {
      DOMElement* element = static_cast<DOMElement*>(inode);
      if (element->hasAttribute(mm.transcode("name"))) {
        text = mm.transcode(element->getAttribute(mm.transcode("name")));
        name = text;
      } else {
        ThrowErrorMissattr_("constants", "attribute", "name", "constant");
      }

      if (element->hasAttribute(mm.transcode("value"))) {
        text = mm.transcode(element->getAttribute(mm.transcode("value")));
        value = text;
      } else {
        ThrowErrorMissattr_("constants", "attribute", "value", "constant");
      }

      if (constants_.find(name) != constants_.end()) {
        Errors::Message msg;
        msg << "An error occurred during parsing node \"constants\"\n";
        msg << "Name \"" << name << "\" is repeated.\n";
        msg << "Please correct and try again \n";
        Exceptions::amanzi_throw(msg);
      } 

      constants_[name] = value;
    }
  }
}


/* ******************************************************************
* Filters out all nodes named "filter" starting with node "parent".
****************************************************************** */
void InputConverter::FilterNodes(DOMNode* parent, const std::string& filter)
{
  DOMNodeList* children = parent->getChildNodes();
  MemoryManager mm;
  
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* child = children->item(i);
  
    if (child->getNodeType() == DOMNode::ELEMENT_NODE) {
      if (mm.transcode(child->getNodeName()) == filter) {
        parent->removeChild(child);
      } else {
        FilterNodes(child, filter);
      }
    }
  }
}


/* ******************************************************************
* Return node described by the list of consequtive names tags 
* separated by commas. It 
****************************************************************** */
DOMNode* InputConverter::GetUniqueElementByTagsString_(
    const std::string& tags, bool& flag)
{
  flag = false;

  MemoryManager mm;
  DOMNode* node = NULL;

  std::vector<std::string> tag_names = CharToStrings_(tags.c_str());
  if (tag_names.size() == 0) return node;

  if (tag_names.size() == 1)
  {
    // We only need the top-level (hopefully unique) tag.
    DOMElement* root = doc_->getDocumentElement();
    return GetChildByName_(root, tag_names[0], flag);
  }

  // get the first node
  DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode(tag_names[0].c_str()));
  if (node_list->getLength() != 1) return node;
  node = node_list->item(0);

  for (int n = 1; n < tag_names.size(); ++n) {
    DOMNodeList* children = node->getChildNodes();
    int ntag(0);
    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
        char* tagname = mm.transcode(inode->getNodeName());   
        if (strcmp(tagname, tag_names[n].c_str()) == 0) {
          node = inode;
          ntag++;
        }
      }
    }
    if (ntag != 1) return node;
  }

  flag = true;
  return node;
}


/* ******************************************************************
* Return node described by the list of consequtive names tags 
* separated by commas.
****************************************************************** */
DOMNode* InputConverter::GetUniqueElementByTagsString_(
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
* Returns the child with the given attribute name and value.
****************************************************************** */
DOMElement* InputConverter::GetUniqueChildByAttribute_(
    xercesc::DOMNode* node, const char* attr_name, const std::string& attr_value,
    bool& flag, bool exception)
{
  flag = false;

  MemoryManager mm;
  int n(0);
  DOMElement* child = NULL;
	
  DOMNodeList* children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    DOMElement* element = static_cast<DOMElement*>(inode);
    if (element->hasAttribute(mm.transcode(attr_name))) {
      char* text = mm.transcode(element->getAttribute(mm.transcode(attr_name)));
      if (strcmp(text, attr_value.c_str()) == 0) {
        child = element;
        n++;
      }
    }
  }
  if (n == 1) flag = true;

  // exception
  if (!flag && exception) {
    Errors::Message msg;
    char* node_name = mm.transcode(node->getNodeName());
    msg << "Node \"" << node_name << "\" has no unique child with attribute \""
        << attr_name << "\" = \"" << attr_value << "\"\n";
    Exceptions::amanzi_throw(msg);
  }

  return child;
}


/* ******************************************************************
* Extracts children of the given node with the given name,
****************************************************************** */
std::vector<xercesc::DOMNode*> InputConverter::GetChildren_(
    xercesc::DOMNode* node, const std::string& name, bool& flag, bool exception)
{
  flag = false;

  MemoryManager mm;
  int n(0), m(0);
  std::vector<DOMNode*> namedChildren;

  DOMNodeList* children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    char* text = mm.transcode(inode->getNodeName());
    if (name == text) {
      namedChildren.push_back(inode);
      flag = true;
    }
  }

  // exception
  if (!flag and exception) {
    char* tagname = mm.transcode(node->getNodeName());
    Errors::Message msg;
    msg << "Amanzi::InputConverter: node \"" << tagname << "\" must have same elements\n";
    if (n) msg << "  The first element is \"" << name << "\".\n";
    msg << "  Please correct and try again.\n";
    Exceptions::amanzi_throw(msg);
  }

  return namedChildren;
}

xercesc::DOMElement* InputConverter::GetChildByName_(
    xercesc::DOMNode* node, const std::string& childName, bool& flag, bool exception)
{
  flag = false;

  MemoryManager mm;
  int n(0), m(0);
  DOMNode* child = NULL;

  DOMNodeList* children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    char* text = mm.transcode(inode->getNodeName());
    if (childName == text) {
      child = inode;
      flag = true;
      break;
    }
  }

  // exception
  if (!flag and exception) {
    Errors::Message msg;
    char* nodeName = mm.transcode(node->getNodeName());
    msg << "Amanzi::InputConverter: child node \"" << childName << "\" was not found for node \"" << nodeName << "\".\n";
    Exceptions::amanzi_throw(msg);
  }

  return static_cast<DOMElement*>(child);

}

/* ******************************************************************
* Extracts children and verifies that their have the common tagname.
* The name is the name of the first element.
****************************************************************** */
std::vector<DOMNode*> InputConverter::GetSameChildNodes_(
    DOMNode* node, std::string& name, bool& flag, bool exception)
{
  flag = false;

  MemoryManager mm;
  int n(0), m(0);
  std::vector<DOMNode*> same;

  DOMNodeList* children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    char* text = mm.transcode(inode->getNodeName());
    if (n == 0) name = text;
    if (strcmp(name.c_str(), text) == 0) {
      same.push_back(inode);
      n++;
    } 
    m++;
  }
  if (n == m) flag = true;

  // exception
  if (!flag && exception) {
    char* tagname = mm.transcode(node->getNodeName());
    Errors::Message msg;
    msg << "Node \"" << tagname << "\" must have similar elements.\n";
    if (n) msg << "The first element is \"" << name << "\".\n";
    msg << "Please correct and try again.\n";
    Exceptions::amanzi_throw(msg);
  }

  return same;
}


/* ******************************************************************
* Extract atribute of type double.
****************************************************************** */
double InputConverter::GetAttributeValueD_(
    DOMElement* elem, const char* attr_name, bool exception, double default_val)
{
  double val;
  MemoryManager mm;

  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    char* text = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    if (constants_.find(text) != constants_.end()) {  // check constants list
      val = TimeStringToValue_(constants_[text]);
    } else {
      val = TimeCharToValue_(text);
    }
  } else if (!exception) {
    val = default_val;
  } else {
    char* tagname = mm.transcode(elem->getNodeName());
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
  MemoryManager mm;

  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    char* text = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    if (constants_.find(text) != constants_.end()) {  // check constants list
      val = std::strtol(constants_[text].c_str(), NULL, 10);
    } else {
      val = std::strtol(text, NULL, 10);
    }
  } else if (! exception) {
    val = default_val;
  } else {
    char* tagname = mm.transcode(elem->getNodeName());
    ThrowErrorMissattr_(tagname, "attribute", attr_name, tagname);
  }

  return val;
}


/* ******************************************************************
* Extract atribute of type std::string.
****************************************************************** */
std::string InputConverter::GetAttributeValueS_(
    DOMElement* elem, const char* attr_name, bool exception, std::string default_val)
{
  std::string val;
  MemoryManager mm;

  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    val = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    boost::algorithm::trim(val);
    // check the list of global constants
    if (constants_.find(val) != constants_.end()) val = constants_[val];
  } else if (!exception) {
    val = default_val;
  } else {
    char* tagname = mm.transcode(elem->getNodeName());
    ThrowErrorMissattr_(tagname, "attribute", attr_name, tagname);
  }

  return val;
}


/* ******************************************************************
* Extract attribute of type vector<double>.
****************************************************************** */
std::vector<double> InputConverter::GetAttributeVector_(
    DOMElement* elem, const char* attr_name, bool exception)
{
  std::vector<double> val;
  MemoryManager mm;

  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    char* text_content = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    val = MakeCoordinates_(text_content);
  } else if (exception) {
    char* tagname = mm.transcode(elem->getNodeName());
    ThrowErrorMissattr_(tagname, "attribute", attr_name, tagname);
  }

  return val;
}

/* ******************************************************************
* Extract attribute of type vector<string>.
****************************************************************** */
std::vector<std::string> InputConverter::GetAttributeVectorS_(DOMElement* elem, const char* attr_name, bool exception)
{
  std::vector<std::string> val;
  MemoryManager mm;

  if (elem->hasAttribute(mm.transcode(attr_name))) {
    char* text_content = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    val = CharToStrings_(text_content);
  } else if (exception) {
    char* tagname = mm.transcode(elem->getNodeName());
    ThrowErrorMissattr_(tagname, "attribute", attr_name, tagname);
  }

  return val;
}

std::string InputConverter::GetChildValueS_(
    xercesc::DOMNode* node, const std::string& childName, bool& flag, bool exception)
{
  MemoryManager mm;

  std::string val;
  flag = false;

  DOMNodeList* children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) 
  {
    DOMNode* inode = children->item(i);
    char* tagname = mm.transcode(inode->getNodeName());   
    if (childName == tagname)
    {
      val = mm.transcode(inode->getTextContent());
      flag = true;
      break;
    }
  }

  if (!flag and exception)
  {
    char* nodeName = mm.transcode(node->getNodeName());
    ThrowErrorMisschild_(nodeName, childName, nodeName);
  }

  return val;
}

std::vector<std::string> InputConverter::GetChildVectorS_(
    xercesc::DOMNode* node, const std::string& childName, bool& flag, bool exception)
{
  MemoryManager mm;

  std::vector<std::string> val;
  flag = false;

  DOMNodeList* children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) 
  {
    DOMNode* inode = children->item(i);
    char* tagname = mm.transcode(inode->getNodeName());   
    if (childName == tagname)
    {
      val = CharToStrings_(mm.transcode(inode->getTextContent()));
      flag = true;
      break;
    }
  }

  if (!flag and exception)
  {
    char* nodeName = mm.transcode(node->getNodeName());
    ThrowErrorMisschild_(nodeName, childName, nodeName);
  }

  return val;
}

/* ******************************************************************
* Extract atribute of type std::string.
****************************************************************** */
std::string InputConverter::GetAttributeValueS_(
    DOMElement* elem, const char* attr_name, const char* options)
{
  std::string val;
  val = GetAttributeValueS_(elem, attr_name);

  std::vector<std::string> names = CharToStrings_(options);
  for (std::vector<std::string>::iterator it = names.begin(); it != names.end(); ++it) {
    if (val == *it) return val;
  }

  MemoryManager mm;
  char* tagname = mm.transcode(elem->getNodeName());
  Errors::Message msg;
  msg << "Validation of attribute \"" << attr_name << "\""
      << " for element \"" << tagname << "\" failed.\n";
  msg << "Available options: \"" << options << "\".\n";
  msg << "Please correct and try again.\n";
  Exceptions::amanzi_throw(msg);

  return val;
}


/* ******************************************************************
* Extract text content and verify it .
****************************************************************** */
std::string InputConverter::GetTextContentS_(
    DOMNode* node, const char* options, bool exception)
{
  std::string val;

  MemoryManager mm;
  val = TrimString_(mm.transcode(node->getTextContent()));

  std::vector<std::string> names = CharToStrings_(options);
  for (std::vector<std::string>::iterator it = names.begin(); it != names.end(); ++it) {
    if (val == *it) return val;
  }

  char* tagname = mm.transcode(node->getNodeName());
  Errors::Message msg;
  msg << "Validation of content \"" << val << "\" for node \"" << tagname << "\" failed.\n";
  msg << "Available options: \"" << options << "\".\n";
  msg << "Please correct and try again.\n";
  Exceptions::amanzi_throw(msg);
}


/* ******************************************************************
* Find positing in the array.
****************************************************************** */
int InputConverter::GetPosition_(const std::vector<std::string>& names, const std::string& name)
{
  for (int i = 0; i < names.size(); ++i) {
    if (strcmp(names[i].c_str(), name.c_str()) == 0) return i;
  }

  Errors::Message msg;
  msg << "Vector of names (e.g. solutes) has no \"" << name << "\".\n";
  msg << "Please correct and try again.\n";
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

  std::vector<std::string> regs;
  char* tmp2;
  tmp2 = strtok(tmp1, ",");
  if (tmp2 == NULL) {
    // No commas in the string means that it's a single name.
    regs.push_back(namelist);
  }

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
double InputConverter::TimeStringToValue_(const std::string& time_value)
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
double InputConverter::TimeCharToValue_(const char* time_value)
{
  double time;
  char* tmp1 = strcpy(new char[strlen(time_value) + 1], time_value);
  char* tmp2 = strtok(tmp1, ";, ");

  time = std::strtod(tmp2, NULL);
  tmp2 = strtok(NULL, ";,");

  if (tmp2 != NULL) {
    if (strcmp(tmp2, "y") == 0) { 
      time *= 365.25 * 24.0 * 3600.0;
    } else if (strcmp(tmp2, "m") == 0) {
      time *= 365.25 * 2.0 * 3600.0;
    } else if (strcmp(tmp2, "d") == 0) {
      time *= 24.0 * 3600.0;
    } else if (strcmp(tmp2, "h") == 0) {
      time *= 3600.0;
    }
  }
  
  delete[] tmp1;
  return time;
}


/* ******************************************************************
* Converts coordinate string to an array of doubles.
****************************************************************** */
std::vector<double> InputConverter::MakeCoordinates_(const std::string& array)
{
  std::vector<double> coords;
  char* tmp1 = strcpy(new char[array.size() + 1], array.c_str());
  char* tmp2 = strtok(tmp1, "(, ");

  while (tmp2 != NULL) {
    std::string str(tmp2);
    boost::algorithm::trim(str);
    coords.push_back(std::strtod(str.c_str(), NULL));
    tmp2 = strtok(NULL, ",");
  }

  delete[] tmp1;
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
* Generate error message when list is empty.
******************************************************************* */
int InputConverter::IsEmpty(DOMNodeList* node_list, const std::string& name, bool exception)
{
  int n = node_list->getLength();
  if (n == 0 && exception) {
    Errors::Message msg;
    msg << "Mandatory element \"" << name << "\" is missing.\n";
    msg << "Please correct and try again.\n";
    Exceptions::amanzi_throw(msg);
  }

  return n;
}


/* *******************************************************************
* Generate unified error message for ill-formed element
******************************************************************* */
void InputConverter::ThrowErrorIllformed_(
    const std::string& section, const std::string& type, const std::string& ill_formed)
{
  Errors::Message msg;
  msg << "An error occurred during parsing node \"" << section << "\"\n";
  msg << "Missing or ill-formed \"" << type << "\" for \"" << ill_formed << "\".\n";
  msg << "Please correct and try again.\n";
  Exceptions::amanzi_throw(msg);
}


/* *****************************************************************************
* Generate unified error message for ill-formed element with options provided
***************************************************************************** */
void InputConverter::ThrowErrorIllformed_(
    const std::string& section, const std::string& type, const std::string& ill_formed, const std::string& options)
{
  Errors::Message msg;
  msg << "An error occurred during parsing node \"" << section << "\"\n";
  msg << "Missing or ill-formed " << type << " for \"" << ill_formed << "\"\n";
  msg << "Valid options are: " << options << "\n";
  msg << "Please correct and try again.\n" ;
  Exceptions::amanzi_throw(msg);
}


/* *******************************************************************
* Generate unified error message for missing item
******************************************************************* */
void InputConverter::ThrowErrorMissattr_(
    const std::string& section, const std::string& type, const std::string& missing, const std::string& name)
{
  Errors::Message msg;
  msg << "An error occurred during parsing node \"" << section << "\"\n";
  msg << "No " << type << " \"" << missing << "\" found for \"" << name << "\".\n";
  msg << "Please correct and try again \n";
  Exceptions::amanzi_throw(msg);
}

/* *******************************************************************
* Generate unified error message for missing child
******************************************************************* */
void InputConverter::ThrowErrorMisschild_(
    const std::string& section, const std::string& missing, const std::string& name)
{
  Errors::Message msg;
  msg << "Amanzi::InputConverter: an error occurred during parsing node \"" << section << "\"\n";
  msg << "  No child \"" << missing << "\" found for \"" << name << "\".\n";
  msg << "  Please correct and try again \n";
  Exceptions::amanzi_throw(msg);
}

}  // namespace AmanziInput
}  // namespace Amanzi
