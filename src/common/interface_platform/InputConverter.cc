/*
  Input Converter

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
#include "xercesc/parsers/XercesDOMParser.hpp"
#include "xercesc/parsers/DOMLSParserImpl.hpp"
#include "xercesc/framework/StdOutFormatTarget.hpp"
#include "xercesc/util/OutOfMemoryException.hpp"

// Amanzi's
#include "ErrorHandler.hpp"
#include "InputConverter.hh"

namespace Amanzi {
namespace AmanziInput {

/* ******************************************************************
* Non-member function: returns parser.
****************************************************************** */
XercesDOMParser* CreateXMLParser()
{
  XMLPlatformUtils::Initialize();

  // Set up an XML DOM parser.
  XercesDOMParser* parser = new XercesDOMParser();
  parser->setExitOnFirstFatalError(true);
  parser->setValidationConstraintFatal(true);
  parser->setValidationScheme(XercesDOMParser::Val_Never);
  parser->setDoNamespaces(true);
  parser->setCreateCommentNodes(false);
  parser->useCachedGrammarInParse(true);

  return parser;
}

/* ******************************************************************
* Non-member function: returns xercesc document.
****************************************************************** */
xercesc::DOMDocument* OpenXMLInput(XercesDOMParser* parser,
                                   const std::string& xml_input)
{
  // Open, parse, and return the document.
  AmanziErrorHandler* errorHandler = new AmanziErrorHandler();
  parser->setErrorHandler(errorHandler);

  try {
    parser->parse(xml_input.c_str());
  }
  catch (const OutOfMemoryException& e) {
    std::cerr << "OutOfMemoryException" << std::endl;
    Exceptions::amanzi_throw(Errors::Message("Ran out of memory while parsing the input file. Aborting."));
  }
  catch (...) {
    Exceptions::amanzi_throw(Errors::Message("Errors occured while parsing the input file. Aborting."));
  }

  parser->setErrorHandler(NULL);
  delete errorHandler;
  xercesc::DOMDocument* doc = parser->getDocument();
  return doc;
}


/* ******************************************************************
* Various constructors.
****************************************************************** */
InputConverter::InputConverter(const std::string& input_filename):
    xmlfilename_(input_filename),
    doc_(NULL),
    parser_(NULL)
{
  parser_ = CreateXMLParser();
  doc_ = OpenXMLInput(parser_, input_filename);
  FilterNodes(doc_, "comments");
}

InputConverter::InputConverter(const std::string& input_filename,
                               xercesc::DOMDocument* input_doc):
    xmlfilename_(input_filename),
    doc_(input_doc),
    parser_(NULL)
{
  FilterNodes(doc_, "comments");
}

InputConverter::~InputConverter()
{
  // if (doc_ != NULL)
  //   delete doc_;
  if (parser_ != NULL)
    delete parser_;
}

XERCES_CPP_NAMESPACE_USE

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
* Check the version number.
****************************************************************** */
void InputConverter::ParseVersion_()
{
  MemoryManager mm;
  
  DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode("amanzi_input"));
  if (node_list->getLength() > 0) {
    std::string version = GetAttributeValueS_(static_cast<DOMElement*>(node_list->item(0)), "version");
    
    int major, minor, micro;
    
    std::stringstream ss;
    ss << version;
    std::string ver;
    
    try {
      getline(ss, ver, '.');
      major = boost::lexical_cast<int>(ver);
      
      getline(ss, ver, '.');
      minor = boost::lexical_cast<int>(ver);
      
      getline(ss,ver);
      micro = boost::lexical_cast<int>(ver);
    }
    catch (...) {
      Errors::Message msg("The version string in the input file '" + version + 
                          "' has the wrong format, use I.J.K.");
      Exceptions::amanzi_throw(msg);
    }

    if ((major != AMANZI_SPEC_VERSION_MAJOR) ||
        (minor != AMANZI_SPEC_VERSION_MINOR) || 
        (micro != AMANZI_SPEC_VERSION_MICRO)) {
      std::stringstream ss;
      ss << AMANZI_SPEC_VERSION_MAJOR << "." << AMANZI_SPEC_VERSION_MINOR << "." << AMANZI_SPEC_VERSION_MICRO;

      Errors::Message msg;
      msg << "The input version " << version << " is not supported. "
          << "Supported versions: "<< ss.str() << ".\n";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    // amanzi input description did not exist, error
    Errors::Message msg("Amanzi input description does not exist <amanzi_input version=...>");
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Populates protected std::map constants_.
****************************************************************** */
void InputConverter::ParseConstants_()
{
  MemoryManager mm;

  char *tagname, *text;
  DOMNode* node;
  DOMNodeList *node_list, *children;
  DOMElement* element;

  // process constants: we ignore type of generic constants.
  node_list = doc_->getElementsByTagName(mm.transcode("constants"));
  if (node_list->getLength() == 0) return;

  children = node_list->item(0)->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    element = static_cast<DOMElement*>(inode);
    
    std::string name, value, type(TYPE_NOT_CONSTANT);
    char* tagname = mm.transcode(inode->getNodeName());   
    if (strcmp(tagname, "constant") == 0) {
      type = GetAttributeValueS_(element, "type");
    } else if (strcmp(tagname, "time_constant") == 0) {
      type = TYPE_TIME;
    } else if (strcmp(tagname, "numerical_constant") == 0) {
      type = TYPE_NUMERICAL;
    } else if (strcmp(tagname, "area_mass_flux_constant") == 0) {
      type = TYPE_AREA_MASS_FLUX;
    }

    if (type != TYPE_NOT_CONSTANT) {
      name = GetAttributeValueS_(element, "name");
      value = GetAttributeValueS_(element, "value");

      if (type == TYPE_TIME) {
        constants_time_[name] = value;
      } else if (type == TYPE_NUMERICAL) {
        constants_numerical_[name] = value;
      } else if (type == TYPE_AREA_MASS_FLUX) {
        constants_area_mass_flux_[name] = value;
      } else if (type == TYPE_NONE) {
        constants_[name] = value;
      }
    }
  }
}


/* ******************************************************************
* Returns node specified by the list of consequtive names tags 
* separated by commas. Only the first tag may be not unique.
****************************************************************** */
DOMNode* InputConverter::GetUniqueElementByTagsString_(
    const std::string& tags, bool& flag)
{
  flag = false;

  MemoryManager mm;
  DOMNode* node = NULL;
  DOMNode* node_good;

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
  int nnodes = node_list->getLength();
  if (nnodes == 0) return node;

  int icnt(0);
  for (int k = 0; k < nnodes; ++k) {
    node = node_list->item(k);

    bool found(true);
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
      if (ntag != 1) {
        found = false;
        break;
      }
    }

    if (found) {
      icnt++;
      node_good = node;
    }
  }

  flag = (icnt == 1);
  return node_good;
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
* Extract attribute of type double.
****************************************************************** */
double InputConverter::GetAttributeValueD_(
    DOMElement* elem, const char* attr_name,
    const std::string& type, bool exception, double default_val)
{
  double val;
  MemoryManager mm;

  Errors::Message msg;
  std::string text, parsed, found_type;

  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    text = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    if (text.size() == 0) {
      msg << "Attribute \"" << attr_name << "\" cannot be empty.\n";
      Exceptions::amanzi_throw(msg);
    }

    // process constants
    found_type = GetConstantType_(text, parsed);
    val = TimeStringToValue_(parsed);

    // no checks for two types
    if (found_type == TYPE_NONE ||
        found_type == TYPE_NOT_CONSTANT) found_type = type;

    if (type != found_type) {
      msg << "Usage of constant \"" << text << "\" of type=" << found_type 
          << ". Expect type=" << type << ".\n";
      Exceptions::amanzi_throw(msg);
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
    DOMElement* elem, const char* attr_name,
    const std::string& type, bool exception, int default_val)
{
  int val;
  MemoryManager mm;

  std::string text, parsed, found_type;
  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    text = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));

    // process constants 
    found_type = GetConstantType_(text, parsed);
    val = std::strtol(parsed.c_str(), NULL, 10);

    // no checks for two types
    if (found_type == TYPE_NONE ||
        found_type == TYPE_NOT_CONSTANT) found_type = type;

    if (type != found_type) {
      Errors::Message msg;
      msg << "Usage of constant \"" << text << "\" of type=" << found_type 
          << ". Expect type=" << type << ".\n";
      Exceptions::amanzi_throw(msg);
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
    DOMElement* elem, const char* attr_name,
    const std::string& type, bool exception, std::string default_val)
{
  std::string val;
  MemoryManager mm;

  std::string text, found_type;
  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    text = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    boost::algorithm::trim(text);

    // check the list of global constants
    found_type = GetConstantType_(text, val);

    // no checks for two types
    if (found_type == TYPE_NONE ||
        found_type == TYPE_NOT_CONSTANT) found_type = type;

    if (type != found_type) {
      Errors::Message msg;
      msg << "Usage of constant \"" << text << "\" of type=" << found_type 
          << ". Expect type=" << type << ".\n";
      Exceptions::amanzi_throw(msg);
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
  std::string text = TrimString_(mm.transcode(node->getTextContent()));
  GetConstantType_(text, val);

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
  return "";
}


/* ******************************************************************
* Find positing in the array.
****************************************************************** */
std::string InputConverter::GetConstantType_(
    const std::string& val, std::string& parsed_val)
{
  std::string type;
  if (constants_time_.find(val) != constants_time_.end()) {
    type = TYPE_TIME;
    parsed_val = constants_time_[val];
  } else if (constants_numerical_.find(val) != constants_numerical_.end()) {
    type = TYPE_NUMERICAL;
    parsed_val = constants_numerical_[val];
  } else if (constants_area_mass_flux_.find(val) != constants_area_mass_flux_.end()) {
    type = TYPE_AREA_MASS_FLUX;
    parsed_val = constants_area_mass_flux_[val];
  } else if (constants_.find(val) != constants_.end()) {
    type = TYPE_NONE;
    parsed_val = constants_[val];
  } else {
    type = TYPE_NOT_CONSTANT;
    parsed_val = val;
  }

  return type;
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
    std::string str(tmp2), parsed_str;
    boost::algorithm::trim(str);

    GetConstantType_(str, parsed_str);
    coords.push_back(std::strtod(parsed_str.c_str(), NULL));
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
  msg << "  No child \"" << missing << "\" found";
  if (!name.empty()) {
    msg << " for \"" << name << "\"";
  }
  msg << ".\n";
  msg << "  Please correct and try again \n";
  Exceptions::amanzi_throw(msg);
}

}  // namespace AmanziInput
}  // namespace Amanzi
