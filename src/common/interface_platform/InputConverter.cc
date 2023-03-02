/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Erin Barker (Erin.Barker@pnnl.gov)
*/

/*
  Input Converter

*/

#include <algorithm>
#include <cfloat>
#include <fstream>
#include <locale>
#include <cfloat>
#include <sstream>
#include <string>
#include <sys/stat.h>

// TPLs
#include <boost/filesystem/operations.hpp>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#define BOOST_FILESYTEM_NO_DEPRECATED

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/XMLString.hpp"
#include "xercesc/util/PlatformUtils.hpp"
#include "xercesc/parsers/XercesDOMParser.hpp"
#include "xercesc/parsers/DOMLSParserImpl.hpp"
#include "xercesc/framework/StdOutFormatTarget.hpp"
#include "xercesc/util/OutOfMemoryException.hpp"

// Amanzi's
#include "ErrorHandler.hpp"
#include "StringExt.hh"
#include "InputConverter.hh"

namespace Amanzi {
namespace AmanziInput {


/* ******************************************************************
* Non-member function: returns parser.
****************************************************************** */
XercesDOMParser*
CreateXMLParser()
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
xercesc::DOMDocument*
OpenXMLInput(XercesDOMParser* parser, const std::string& xml_input)
{
  // Open, parse, and return the document.
  AmanziErrorHandler* errorHandler = new AmanziErrorHandler();
  parser->setErrorHandler(errorHandler);

  try {
    parser->parse(xml_input.c_str());
  } catch (const OutOfMemoryException& e) {
    std::cerr << "OutOfMemoryException" << std::endl;
    Exceptions::amanzi_throw(
      Errors::Message("Ran out of memory while parsing the input file. Aborting."));
  } catch (...) {
    Exceptions::amanzi_throw(
      Errors::Message("Errors occured while parsing the input file. Aborting."));
  }

  parser->setErrorHandler(NULL);
  delete errorHandler;
  xercesc::DOMDocument* doc = parser->getDocument();
  return doc;
}


/* ******************************************************************
* Various constructors.
****************************************************************** */
InputConverter::InputConverter(const std::string& input_filename)
  : units_("molar"), xmlfilename_(input_filename), doc_(NULL), parser_(NULL)
{
  parser_ = CreateXMLParser();
  doc_ = OpenXMLInput(parser_, input_filename);
  FilterNodes("comments");
}

InputConverter::InputConverter(const std::string& input_filename, xercesc::DOMDocument* input_doc)
  : units_("molar"), xmlfilename_(input_filename), doc_(input_doc), parser_(NULL)
{
  FilterNodes("comments");
}

InputConverter::~InputConverter()
{
  // if (doc_ != NULL)
  //   delete doc_;
  if (parser_ != NULL) delete parser_;
}

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Filters out all nodes named "filter" starting with node "parent".
****************************************************************** */
void
InputConverter::FilterNodes(const std::string& filter)
{
  MemoryManager mm;

  DOMNodeList* remove = doc_->getElementsByTagName(mm.transcode(filter.c_str()));
  while (remove->getLength() > 0) {
    DOMNode* child = remove->item(0);

    if (child->getNodeType() == DOMNode::ELEMENT_NODE) {
      child->getParentNode()->removeChild(child);
    }
  }
}


/* ******************************************************************
* Check the version number.
****************************************************************** */
void
InputConverter::ParseVersion_()
{
  MemoryManager mm;

  DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode("amanzi_input"));
  if (node_list->getLength() > 0) {
    std::string version = GetAttributeValueS_(node_list->item(0), "version");

    int major, minor, micro;

    std::stringstream ss;
    ss << version;
    std::string ver;

    try {
      getline(ss, ver, '.');
      major = std::stoi(ver);

      getline(ss, ver, '.');
      minor = std::stoi(ver);

      getline(ss, ver);
      micro = std::stoi(ver);
    } catch (...) {
      Errors::Message msg("The version string in the input file '" + version +
                          "' has the wrong format, use I.J.K.");
      Exceptions::amanzi_throw(msg);
    }

    if ((major != AMANZI_SPEC_VERSION_MAJOR) || (minor != AMANZI_SPEC_VERSION_MINOR) ||
        (micro < AMANZI_SPEC_VERSION_MICRO)) {
      std::stringstream ss1;
      ss1 << AMANZI_SPEC_VERSION_MAJOR << "." << AMANZI_SPEC_VERSION_MINOR << "."
          << AMANZI_SPEC_VERSION_MICRO;

      Errors::Message msg;
      msg << "The input version " << version << " is not supported. "
          << "Supported versions: " << ss1.str() << " and higher.\n";
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
void
InputConverter::ParseConstants_()
{
  MemoryManager mm;

  DOMNodeList *node_list, *children;
  DOMElement* element;

  // process constants: we ignore type of generic constants.
  node_list = doc_->getElementsByTagName(mm.transcode("constants"));
  if (node_list->getLength() == 0) return;

  children = node_list->item(0)->getChildNodes();
  int nchildren = children->getLength();

  for (int i = 0; i < nchildren; ++i) {
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
* Verify miscallenous assumptions for specific engines.
****************************************************************** */
void
InputConverter::ParseGeochemistry_()
{
  bool flag;
  DOMNode* node = GetUniqueElementByTagsString_("process_kernels, chemistry", flag);
  if (!flag) return;

  std::string engine = GetAttributeValueS_(node, "engine");

  node = GetUniqueElementByTagsString_("geochemistry, constraints", flag);
  if (flag && engine.substr(0, 10) == "crunchflow") {
    std::vector<DOMNode*> children = GetChildren_(node, "constraint", flag);

    for (int i = 0; i < children.size(); ++i) {
      std::string name = GetAttributeValueS_(children[i], "name");
      std::string tmp(name);
      transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
      if (tmp != name) {
        Errors::Message msg;
        msg << "CrunchFlow uses lower-case for geochemical constraint name, \"" << name
            << "\" does not comply.\n";
        Exceptions::amanzi_throw(msg);
      }
    }
  }
}


/* ******************************************************************
* Returns node specified by the list of consequitive namestags
* separated by commas. Only the first tag may be not unique.
****************************************************************** */
DOMNode*
InputConverter::GetUniqueElementByTagsString_(const std::string& tags, bool& flag)
{
  flag = false;

  MemoryManager mm;
  DOMNode* node = NULL;
  DOMNode* node_good;

  std::vector<std::string> tag_names = CharToStrings_(tags.c_str());
  if (tag_names.size() == 0) return node;

  if (tag_names.size() == 1) {
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
      int ntag(0);
      DOMNode* inode = node->getFirstChild();

      while (inode != NULL) {
        if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
          char* tagname = mm.transcode(inode->getNodeName());
          if (strcmp(tagname, tag_names[n].c_str()) == 0) {
            node = inode;
            ntag++;
          }
        }
        inode = inode->getNextSibling();
      }

      if (ntag > 1) {
        Errors::Message msg;
        msg << "Tag \"" << tag_names[n] << "\" in \"" << tags << "\"\n";
        msg << "  must be either unique or absent (default value will be used).\n";
        msg << "  Please correct and try again.\n";
        Exceptions::amanzi_throw(msg);
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
* Returns node specified by the list of consequtive names tags
* separated by commas. Only the first tag may be not unique.
****************************************************************** */
xercesc::DOMNode*
InputConverter::GetUniqueElementByTagsString_(const std::string& tags, bool& flag, bool exception)
{
  xercesc::DOMNode* node = GetUniqueElementByTagsString_(tags, flag);

  if (!flag && exception) {
    Errors::Message msg;
    msg << "No unique element for tags \"" << tags << "\"\n";
    Exceptions::amanzi_throw(msg);
  }

  return node;
}


/* ******************************************************************
* Return node described by the list of consequtive names tags
* separated by commas.
****************************************************************** */
DOMNode*
InputConverter::GetUniqueElementByTagsString_(const DOMNode* node1,
                                              const std::string& tags,
                                              bool& flag)
{
  DOMNode* node = nullptr;

  flag = false;
  std::vector<std::string> tag_names = CharToStrings_(tags.c_str());
  if (tag_names.size() == 0) return node;

  // get the first node
  node = const_cast<DOMNode*>(node1);

  for (int n = 0; n < tag_names.size(); ++n) {
    int ntag(0);
    DOMNode* inode = node->getFirstChild();

    while (inode != NULL) {
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
        char* tagname = XMLString::transcode(inode->getNodeName());
        if (strcmp(tagname, tag_names[n].c_str()) == 0) {
          node = inode;
          ntag++;
        }
        XMLString::release(&tagname);
      }
      inode = inode->getNextSibling();
    }
    if (ntag != 1) return node;
  }

  flag = true;
  return node;
}


/* ******************************************************************
* Returns the child with the given attribute name and value.
****************************************************************** */
DOMElement*
InputConverter::GetUniqueChildByAttribute_(xercesc::DOMNode* node,
                                           const char* attr_name,
                                           const std::string& attr_value,
                                           bool& flag,
                                           bool exception)
{
  flag = false;

  MemoryManager mm;
  int n(0);
  DOMElement* child = NULL;

  DOMNodeList* children = node->getChildNodes();
  int nchildren = children->getLength();

  for (int i = 0; i < nchildren; ++i) {
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
    char* parname = mm.transcode(node->getParentNode()->getNodeName());
    msg << "Node \"" << parname << "->" << node_name << "\" has no unique child with attribute \""
        << attr_name << "\" = \"" << attr_value << "\"\n";
    Exceptions::amanzi_throw(msg);
  }

  return child;
}


/* ******************************************************************
* Extracts children of the given node with the given name,
****************************************************************** */
std::vector<xercesc::DOMNode*>
InputConverter::GetChildren_(DOMNode* node, const std::string& name, bool& flag, bool exception)
{
  flag = false;

  MemoryManager mm;
  int n(0);
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
  if (!flag && exception) {
    char* tagname = mm.transcode(node->getNodeName());
    char* parname = mm.transcode(node->getParentNode()->getNodeName());
    Errors::Message msg;
    msg << "Node \"" << parname << "->" << tagname << "\" must have same elements\n";
    if (n) msg << "  The first element is \"" << name << "\".\n";
    msg << "  Please correct and try again.\n";
    Exceptions::amanzi_throw(msg);
  }

  return namedChildren;
}


/* ******************************************************************
* Returns any child with the given name.
****************************************************************** */
DOMElement*
InputConverter::GetChildByName_(xercesc::DOMNode* node,
                                const std::string& childName,
                                bool& flag,
                                bool exception)
{
  flag = false;

  MemoryManager mm;
  DOMNode* child = NULL;

  DOMNodeList* children = node->getChildNodes();
  int nchildren = children->getLength();

  for (int i = 0; i < nchildren; ++i) {
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
  if (!flag && exception) {
    Errors::Message msg;
    char* node_name = mm.transcode(node->getNodeName());
    char* parname = mm.transcode(node->getParentNode()->getNodeName());
    msg << "Child node \"" << childName << "\" was not found for node \"" << parname << "->"
        << node_name << "\".\n";
    Exceptions::amanzi_throw(msg);
  }

  return static_cast<DOMElement*>(child);
}


/* ******************************************************************
* Extracts children and verifies that their have the common tagname.
* Returned name is the tagname of the first element.
****************************************************************** */
std::vector<DOMNode*>
InputConverter::GetSameChildNodes_(DOMNode* node, std::string& name, bool& flag, bool exception)
{
  flag = false;

  MemoryManager mm;
  int n(0), m(0);
  std::vector<DOMNode*> same;

  name = "";
  DOMNode* inode = node->getFirstChild();
  while (inode != NULL) {
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) {
      inode = inode->getNextSibling();
      continue;
    }

    char* text = mm.transcode(inode->getNodeName());
    if (n == 0) name = text;
    if (strcmp(name.c_str(), text) == 0) {
      same.push_back(inode);
      n++;
    }
    m++;
    inode = inode->getNextSibling();
  }
  if (n == m && n > 0) flag = true;

  // exception
  if (!flag && exception) {
    char* tagname = mm.transcode(node->getNodeName());
    char* parname = mm.transcode(node->getParentNode()->getNodeName());
    Errors::Message msg;
    msg << "Node \"" << parname << "->" << tagname << "\" must have similar elements.\n";
    if (n) msg << "The first element is \"" << name << "\".\n";
    msg << "Please correct and try again.\n";
    Exceptions::amanzi_throw(msg);
  }

  return same;
}


/* ******************************************************************
* Extract attribute of type double.
****************************************************************** */
double
InputConverter::GetAttributeValueD_(DOMElement* elem,
                                    const char* attr_name,
                                    const std::string& type,
                                    double valmin,
                                    double valmax,
                                    std::string unit,
                                    bool exception,
                                    double default_val)
{
  double val;
  MemoryManager mm;

  Errors::Message msg;
  std::string text, parsed, found_type, unit_in;

  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    text = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    if (text.size() == 0) {
      msg << "Attribute \"" << attr_name << "\" cannot be empty.\n";
      Exceptions::amanzi_throw(msg);
    }

    // process constants and known units
    // -- extract units
    found_type = GetConstantType_(text, parsed);
    val = ConvertUnits_(parsed, unit_in);

    // -- replace empty input unit by the default unit
    if (unit_in == "") {
      bool flag;
      unit_in = units_.ConvertUnitS(unit, units_.system());
      val = units_.ConvertUnitD(val, unit_in, "SI", -1.0, flag);
    }

    // no checks for two types
    if (found_type == TYPE_NONE || found_type == TYPE_NOT_CONSTANT) found_type = type;

    if (type != found_type) {
      msg << "Usage of constant \"" << text << "\" of type=" << found_type
          << ". Expect type=" << type << ".\n";
      Exceptions::amanzi_throw(msg);
    }

    if (!(val >= valmin && val <= valmax)) {
      msg << "Value of attribute \"" << attr_name << "\"=" << val
          << "\" is out of range: " << valmin << " " << valmax << " [" << unit << "].\n";
      Exceptions::amanzi_throw(msg);
    }

    if ((unit != "" && unit_in != "") || (unit == "-" && unit_in != "")) {
      if (!units_.CompareUnits(unit, unit_in)) {
        msg << "Input unit [" << unit_in << "] for attribute \"" << attr_name
            << "\" does not match the expected unit [" << unit << "].\n";
        Exceptions::amanzi_throw(msg);
      }
    }
  } else if (!exception) {
    val = default_val;
  } else {
    char* tagname = mm.transcode(elem->getNodeName());
    ThrowErrorMissing_(tagname, "attribute", attr_name, tagname);
    val = default_val;
  }

  return val;
}


/* ******************************************************************
* Extract atribute of type int.
****************************************************************** */
int
InputConverter::GetAttributeValueL_(DOMElement* elem,
                                    const char* attr_name,
                                    const std::string& type,
                                    int valmin,
                                    int valmax,
                                    bool exception,
                                    int default_val)
{
  int val(INT_MIN);
  MemoryManager mm;

  std::string text, parsed, found_type;
  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    text = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));

    // process constants
    found_type = GetConstantType_(text, parsed);
    val = std::strtol(parsed.c_str(), NULL, 10);

    // no checks for two types
    if (found_type == TYPE_NONE || found_type == TYPE_NOT_CONSTANT) found_type = type;

    if (type != found_type) {
      Errors::Message msg;
      msg << "Usage of constant \"" << text << "\" of type=" << found_type
          << ". Expect type=" << type << ".\n";
      Exceptions::amanzi_throw(msg);
    }

    if (!(val >= valmin && val <= valmax)) {
      Errors::Message msg;
      msg << "Value of attribute \"" << attr_name << "\"=" << val << " is out of range: " << valmin
          << " " << valmax << ".\n";
      Exceptions::amanzi_throw(msg);
    }
  } else if (!exception) {
    val = default_val;
  } else {
    char* tagname = mm.transcode(elem->getNodeName());
    ThrowErrorMissing_(tagname, "attribute", attr_name, tagname);
  }

  return val;
}


/* ******************************************************************
* Extract atribute of type std::string.
****************************************************************** */
std::string
InputConverter::GetAttributeValueS_(DOMElement* elem,
                                    const char* attr_name,
                                    const std::string& type,
                                    bool exception,
                                    std::string default_val)
{
  std::string val;
  MemoryManager mm;

  std::string text, found_type;
  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    text = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    trim(text);

    // check the list of global constants
    found_type = GetConstantType_(text, val);

    // no checks for two types
    if (found_type == TYPE_NONE || found_type == TYPE_NOT_CONSTANT) found_type = type;

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
    ThrowErrorMissing_(tagname, "attribute", attr_name, tagname);
  }

  return val;
}


/* ******************************************************************
* Extract attribute of type vector<double>. Unit of each component
* must match the expected unit.
****************************************************************** */
std::vector<double>
InputConverter::GetAttributeVectorD_(DOMElement* elem,
                                     const char* attr_name,
                                     int length,
                                     std::string unit,
                                     bool exception,
                                     double mol_mass)
{
  std::vector<double> val;
  MemoryManager mm;

  if (elem != NULL && elem->hasAttribute(mm.transcode(attr_name))) {
    std::vector<std::string> unit_in;
    char* text_content = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    val = MakeVector_(text_content, unit_in, mol_mass);

    for (int i = 0; i < unit_in.size(); ++i) {
      if ((unit != "" && unit_in[i] != "") || (unit == "-" && unit_in[i] != "")) {
        if (!units_.CompareUnits(unit, unit_in[i])) {
          Errors::Message msg;
          msg << "Input unit [" << unit_in[i] << "] for attribute \"" << attr_name
              << "\" does not match the expected unit [" << unit << "].\n";
          Exceptions::amanzi_throw(msg);
        }
      }
    }

    if (length > 0 && val.size() != length) {
      Errors::Message msg;
      msg << "Attribute \"" << attr_name << "\" has too few parameters: " << (int)val.size()
          << ", expected: " << length << ". Hint: check \"mesh->dimension\".\n";
      Exceptions::amanzi_throw(msg);
    }

  } else if (exception) {
    char* tagname = mm.transcode(elem->getNodeName());
    ThrowErrorMissing_(tagname, "attribute", attr_name, tagname);
  }

  return val;
}


/* ******************************************************************
* Extract attribute of type vector<string>.
****************************************************************** */
std::vector<std::string>
InputConverter::GetAttributeVectorS_(DOMElement* elem, const char* attr_name, bool exception)
{
  std::vector<std::string> val, tmp;
  MemoryManager mm;

  if (elem->hasAttribute(mm.transcode(attr_name))) {
    char* text_content = mm.transcode(elem->getAttribute(mm.transcode(attr_name)));
    tmp = CharToStrings_(text_content);
    for (int i = 0; i < tmp.size(); ++i) {
      std::string parsed_val;
      GetConstantType_(tmp[i], parsed_val);
      val.push_back(parsed_val);
    }
  } else if (exception) {
    char* tagname = mm.transcode(elem->getNodeName());
    ThrowErrorMissing_(tagname, "attribute", attr_name, tagname);
  }

  return val;
}


/* ******************************************************************
* Returns string value of <node> -> <childName>value</childName>.
****************************************************************** */
std::string
InputConverter::GetChildValueS_(DOMNode* node,
                                const std::string& childName,
                                bool& flag,
                                bool exception)
{
  MemoryManager mm;

  std::string val;
  flag = false;

  DOMNodeList* children = node->getChildNodes();
  int nchildren = children->getLength();
  for (int i = 0; i < nchildren; ++i) {
    DOMNode* inode = children->item(i);
    char* tagname = mm.transcode(inode->getNodeName());
    if (childName == tagname) {
      val = mm.transcode(inode->getTextContent());
      flag = true;
      break;
    }
  }

  if (!flag and exception) {
    char* nodeName = mm.transcode(node->getNodeName());
    ThrowErrorMisschild_(nodeName, childName, nodeName);
  }

  return val;
}


/* ******************************************************************
* Returns array of values in <node> -> <childName>val1,val2,val3</childName>
****************************************************************** */
std::vector<std::string>
InputConverter::GetChildVectorS_(DOMNode* node,
                                 const std::string& childName,
                                 bool& flag,
                                 bool exception)
{
  MemoryManager mm;

  std::vector<std::string> val;
  flag = false;

  DOMNodeList* children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    char* tagname = mm.transcode(inode->getNodeName());
    if (childName == tagname) {
      val = CharToStrings_(mm.transcode(inode->getTextContent()));
      flag = true;
      break;
    }
  }

  if (!flag and exception) {
    char* nodeName = mm.transcode(node->getNodeName());
    ThrowErrorMisschild_(nodeName, childName, nodeName);
  }

  return val;
}


/* ******************************************************************
* Extract atribute of type std::string.
****************************************************************** */
std::string
InputConverter::GetAttributeValueS_(DOMNode* node, const char* attr_name, const char* options)
{
  DOMElement* element = static_cast<DOMElement*>(node);

  std::string val;
  val = GetAttributeValueS_(element, attr_name);

  std::vector<std::string> names = CharToStrings_(options);
  for (std::vector<std::string>::iterator it = names.begin(); it != names.end(); ++it) {
    if (val == *it) return val;
  }

  MemoryManager mm;
  char* tagname = mm.transcode(element->getNodeName());
  Errors::Message msg;
  msg << "Validation of attribute \"" << attr_name << "\""
      << " for element \"" << tagname << "\" failed.\n";
  msg << "Available options: \"" << options << "\".\n";
  msg << "Please correct and try again.\n";
  Exceptions::amanzi_throw(msg);

  return val;
}


/* ******************************************************************
* Extract text content and convert it for double
****************************************************************** */
double
InputConverter::GetTextContentD_(DOMNode* node,
                                 std::string unit,
                                 bool exception,
                                 double default_val)
{
  double val(DBL_MIN);
  std::string parsed_text, unit_in;

  MemoryManager mm;

  if (node != NULL) {
    std::string text = TrimString_(mm.transcode(node->getTextContent()));
    GetConstantType_(text, parsed_text);
    val = ConvertUnits_(parsed_text, unit_in);

    // replace empty input unit by the default unit
    if (unit_in == "") {
      bool flag;
      unit_in = units_.ConvertUnitS(unit, units_.system());
      val = units_.ConvertUnitD(val, unit_in, "SI", -1.0, flag);
    }

    if ((unit != "" && unit_in != "") || (unit == "-" && unit_in != "")) {
      if (!units_.CompareUnits(unit, unit_in)) {
        char* tagname = mm.transcode(node->getNodeName());
        Errors::Message msg;
        msg << "Input unit [" << unit_in << "] for element \"" << tagname
            << "\" does not match the expected unit [" << unit << "].\n";
        Exceptions::amanzi_throw(msg);
        val = 0.0;
      }
    }
  } else if (!exception) {
    val = default_val;
  } else {
    char* tagname = mm.transcode(node->getNodeName());
    ThrowErrorMisschild_("unknown", tagname, "unknwon");
  }

  return val;
}


/* ******************************************************************
* Extract text content and verify it against list of options.
****************************************************************** */
std::string
InputConverter::GetTextContentS_(DOMNode* node, const char* options, bool exception)
{
  std::string val;

  MemoryManager mm;
  std::string text = TrimString_(mm.transcode(node->getTextContent()));
  GetConstantType_(text, val);

  std::vector<std::string> names = CharToStrings_(options);
  for (std::vector<std::string>::iterator it = names.begin(); it != names.end(); ++it) {
    if (val == *it) return val;
  }

  if (exception) {
    char* tagname = mm.transcode(node->getNodeName());
    Errors::Message msg;
    msg << "Validation of content \"" << val << "\" for node \"" << tagname << "\" failed.\n";
    msg << "Available options: \"" << options << "\".\n";
    msg << "Please correct and try again.\n";
    Exceptions::amanzi_throw(msg);
  }

  return val;
}


/* ******************************************************************
* Parse input string using lists of available constants.
****************************************************************** */
std::string
InputConverter::GetConstantType_(const std::string& val_in, std::string& parsed_val)
{
  std::string val(val_in);
  trim(val);

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
int
InputConverter::GetPosition_(const std::vector<std::string>& names, const std::string& name)
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
* Useful function to get a root of domain subtree
****************************************************************** */
DOMNode*
InputConverter::GetRoot_(const std::string& domain, bool& flag)
{
  DOMNode* root;
  if (domain == "fracture") {
    root = GetUniqueElementByTagsString_("fracture_network", flag);
  } else {
    root = doc_->getDocumentElement();
    flag = true;
  }

  return root;
}


/* ******************************************************************
* Converts string of names separated by comma to array of strings.
****************************************************************** */
std::vector<std::string>
InputConverter::CharToStrings_(const char* namelist)
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
    trim(str);
    regs.push_back(str);
    tmp2 = strtok(NULL, ",");
  }

  delete[] tmp1;
  return regs;
}


/* ******************************************************************
* Extract unit and convert values. We assume that units are unique.
* If no unit specified, returned unit is the empty string.
****************************************************************** */
double
InputConverter::ConvertUnits_(const std::string& val, std::string& unit, double mol_mass)
{
  char* copy = strcpy(new char[val.size() + 1], val.c_str());
  char* data = strtok(copy, ";, ");

  double out;
  try {
    out = std::stod(std::string(data));
  } catch (...) {
    Errors::Message msg;
    msg << "\nInput string \"" << val << "\" cannot be converted to double + optional unit."
        << "\nThe string was parsed as \"" << data << "\".\n";
    Exceptions::amanzi_throw(msg);
  }
  data = strtok(NULL, ";,");

  // if units were found
  bool flag;
  unit = "";
  if (data != NULL) {
    unit = std::string(data);
    trim(unit);
    out = units_.ConvertUnitD(out, unit, "SI", mol_mass, flag);

    if (!flag) {
      Errors::Message msg;
      msg << "\nPrototype code for units of measurement cannot parse unit \"" << data << "\"."
          << "\n  (1) Try to convert a derived unit to atomic units, e.g. Pa=kg/m/s^2."
          << "\n  (2) Fix possible format errors, e.g. no space is allowed inside unit."
          << "\n  (3) Missing commas between coordinates.\n";
      Exceptions::amanzi_throw(msg);
    }
    found_units_.insert(data);
  }

  delete[] copy;
  return out;
}


/* ******************************************************************
* Extract unit and convert time values to seconds.
****************************************************************** */
double
InputConverter::TimeCharToValue_(const char* time_value)
{
  double time;
  char* tmp1 = strcpy(new char[strlen(time_value) + 1], time_value);
  char* tmp2 = strtok(tmp1, ";, ");

  time = std::strtod(tmp2, NULL);
  tmp2 = strtok(NULL, ";, ");

  // if units were found
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
std::vector<double>
InputConverter::MakeCoordinates_(const std::string& array)
{
  std::vector<double> coords;
  char* tmp1 = strcpy(new char[array.size() + 1], array.c_str());
  char* tmp2 = strtok(tmp1, "(,;");

  while (tmp2 != NULL) {
    std::string str(tmp2), parsed_str;
    trim(str);

    GetConstantType_(str, parsed_str);
    coords.push_back(std::strtod(parsed_str.c_str(), NULL));
    tmp2 = strtok(NULL, ",");
  }

  delete[] tmp1;
  return coords;
}


/* ******************************************************************
* Converts string of data to real numbers using units.
****************************************************************** */
std::vector<double>
InputConverter::MakeVector_(const std::string& array,
                            std::vector<std::string>& unit,
                            double mol_mass)
{
  std::vector<double> data;
  std::vector<std::string> tmp;

  tmp = CharToStrings_(array.c_str());

  unit.clear();
  for (int i = 0; i < tmp.size(); ++i) {
    std::string parsed_str, unit_in;
    GetConstantType_(tmp[i], parsed_str);
    data.push_back(ConvertUnits_(parsed_str, unit_in, mol_mass));
    unit.push_back(unit_in);
  }

  return data;
}


/* ******************************************************************
* Remove white spaces from both sides.
****************************************************************** */
std::string
InputConverter::TrimString_(char* tmp)
{
  std::string str(tmp);
  trim(str);
  return str;
}


/* *******************************************************************
* Generate error message when list is empty.
******************************************************************* */
int
InputConverter::IsEmpty(DOMNodeList* node_list, const std::string& name, bool exception)
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
void
InputConverter::ThrowErrorIllformed_(const std::string& section,
                                     const std::string& type,
                                     const std::string& ill_formed)
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
void
InputConverter::ThrowErrorIllformed_(const std::string& section,
                                     const std::string& type,
                                     const std::string& ill_formed,
                                     const std::string& options)
{
  Errors::Message msg;
  msg << "An error occurred during parsing node \"" << section << "\"\n";
  msg << "Missing or ill-formed " << type << " for \"" << ill_formed << "\"\n";
  msg << "Valid options are: " << options << "\n";
  msg << "Please correct and try again.\n";
  Exceptions::amanzi_throw(msg);
}


/* *******************************************************************
* Generate unified error message for missing item
******************************************************************* */
void
InputConverter::ThrowErrorMissing_(const std::string& node,
                                   const std::string& type,
                                   const std::string& key,
                                   const std::string& subnode)
{
  Errors::Message msg;
  msg << "An error occurred during parsing high-level node \"" << node << "\"\n";
  msg << "No " << type << " \"" << key << "\" found for sub-node \"" << subnode << "\".\n";
  Exceptions::amanzi_throw(msg);
}


/* *******************************************************************
* Generate unified error message for missing child
******************************************************************* */
void
InputConverter::ThrowErrorMisschild_(const std::string& section,
                                     const std::string& missing,
                                     const std::string& name)
{
  Errors::Message msg;
  msg << "Amanzi::InputConverter: an error occurred during parsing node \"" << section << "\"\n";
  msg << "  No child \"" << missing << "\" found";
  if (!name.empty()) { msg << " for \"" << name << "\""; }
  msg << ".\n";
  msg << "  Please correct and try again \n";
  Exceptions::amanzi_throw(msg);
}


/* ***************************************************************************
* Extracts information for and write the Pflotran Chemistry Engine input file
* Returns the name of this file.
*************************************************************************** */
std::string
InputConverter::CreateINFile_(std::string& filename, int rank)
{
  MemoryManager mm;
  DOMNode* node;
  DOMNode* base;
  DOMElement* element;

  std::ofstream in_file;
  std::stringstream primaries;
  std::stringstream secondaries;
  std::stringstream redoxes;
  std::stringstream minerals;
  std::stringstream mineral_list;
  std::stringstream gases;
  std::stringstream mineral_kinetics;
  std::stringstream isotherms;
  std::stringstream cations;
  std::stringstream complexes;
  std::stringstream constraints;
  std::stringstream reactionrates;
  std::stringstream decayrates;
  std::stringstream controls;

  int first_cation;

  // database filename and controls
  bool flag;
  node = GetUniqueElementByTagsString_("process_kernels, chemistry", flag);
  element = static_cast<DOMElement*>(node);
  std::string datfilename = GetAttributeValueS_(element, "database", TYPE_NONE, true, "");

  struct stat buffer;
  int status = stat(datfilename.c_str(), &buffer);
  if (status == -1) {
    Errors::Message msg("The database file '" + datfilename + "' is missing.");
    Exceptions::amanzi_throw(msg);
  }

  // add relative path from xmlfilename_ to datfilename (simplified code)
  size_t pos0;
  std::string path(xmlfilename_);
  if ((pos0 = path.find_last_of('/')) == std::string::npos)
    path = "./";
  else
    path.erase(path.begin() + pos0, path.end());

  controls << "  DATABASE " << boost::filesystem::relative(datfilename, path).string().c_str()
           << "\n";

  base = GetUniqueElementByTagsString_(
    "numerical_controls, unstructured_controls, unstr_chemistry_controls", flag);
  if (flag) {
    node = GetUniqueElementByTagsString_(base, "activity_coefficients", flag);
    std::string tmp("  ACTIVITY_COEFFICIENTS TIMESTEP");
    if (flag) {
      std::string value = TrimString_(mm.transcode(node->getTextContent()));
      if (value == "off") { tmp = "  ACTIVITY_COEFFICIENTS OFF"; }
    }
    controls << tmp << "\n";

    node = GetUniqueElementByTagsString_(base, "log_formulation", flag);
    if (flag) {
      std::string value = TrimString_(mm.transcode(node->getTextContent()));
      if (value == "on") { controls << "  LOG_FORMULATION \n"; }
    } else {
      controls << "  LOG_FORMULATION \n";
    }

    node = GetUniqueElementByTagsString_(base, "max_relative_change_tolerance", flag);
    if (flag) {
      std::string value = TrimString_(mm.transcode(node->getTextContent()));
      controls << "  MAX_RELATIVE_CHANGE_TOLERANCE " << value << "\n";
    }

    node = GetUniqueElementByTagsString_(base, "max_residual_tolerance", flag);
    if (flag) {
      std::string value = TrimString_(mm.transcode(node->getTextContent()));
      controls << "  MAX_RESIDUAL_TOLERANCE " << value << "\n";
    }

    node = GetUniqueElementByTagsString_(base, "use_full_geochemistry", flag);
    if (flag) {
      std::string value = TrimString_(mm.transcode(node->getTextContent()));
      if (value == "on") { controls << "  USE_FULL_GEOCHEMISTRY \n"; }
    } else {
      controls << "  USE_FULL_GEOCHEMISTRY \n";
    }
  }

  // set up Chemistry Options ParameterList
  Teuchos::ParameterList ChemOptions;
  ChemOptions.sublist("isotherms");
  ChemOptions.sublist("ion_exchange");
  ChemOptions.sublist("surface_complexation");

  // get primary species names
  // and check for forward/backward rates (for solutes/non-reactive primaries only)
  std::map<std::string, int> incomplete_tracers;

  node = GetUniqueElementByTagsString_("liquid_phase, dissolved_components, primaries", flag);
  if (flag) {
    std::string name, primary;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, primary, flag, false);

    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      name = TrimString_(mm.transcode(inode->getTextContent()));
      primaries << "    " << name << "\n";
      element = static_cast<DOMElement*>(inode);

      if (element->hasAttribute(mm.transcode("first_order_decay_constant"))) {
        double decay = GetAttributeValueD_(element, "first_order_decay_constant");
        name = TrimString_(mm.transcode(inode->getTextContent()));
        decayrates << "    REACTION " << name << " <-> \n";
        decayrates << "    RATE_CONSTANT " << decay << "\n";
      } else {
        // verify that species is non-reactive (only one record (1 line) in database file)
        int ncount = CountFileLinesWithWord_(datfilename, name);
        if (ncount == 1 && element->hasAttribute(mm.transcode("forward_rate"))) {
          double frate = GetAttributeValueD_(element, "forward_rate");
          double brate = GetAttributeValueD_(element, "backward_rate");
          name = TrimString_(mm.transcode(inode->getTextContent()));
          reactionrates << "    REACTION " << name << " <->\n";
          reactionrates << "    FORWARD_RATE " << frate << "\n";
          reactionrates << "    BACKWARD_RATE " << brate << "\n";
        } else {
          // not yet a bug: species may be involved in sorption process
          incomplete_tracers[name];
        }
      }
    }
  }

  // get secondaries
  node = GetUniqueElementByTagsString_("liquid_phase, dissolved_components, secondaries", flag);
  if (flag) {
    std::string secondary;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, secondary, flag, false);

    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      std::string name = TrimString_(mm.transcode(inode->getTextContent()));
      secondaries << "    " << name << "\n";
    }
  }

  // get redox species
  node = GetUniqueElementByTagsString_("liquid_phase, dissolved_components, redox", flag);
  if (flag) {
    std::string redox;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, redox, flag, false);

    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      std::string name = TrimString_(mm.transcode(inode->getTextContent()));
      redoxes << "    " << name << "\n";
    }
  }

  // get minerals and mineral kinetics
  node = GetUniqueElementByTagsString_("phases, solid_phase, minerals", flag);
  if (flag) {
    std::string mineral;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, mineral, flag, false);

    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      std::string name = TrimString_(mm.transcode(inode->getTextContent()));
      minerals << "    " << name << "\n";
      mineral_list << name << ", ";

      element = static_cast<DOMElement*>(inode);
      double rate = GetAttributeValueD_(
        element, "rate_constant", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "", false, 0.0);
      // double val = GetAttributeValueD_(element, "rate_constant", TYPE_NUMERICAL, 0.0, DVAL_MAX, "mol/m^2/s", false, 0.0);
      // double rate = units_.ConvertUnitD(val, "mol/m^2/s", "mol/cm^2/s", -1.0, flag);

      mineral_kinetics << "    " << name << "\n";
      mineral_kinetics << "      RATE_CONSTANT " << rate << "\n";
      //<< " mol/cm^2-sec\n";

      if (element->hasAttribute(mm.transcode("prefactor_species"))) {
        std::string species = GetAttributeValueS_(element, "prefactor_species");
        double alpha = GetAttributeValueD_(element, "alpha");
        mineral_kinetics << "      PREFACTOR\n";
        mineral_kinetics << "        RATE_CONSTANT " << rate << "\n";
        //" mol/cm^2-sec\n";
        mineral_kinetics << "        PREFACTOR_SPECIES " << species << "\n";
        mineral_kinetics << "          ALPHA " << alpha << "\n";
        mineral_kinetics << "        /\n";
        mineral_kinetics << "      /\n";
      }

      mineral_kinetics << "    /\n";
    }
  }

  // get gases
  node = GetUniqueElementByTagsString_("phases, gas_phase, gases", flag);
  if (flag) {
    std::string gas;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, gas, flag, false);

    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      std::string name = TrimString_(mm.transcode(inode->getTextContent()));
      gases << "    " << name << "\n";
    }
  }

  // gather material specific chemistry options and put into ParameterList
  node = GetUniqueElementByTagsString_("materials", flag);
  if (flag) {
    // loop over materials to grab Kd, cations, and surface complexation
    std::string name;
    std::vector<DOMNode*> mat_list = GetSameChildNodes_(node, name, flag, false);
    for (int i = 0; i < mat_list.size(); ++i) {
      bool flag2;
      DOMElement* child_elem;

      DOMNode* inode = mat_list[i];
      element = static_cast<DOMElement*>(inode);
      name = GetAttributeValueS_(element, "name");

      // look for sorption_isotherms
      child_elem = GetChildByName_(inode, "sorption_isotherms", flag2, false);
      if (flag2) {
        // loop over sublist of primaries to get Kd information
        std::vector<DOMNode*> primary_list =
          GetSameChildNodes_(static_cast<DOMNode*>(child_elem), name, flag2, false);
        if (flag2) {
          for (int j = 0; j < primary_list.size(); ++j) {
            DOMNode* jnode = primary_list[j];
            std::string primary_name = GetAttributeValueS_(jnode, "name");
            DOMNode* kd_node = GetUniqueElementByTagsString_(jnode, "kd_model", flag2);
            DOMElement* kd_elem = static_cast<DOMElement*>(kd_node);

            if (flag2) {
              Teuchos::ParameterList kd_list;
              double kd = GetAttributeValueD_(kd_elem, "kd");
              std::string model = GetAttributeValueS_(kd_elem, "model");
              kd_list.set<std::string>("model", model);
              kd_list.set<double>("kd", kd);

              if (model == "langmuir") {
                double b = GetAttributeValueD_(kd_elem, "b");
                kd_list.set<double>("b", b);
              } else if (model == "freundlich") {
                double n = GetAttributeValueD_(kd_elem, "n");
                kd_list.set<double>("n", n);
              }

              ChemOptions.sublist("isotherms").sublist(primary_name) = kd_list;
            }
          }
        }
      }

      // look for ion exchange
      child_elem = GetChildByName_(inode, "ion_exchange", flag2, false);
      if (flag2) {
        // loop over sublist of cations to get information
        std::vector<DOMNode*> mineral_nodes =
          GetSameChildNodes_(static_cast<DOMNode*>(child_elem), name, flag2, false);
        if (flag2) {
          for (int j = 0; j < mineral_nodes.size(); ++j) {
            DOMNode* jnode = mineral_nodes[j];
            Teuchos::ParameterList ion_list;
            double cec = GetAttributeValueD_(jnode, "cec");
            ion_list.set<double>("cec", cec);

            // loop over list of cation names/selectivity pairs
            std::vector<DOMNode*> cations_list = GetSameChildNodes_(jnode, name, flag2, false);
            if (flag2) {
              std::vector<std::string> cation_names;
              std::vector<double> cation_selectivity;

              for (int k = 0; k < cations_list.size(); ++k) {
                DOMNode* knode = cations_list[k];
                DOMElement* kelement = static_cast<DOMElement*>(knode);
                cation_names.push_back(GetAttributeValueS_(kelement, "name"));

                double value = GetAttributeValueD_(kelement, "value");
                cation_selectivity.push_back(value);
                if (value == 1.0) { first_cation = k; }
              }

              ion_list.set<Teuchos::Array<std::string>>("cations", cation_names);
              ion_list.set<Teuchos::Array<double>>("values", cation_selectivity);
            }

            ChemOptions.sublist("ion_exchange").sublist("bulk") = ion_list;
          }
        }
      }

      // look for surface complexation
      child_elem = GetChildByName_(inode, "surface_complexation", flag2, false);
      if (flag2) {
        // loop over sublist of cations to get information
        std::vector<DOMNode*> site_list =
          GetSameChildNodes_(static_cast<DOMNode*>(child_elem), name, flag2, false);
        if (flag2) {
          for (int j = 0; j < site_list.size(); ++j) {
            DOMNode* jnode = site_list[j];
            Teuchos::ParameterList surface_list;
            std::string site = GetAttributeValueS_(jnode, "name");
            double density = GetAttributeValueD_(jnode, "density");
            surface_list.set<double>("density", density);

            std::string name2;
            std::vector<DOMNode*> children = GetSameChildNodes_(jnode, name2, flag2, false);
            Teuchos::Array<std::string> complexe_names =
              CharToStrings_(mm.transcode(children[0]->getTextContent()));
            surface_list.set<Teuchos::Array<std::string>>("complexes", complexe_names);

            ChemOptions.sublist("surface_complexation").sublist(site) = surface_list;
          }
        }
      }

      // look for minerals
      child_elem = GetChildByName_(inode, "minerals", flag2, false);
      if (flag2) {
        // loop over sublist of minerals to get information
        std::vector<DOMNode*> minerals_list =
          GetSameChildNodes_(static_cast<DOMNode*>(child_elem), name, flag2, false);
        if (flag2) {
          for (int j = 0; j < minerals_list.size(); ++j) {
            DOMNode* jnode = minerals_list[j];
            Teuchos::ParameterList mineral_out;
            std::string mineral_name = GetAttributeValueS_(jnode, "name");
            double volume_fraction = GetAttributeValueD_(jnode, "volume_fraction");
            double specific_surface_area = GetAttributeValueD_(jnode, "specific_surface_area");
            mineral_out.set<double>("volume_fraction", volume_fraction);
            mineral_out.set<double>("specific_surface_area", specific_surface_area);

            ChemOptions.sublist("minerals").sublist(mineral_name) = mineral_out;
          }
        }
      }
    }
  }

  // create text for chemistry options - isotherms
  Teuchos::ParameterList& iso = ChemOptions.sublist("isotherms");

  for (auto iter = iso.begin(); iter != iso.end(); ++iter) {
    std::string primary = iso.name(iter);
    Teuchos::ParameterList& curprimary = iso.sublist(primary);
    std::string model = curprimary.get<std::string>("model");
    double kd = curprimary.get<double>("kd");

    if (model == "linear") {
      isotherms << "      " << primary << "\n";
      isotherms << "        TYPE LINEAR\n";
      isotherms << "        DISTRIBUTION_COEFFICIENT " << kd << "\n";
    } else if (model == "langmuir") {
      isotherms << "      " << primary << "\n";
      isotherms << "        TYPE LANGMUIR\n";
      isotherms << "        DISTRIBUTION_COEFFICIENT " << kd << "\n";
      double b = curprimary.get<double>("b");
      isotherms << "        LANGMUIR_B " << b << "\n";
    } else if (model == "freundlich") {
      isotherms << "      " << primary << "\n";
      isotherms << "        TYPE FREUNDLICH\n";
      isotherms << "        DISTRIBUTION_COEFFICIENT " << kd << "\n";
      double n = curprimary.get<double>("n");
      isotherms << "        FREUNDLICH_N " << n << "\n";
    }
    isotherms << "      /\n";
  }

  // create text for chemistry options - ion exchange
  Teuchos::ParameterList& ion = ChemOptions.sublist("ion_exchange");
  for (auto iter = ion.begin(); iter != ion.end(); ++iter) {
    std::string mineral = ion.name(iter);
    Teuchos::ParameterList& curmineral = ion.sublist(mineral);
    double cec = curmineral.get<double>("cec");
    Teuchos::Array<std::string> names = curmineral.get<Teuchos::Array<std::string>>("cations");
    Teuchos::Array<double> values = curmineral.get<Teuchos::Array<double>>("values");

    cations << "      CEC " << cec << "\n";
    cations << "      CATIONS\n";
    // print primary cation first (this matters to pflotran)
    cations << "        " << names[first_cation] << " " << values[first_cation] << " REFERENCE\n";
    for (int i = 0; i < names.size(); i++) {
      if (i != first_cation) cations << "        " << names[i] << " " << values[i] << "\n";
    }
  }

  // create text for chemistry options - surface complexation
  Teuchos::ParameterList& surf = ChemOptions.sublist("surface_complexation");
  for (auto iter = surf.begin(); iter != surf.end(); ++iter) {
    std::string site = surf.name(iter);
    Teuchos::ParameterList& cursite = surf.sublist(site);
    Teuchos::Array<std::string> complexe_names =
      cursite.get<Teuchos::Array<std::string>>("complexes");
    double density = cursite.get<double>("density");
    complexes << "    SURFACE_COMPLEXATION_RXN\n";
    complexes << "      EQUILIBRIUM\n";
    complexes << "      SITE " << site << " " << density << "\n";
    complexes << "      COMPLEXES\n";
    for (Teuchos::Array<std::string>::iterator it = complexe_names.begin();
         it != complexe_names.end();
         it++) {
      complexes << "        " << *it << "\n";
    }
    complexes << "      /\n";
    complexes << "    /\n";
  }

  // create text for minerals (read from materials section and written out to constraints section
  std::stringstream mineral;
  Teuchos::ParameterList& mins = ChemOptions.sublist("minerals");
  for (auto iter = mins.begin(); iter != mins.end(); ++iter) {
    std::string mineral_name = mins.name(iter);
    Teuchos::ParameterList& curmineral = mins.sublist(mineral_name);
    double constvf = curmineral.get<double>("volume_fraction");
    double constsa = curmineral.get<double>("specific_surface_area");
    mineral << "    " << mineral_name << " " << constvf << " " << constsa << "\n";
  }

  // constraints
  node = GetUniqueElementByTagsString_("geochemistry, constraints", flag);
  if (flag) {
    // loop over constraints
    std::string name;
    std::vector<DOMNode*> constraint_list = GetSameChildNodes_(node, name, flag, false);
    for (int i = 0; i < constraint_list.size(); ++i) {
      std::stringstream primary;
      std::stringstream gas;

      DOMNode* inode = constraint_list[i];
      element = static_cast<DOMElement*>(inode);
      std::string const_name = GetAttributeValueS_(element, "name");

      // loop over individual items in current constraint
      std::string constraint;
      bool found;
      std::vector<DOMNode*> children = GetChildren_(inode, "primary", found);
      if (found) {
        for (int j = 0; j < children.size(); ++j) {
          DOMNode* jnode = children[j];
          element = static_cast<DOMElement*>(jnode);

          std::string constname = GetAttributeValueS_(element, "name");
          std::string ctype = GetAttributeValueS_(element, "type");
          double constvalue = GetAttributeValueD_(element, "value");
          std::string typeletter;
          if (ctype == "free_ion") {
            typeletter = "F";
          } else if (ctype == "pH" || ctype == "ph") {
            typeletter = "P";
          } else if (ctype == "total") {
            typeletter = "T";
          } else if (ctype == "total+sorbed") {
            typeletter = "S";
          } else if (ctype == "charge") {
            typeletter = "Z";
          } else if (ctype == "mineral") {
            std::string mname = GetAttributeValueS_(element, "mineral");
            typeletter = "M " + mname;
          } else if (ctype == "gas") {
            std::string gname = GetAttributeValueS_(element, "gas");
            typeletter = "G " + gname;
          }
          primary << "    " << constname << " " << constvalue << " " << typeletter << "\n";
        }
      }

      children = GetChildren_(inode, "mineral", found);
      if (found) {
        for (int j = 0; j < children.size(); ++j) {
          DOMNode* jnode = children[j];
          element = static_cast<DOMElement*>(jnode);
          std::string constname = GetAttributeValueS_(element, "name");
          double constvf = GetAttributeValueD_(element, "volume_fraction");
          double constsa = GetAttributeValueD_(element, "specific_surface_area");
          mineral << "    " << constname << " " << constvf << " " << constsa << "\n";
        }
      }

      children = GetChildren_(inode, "gas", found);
      if (found) {
        for (int j = 0; j < children.size(); ++j) {
          DOMNode* jnode = children[j];
          element = static_cast<DOMElement*>(jnode);
          std::string constname = GetAttributeValueS_(element, "name");
          gas << "    " << constname << "\n";
        }
      }

      constraints << "CONSTRAINT " << const_name << "\n";
      constraints << "  CONCENTRATIONS\n";
      constraints << primary.str();
      constraints << "  /\n";
      if (!mineral.str().empty()) {
        constraints << "  MINERALS\n";
        constraints << mineral.str();
        constraints << "  /\n";
      }
      constraints << "END\n";
    }
  }


  // build in filename
  std::string infilename(filename);
  std::string new_extension(".in");
  size_t pos = infilename.find(".xml");
  infilename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)4);

  // open output in file
  if (rank == 0) {
    // struct stat buffer;
    status = stat(infilename.c_str(), &buffer);
    // TODO EIB - fix this
    /*
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "Writing PFloTran partial input file \"" <<
          infilename.c_str() << "\"" << std::endl;
    }
    */
    std::cout << "Writing PFloTran partial input file " << infilename.c_str() << std::endl;

    in_file.open(infilename.c_str());

    in_file << "CHEMISTRY\n";
    // Chemistry - Species
    in_file << "  PRIMARY_SPECIES\n";
    in_file << primaries.str();
    in_file << "  /\n";
    if (!reactionrates.str().empty()) {
      in_file << "  GENERAL_REACTION\n";
      in_file << reactionrates.str();
      in_file << "  /\n";
    }
    if (!decayrates.str().empty()) {
      in_file << "  RADIOACTIVE_DECAY_REACTION\n";
      in_file << decayrates.str();
      in_file << "  /\n";
    }
    if (!secondaries.str().empty()) {
      in_file << "  SECONDARY_SPECIES\n";
      in_file << secondaries.str();
      in_file << "  /\n";
    }
    if (!redoxes.str().empty()) {
      in_file << "  DECOUPLED_EQUILIBRIUM_REACTIONS\n";
      in_file << redoxes.str();
      in_file << "  /\n";
    }
    if (!gases.str().empty()) {
      in_file << "  PASSIVE_GAS_SPECIES\n";
      in_file << gases.str();
      in_file << "  /\n";
    }
    if (!minerals.str().empty()) {
      in_file << "  MINERALS\n";
      in_file << minerals.str();
      in_file << "  /\n";
    }

    // Chemistry - Mineral Kinetics
    if (!mineral_kinetics.str().empty()) {
      in_file << "  MINERAL_KINETICS\n";
      in_file << mineral_kinetics.str();
      in_file << "  /\n";
    }

    // Chemistry - Sorption
    if (!isotherms.str().empty() || !complexes.str().empty() || !cations.str().empty()) {
      in_file << "  SORPTION\n";
      if (!isotherms.str().empty()) {
        in_file << "    ISOTHERM_REACTIONS\n";
        in_file << isotherms.str();
        in_file << "    /\n";
      }
      if (!complexes.str().empty()) { in_file << complexes.str(); }
      if (!cations.str().empty()) {
        in_file << "    ION_EXCHANGE_RXN\n";
        in_file << cations.str();
        in_file << "      /\n";
        in_file << "    /\n";
      }
      in_file << "  /\n";
    }

    // Chemistry - Controls
    in_file << controls.str();
    in_file << "END\n";

    // Constraints
    in_file << constraints.str();
    //in_file << "END\n";

    in_file.close();
  }

  // prints warnings
  if (incomplete_tracers.size() > 0) {
    std::cout << "WARNING: primary species that may have missing atributes: ";
    for (auto it = incomplete_tracers.begin(); it != incomplete_tracers.end(); ++it) {
      std::cout << it->first << ", ";
    }
    std::cout << std::endl;
  }

  return infilename;
}


/* ******************************************************************
* Returns number of lines that use word in the file.
****************************************************************** */
int
InputConverter::CountFileLinesWithWord_(const std::string& filename, const std::string& word)
{
  std::ifstream file;
  file.open(filename);

  int n = 0;
  while (!file.eof()) {
    std::string line;
    file >> line;
    if (line.find(word) != std::string::npos) n++;
  }
  file.close();

  return n;
}

} // namespace AmanziInput
} // namespace Amanzi
