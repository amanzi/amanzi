/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_INPUT_CONVERTER_HH_
#define AMANZI_INPUT_CONVERTER_HH_

#include <set>
#include <climits>
#include <list>

#include "boost/lambda/lambda.hpp"
#include "boost/bind.hpp"
#include "boost/lexical_cast.hpp"

#define  BOOST_FILESYTEM_NO_DEPRECATED
#include "boost/format.hpp"

// Xerces
#include "xercesc/dom/DOM.hpp"
#include "xercesc/parsers/XercesDOMParser.hpp"

// Amanzi's
#include "Units.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziInput {

// Amanzi version
#define AMANZI_SPEC_VERSION_MAJOR 2
#define AMANZI_SPEC_VERSION_MINOR 3
#define AMANZI_SPEC_VERSION_MICRO 0

// constants
const std::string TYPE_TIME = "time";
const std::string TYPE_NUMERICAL = "numerical";
const std::string TYPE_AREA_MASS_FLUX = "area_mass_flux";
const std::string TYPE_NONE = "none";
const std::string TYPE_NOT_CONSTANT = "not_constant";

const double DVAL_MIN = -1.0e+99;
const double DVAL_MAX = 1.0e+99;

XERCES_CPP_NAMESPACE_USE

/* 
* A simple wrapper for XMLString class. It collects memory pointers
* and destroys them later. The focus is on simplicity of its using,
* so that release is never called.
*/
class MemoryManager {
 public:
  MemoryManager() {};
  ~MemoryManager() { Destroy_(); }

  XMLCh* transcode(const char* str) {
    XMLCh* xstr = xercesc::XMLString::transcode(str);
    xchar.push_back(xstr);
    return xstr;
  }

  char* transcode(const XMLCh* xstr) {
    char* str = xercesc::XMLString::transcode(xstr);
    pchar.push_back(str);
    return str;
  }

 private:
  void Destroy_() {
    for (auto it = pchar.begin(); it != pchar.end(); ++it)
      xercesc::XMLString::release(&*it);
    for (auto it = xchar.begin(); it != xchar.end(); ++it)
      xercesc::XMLString::release(&*it);
  }

 private:
  std::list<char*> pchar;
  std::list<XMLCh*> xchar;
};


//------------------------------------------------------------------------
// XML helper methods for parsing outside of an InputConverter:

// Creates an XML parser with our desired settings. This parser must be deleted
// after the document has been used.
XercesDOMParser* CreateXMLParser();

// Using the given XML parser, parses the document contained in the file with 
// the given name.
DOMDocument* OpenXMLInput(XercesDOMParser* parser,
                          const std::string& xml_input);
//------------------------------------------------------------------------


class InputConverter {
 public:
  // This constructor opens up the file with the given name, sets up a parser,
  // and parses the file.
  explicit InputConverter(const std::string& input_filename);

  // This constructor uses an already-parsed XML document, and does not 
  // manage the parser.
  InputConverter(const std::string& input_filename, DOMDocument* input_doc);

  virtual ~InputConverter();

  // parse various nodes
  void ParseVersion_();
  void ParseConstants_();
  void FilterNodes(const std::string& filter);

  // auto-generated input files
  // -- native chemistry
  std::string CreateBGDFile_(std::string& filename, int rank, int& status);
  // -- Pflotran input file
  std::string CreateINFile_(std::string& filename, int rank);

 protected:
  // Useful tools wrapping low-level DOM commands

  // -- generalization of getElementsByTagNames(): returns node
  //    where all tags (list of names in the given string) are unique 
  //    leaves of the tree.
  xercesc::DOMNode* GetUniqueElementByTagsString_(
      const std::string& tags, bool& flag);

  // -- modification of the previous routine where the first tag 
  //    is replaced by a pointer to document's element
  xercesc::DOMNode* GetUniqueElementByTagsString_(
      const xercesc::DOMNode* node1, const std::string& tags, bool& flag);

  // -- extract child with the given attribute
  xercesc::DOMElement* GetUniqueChildByAttribute_(
      xercesc::DOMNode* node, const char* attr_name, const std::string& attr_value,
      bool& flag, bool exception = false);

  // -- extracts the child of the given node with the given name.
  xercesc::DOMElement* GetChildByName_(
      xercesc::DOMNode* node, const std::string& childName, bool& flag, bool exception = false);

  // -- extract and verify children
  // -- extract existing attribute value and verify it optionally against expected units
  //    Consider two examples <parameters alpha="2.06e-03 Pa^-1"/>  and
  //                          <parameters alpha=ALPHA/>
  //    elem      = pointer to an element <.../> in an XML document
  //    attr_name = name of the attribute  ("alpha" in the examples)
  //    type      = value type. If the value (ALPHA) is defined in the XML section "constants",
  //                then type should match the type of constant as follows: TYPE_TIME for 
  //                time_constant TYPE_NUMERICAL for numerical_constant and TYPE_AREA_MASS_FLUX 
  //                for area_mass_flux_constant. Use default value, otherwise.
  //    valmin    = minimum allowed value in SI units (2.06e-03 >= valmin)
  //    valmax    = maximum allowed value in SI units (2.06e-03 <= valmax) 
  //    unit      = check expected units if provided. Derived and non-SI units are allowed.
  //    length    = expected length of vector. If size is unknown, use -1.
  //    exception = if false, use default value when either element or attribute is missing.
  //                if true, throw an exception.
  //    default_val = defult value
  int GetAttributeValueL_(
      xercesc::DOMElement* elem, const char* attr_name, const std::string& type = TYPE_NUMERICAL,
      int valmin = INT_MIN, int valmax = INT_MAX, bool exception = true, int defaul_val = 0);
  double GetAttributeValueD_(  // supports units except for ppbm
      xercesc::DOMElement* elem, const char* attr_name, const std::string& type = TYPE_NUMERICAL,
      double valmin = DVAL_MIN, double valmax = DVAL_MAX, std::string unit = "",
      bool exception = true, double default_val = 0.0);
  std::string GetAttributeValueS_(
      xercesc::DOMElement* elem, const char* attr_name, const std::string& type = TYPE_NUMERICAL,
      bool exception = true, std::string default_val = "");
  std::vector<double> GetAttributeVectorD_(  // supports units except ppbm
      xercesc::DOMElement* elem, const char* attr_name, int length = -1,
      std::string unit = "", bool exception = true, double mol_mass = -1.0);
  std::vector<std::string> GetAttributeVectorS_(
      xercesc::DOMElement* elem, const char* attr_name, bool exception = true);

  // -- node is used more often then element
  int GetAttributeValueL_(
      xercesc::DOMNode* node, const char* attr_name,
      const std::string& type = TYPE_NUMERICAL, int valmin = INT_MIN, int valmax = INT_MAX,
       bool exception = true, int val = 0) {
    xercesc::DOMElement* element = static_cast<xercesc::DOMElement*>(node);
    return GetAttributeValueL_(element, attr_name, type, valmin, valmax, exception, val);
  }
  double GetAttributeValueD_(  // supports units except for ppbm
      xercesc::DOMNode* node, const char* attr_name, const std::string& type = TYPE_NUMERICAL,
      double valmin = DVAL_MIN, double valmax = DVAL_MAX, std::string unit = "",
      bool exception = true, double val = 0.0) {
    xercesc::DOMElement* element = static_cast<xercesc::DOMElement*>(node);
    return GetAttributeValueD_(element, attr_name, type, valmin, valmax, unit, exception, val);
  }
  std::string GetAttributeValueS_(
      xercesc::DOMNode* node, const char* attr_name,
      const std::string& type = TYPE_NUMERICAL, bool exception = true, std::string val = "") {
    xercesc::DOMElement* element = static_cast<xercesc::DOMElement*>(node);
    return GetAttributeValueS_(element, attr_name, type, exception, val);
  }
  std::vector<double> GetAttributeVectorD_(  // supports units except ppbm
      xercesc::DOMNode* node, const char* attr_name, int length = -1,
      std::string unit = "", bool exception = true, double mol_mass = -1.0) {
    xercesc::DOMElement* element = static_cast<xercesc::DOMElement*>(node);
    return GetAttributeVectorD_(element, attr_name, length, unit, exception, mol_mass);
  }
  std::vector<std::string> GetAttributeVectorS_(
      xercesc::DOMNode* node, const char* attr_name, bool exception = true) {
    xercesc::DOMElement* element = static_cast<xercesc::DOMElement*>(node);
    return GetAttributeVectorS_(element, attr_name, exception);
  }
 
  // -- extract the text content of the child of the given node with the given name.
  std::string GetChildValueS_(
      xercesc::DOMNode* node, const std::string& childName, bool& flag, bool exception = false);
  std::vector<std::string> GetChildVectorS_(
      xercesc::DOMNode* node, const std::string& childName, bool& flag, bool exception = false);

  // -- extract existing attribute value and verify it
  std::string GetAttributeValueS_(
      xercesc::DOMNode* node, const char* attr_name, const char* options);

  // -- extract all children of the given node that share the given common name.
  std::vector<xercesc::DOMNode*> GetChildren_(
      xercesc::DOMNode* node, const std::string& childrenName, bool& flag, bool exception = false);

  //    the name of identical nodes will be extracted too
  std::vector<xercesc::DOMNode*> GetSameChildNodes_(
      xercesc::DOMNode* node, std::string& name, bool& flag, bool exception = false);

  // -- extract text content and verify it against list of options
  double GetTextContentD_(  // supports units
      xercesc::DOMNode* node, std::string unit = "", bool exception = false, double val = 0.0);
  std::string GetTextContentS_(
      xercesc::DOMNode* node, const char* options, bool exception = true);

  // data streaming/trimming/converting
  // -- units. Molar mass is required for converting ppm and ppb units.
  double ConvertUnits_(const std::string& val, std::string& unit, double mol_mass = -1.0);

  // -- times
  double TimeCharToValue_(const char* time_value);
  std::string GetConstantType_(const std::string& val, std::string& parsed_val);

  // -- coordinates and vectors.
  std::vector<double> MakeCoordinates_(const std::string& array);
  std::vector<double> MakeVector_(  // support units
      const std::string& array, std::vector<std::string>& unit, double mol_mass = -1.0);

  // -- string modifications
  std::vector<std::string> CharToStrings_(const char* namelist);
  std::string TrimString_(char* tmp);

  // -- vector parsing
  int GetPosition_(const std::vector<std::string>& names, const std::string& name);

  // -- test file parsing 
  int CountFileLinesWithWord_(const std::string& filename, const std::string& word);

  // is spec structurally sound?
  // -- mandatory objects
  int IsEmpty(xercesc::DOMNodeList* node_list, const std::string& name, bool exception = true);
   
  // -- error messages
  void ThrowErrorIllformed_(
      const std::string& section, const std::string& type, const std::string& ill_formed);
  void ThrowErrorIllformed_(
      const std::string& section, const std::string& type, const std::string& ill_formed, const std::string& options);
  void ThrowErrorMissing_(
      const std::string& node, const std::string& type, const std::string& key, const std::string& subnode);
  void ThrowErrorMisschild_(
    const std::string& section, const std::string& missing, const std::string& name = std::string());

 protected:
  // various constants defined by the users
  // consistency check is performed for all but constants_
  std::map<std::string, std::string> constants_time_;  
  std::map<std::string, std::string> constants_numerical_; 
  std::map<std::string, std::string> constants_area_mass_flux_; 
  std::map<std::string, std::string> constants_;  // no check

  Utils::Units units_;

  std::string xmlfilename_;
  DOMDocument* doc_;
  XercesDOMParser* parser_;

 private:
  // Disallowed deep-copy-related methods.
  InputConverter(const InputConverter&);
  InputConverter& operator=(const InputConverter&);

 protected:
  // statistics
  std::set<std::string> found_units_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
