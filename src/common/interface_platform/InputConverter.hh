/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_INPUT_CONVERTER_HH_
#define AMANZI_INPUT_CONVERTER_HH_

#include "boost/lambda/lambda.hpp"
#include "boost/bind.hpp"
#include "boost/lexical_cast.hpp"

#define  BOOST_FILESYTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/format.hpp"

// TPLs
#include "xercesc/dom/DOM.hpp"
#include "xercesc/parsers/XercesDOMParser.hpp"

// Amanzi's
#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziInput {

/* 
* A simple wrapper for XMLString class. It collects memory pointers
* and destroys them later. The focus is on simplicity of its using,
* so that release is never called.
*/
class XString {
 public:
  XString() {};
  ~XString() { Destroy_(); }

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
    for (std::vector<char*>::iterator it = pchar.begin(); it != pchar.end(); ++it)
      xercesc::XMLString::release(&*it);
    for (std::vector<XMLCh*>::iterator it = xchar.begin(); it != xchar.end(); ++it)
      xercesc::XMLString::release(&*it);
  }

 private:
  std::vector<char*> pchar;
  std::vector<XMLCh*> xchar;
};

class InputConverter {
 public:
  InputConverter() {
    xercesc::XMLPlatformUtils::Initialize();
  }

  ~InputConverter() {
    delete parser;
    xercesc::XMLPlatformUtils::Terminate();
  }

  // main member: creates xerces document using the file name
  void Init(const std::string& xmlfilename);

  // parse various nodes
  void ParseConstants_();

 protected:
  // Useful tools wrapping low-level DOM commands
  // -- generalization of getElementsByTagNames(): returns node
  //    tag1->tag2 or tag1->tag2-tag3 where all tags are unique 
  //    leaves of the tree.
  xercesc::DOMNode* getUniqueElementByTagNames_(
      const std::string& tag1, const std::string& tag2, bool& flag);
  xercesc::DOMNode* getUniqueElementByTagNames_(
      const std::string& tag1, const std::string& tag2, const std::string& tag3, bool& flag);
  // -- tags contains list of names separated by commas. It 
  //    will replace eventually the previous routine.
  xercesc::DOMNode* getUniqueElementByTagsString_(
      const std::string& tags, bool& flag);

  // -- modification of the previous routines where the first tag 
  //    is replaced by a pointer to document's element
  xercesc::DOMNode* getUniqueElementByTagNames_(
      const xercesc::DOMNode* node1, const std::string& tag2, bool& flag);
  xercesc::DOMNode* getUniqueElementByTagNames_(
      const xercesc::DOMNode* node1, const std::string& tag2, const std::string& tag3, bool& flag);
  // -- tags contains list of names separated by commas. It 
  //    will replace eventually the previous routine.
  xercesc::DOMNode* getUniqueElementByTagsString_(
      const xercesc::DOMNode* node1, const std::string& tags, bool& flag);

  // -- extract existing atrribute values
  int GetAttributeValueL_(
      xercesc::DOMElement* elem, const char* attr_name, bool exception = true, int val = 0);
  double GetAttributeValueD_(
      xercesc::DOMElement* elem, const char* attr_name, bool exception = true, double val = 0.0);
  std::string GetAttributeValueS_(
      xercesc::DOMElement* elem, const char* attr_name, bool exception = true, std::string val = "");
  std::vector<double> GetAttributeVector_(
      xercesc::DOMElement* elem, const char* attr_name);

  // -- extract and verify children
  //    the name of identical nodes will be extracted too
  std::vector<xercesc::DOMNode*> getSameChildNodes_(
      xercesc::DOMNode* node, std::string& name, bool& flag, bool exception = false);

  // data streaming/trimming/converting
  // -- times
  double TimeStringToValue_(const std::string& time_value);
  double TimeCharToValue_(const char* time_value);

  // -- coordinates
  std::vector<double> MakeCoordinates_(const std::string& array);

  // -- string modifications
  std::vector<std::string> CharToStrings_(const char* namelist);
  std::string TrimString_(char* tmp);

  // -- vector parsing
  int GetPosition_(const std::vector<std::string>& names, const std::string& name);

  // error messages
  void ThrowErrorIllformed_(
      const std::string& section, const std::string& type, const std::string& ill_formed);
  void ThrowErrorIllformed_(
      const std::string& section, const std::string& type, const std::string& ill_formed, const std::string& options);
  void ThrowErrorMissattr_(
      const std::string& section, const std::string& type, const std::string& missing, const std::string& name);

 protected:
  // variois constants defined by the users
  std::map<std::string, std::string> constants_; 

  xercesc::XercesDOMParser* parser;
  xercesc::DOMDocument* doc_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
