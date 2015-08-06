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

// Amanzi's
#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziInput {

class InputConverter {
 public:
  InputConverter() {
    xercesc::XMLPlatformUtils::Initialize();
  }

  ~InputConverter() {
    xercesc::XMLPlatformUtils::Terminate();
  }

  // main member: creates xerces document using the file name
  void Init(const std::string& xmlfilename);

 protected:
  // DOM useful tool
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

  int GetAttributeValueL_(xercesc::DOMElement* elem, const char* attr_name);
  double GetAttributeValueD_(xercesc::DOMElement* elem, const char* attr_name);
  std::vector<double> GetAttributeVector_(xercesc::DOMElement* elem, const char* attr_name);

  // data streaming/trimming/converting
  // -- times
  double GetTimeValue_(std::string time_value);
  double ConvertTimeValue_(char* time_value);

  // -- coordinates
  std::vector<double> MakeCoordinates_(char* char_array);

  // -- string modifications
  std::vector<std::string> CharToStrings_(const char* namelist);
  std::string TrimString_(char* tmp);

  // error messages
  void ThrowErrorIllformed_(
      const std::string& section, const std::string& type, const std::string& ill_formed);
  void ThrowErrorIllformed_(
      const std::string& section, const std::string& type, const std::string& ill_formed, const std::string& options);
  void ThrowErrorMissattr_(
      const std::string& section, const std::string& type, const std::string& missing, const std::string& name);

 protected:
  xercesc::DOMDocument* doc_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
