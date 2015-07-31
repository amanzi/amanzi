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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

// Amanzi's
#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziInput {

/*
  A simple wrapper that delegates transcode to Xercecs and deallocates 
  automatically memory.
*/
/*
class XChar {
 public:
  XChar();
  XMLCh* transcode(const char* name) {
    str_ = xercesc::XMLString::transcode(name);
    return str_;
  } 
  ~XChar() { xercesc::XMLString::release(&str_); }

 private:
  XMLCh* str_;
};
*/

class InputConverter {
 public:
  InputConverter() : vo_(NULL) {
    xercesc::XMLPlatformUtils::Initialize();
  }

  ~InputConverter() {
    if (vo_ != NULL) delete vo_;
    xercesc::XMLPlatformUtils::Terminate();
  }

  // main member: creates xerces document using the file name
  void Init(const std::string& xmlfilename);

 protected:
  // verbosity XML
  Teuchos::ParameterList GetVerbosity_();

  // data streaming/trimming/converting
  Teuchos::Array<double> MakeCoordinates_(char* char_array);
  std::string TrimString_(char* tmp);

  // error messages
  void ThrowErrorIllformed_(std::string section, std::string element_type, std::string ill_formed);
  void ThrowErrorIllformed_(std::string section, std::string element_type, std::string ill_formed, std::string options);
  void ThrowErrorMissattr_(std::string section, std::string att_elem_type, std::string missing, std::string elem_name);

 protected:
  xercesc::DOMDocument* doc_;

  Teuchos::ParameterList verb_list_;
  VerboseObject* vo_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
