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

typedef std::pair<std::string, Teuchos::RCP<Teuchos::ParameterList> > PK;
typedef std::pair<std::string, std::vector<std::string> > Solutes;

class InputConverter {
 public:
  InputConverter() : vo_(NULL) {};
  ~InputConverter() { if (vo_ != NULL) delete vo_; }

  // main members
  Teuchos::ParameterList Translate(const std::string& xmlfilename);

 private:
  Teuchos::ParameterList TranslateMesh_(xercesc::DOMDocument* doc);

  Teuchos::ParameterList GetVerbosity_(xercesc::DOMDocument* doc);

  Teuchos::Array<double> MakeCoordinates_(char* char_array);
  std::string TrimString_(char* tmp);

  void ThrowErrorIllformed_(std::string section, std::string element_type, std::string ill_formed);

 private:
  int dim_;
  std::vector<Solutes> solutes_;

 protected:
  VerboseObject* vo_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
