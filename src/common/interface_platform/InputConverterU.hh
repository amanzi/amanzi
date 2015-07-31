/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_INPUT_CONVERTER_UNSTRUCTURED_HH_
#define AMANZI_INPUT_CONVERTER_UNSTRUCTURED_HH_

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
#include "InputConverter.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziInput {

typedef std::map<std::string, Teuchos::RCP<Teuchos::ParameterList> > PK;
typedef std::map<std::string, std::vector<std::string> > Tree;

class InputConverterU : public InputConverter {
 public:
  InputConverterU() {};
  ~InputConverterU() {};

  // main members
  Teuchos::ParameterList Translate();

 private:
  Teuchos::ParameterList TranslateMesh_();
  Teuchos::ParameterList TranslateRegions_();

 private:
  int dim_;
  Tree tree_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
