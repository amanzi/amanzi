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
  InputConverterU() :
      vo_(NULL),
      flow_single_phase(false),
      compressibility_(false) {};
  ~InputConverterU() { if (vo_ != NULL) delete vo_; }

  // main members
  Teuchos::ParameterList Translate();

 private:
  void ParseSolutes_();

  Teuchos::ParameterList TranslateMesh_();
  Teuchos::ParameterList TranslateRegions_();
  Teuchos::ParameterList TranslateOutput_();
  Teuchos::ParameterList TranslatePreconditioners_();
  Teuchos::ParameterList TranslateTrilinosML_();
  Teuchos::ParameterList TranslateHypreAMG_();
  Teuchos::ParameterList TranslateBILU_();
  Teuchos::ParameterList TranslateSolvers_();
  Teuchos::ParameterList TranslateState_();
  Teuchos::ParameterList TranslateMaterialsPartition_();
  Teuchos::ParameterList TranslateCycleDriver_();

  void ProcessMacros_(const std::string& prefix, char* text_content,
                      Teuchos::ParameterList& mPL, Teuchos::ParameterList& outPL);

  Teuchos::ParameterList GetVerbosity_();

 private:
  int dim_;
  Tree tree_;
  Tree phases_;

  bool flow_single_phase;
  bool compressibility_;
  std::vector<std::string> comp_names_all_;

  Teuchos::ParameterList verb_list_;
  VerboseObject* vo_;
};

}  // namespace AmanziInput
}  // namespace Amanzi

#endif
