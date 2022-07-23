#ifndef AMANZI_INPUT_ANALYSIS_HH_
#define AMANZI_INPUT_ANALYSIS_HH_

#include <string>

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "InputAnalysis.hh"
#include "Mesh.hh"
#include "VerboseObject.hh"

namespace Amanzi {

class InputAnalysis {
 public:
  InputAnalysis(Teuchos::RCP<const AmanziMesh::Mesh> mesh, const std::string& domain) 
    : mesh_(mesh), domain_(domain), vo_(NULL) {};
  ~InputAnalysis() {
    if (vo_ != NULL) delete vo_;
  };

  // main members
  void Init(Teuchos::ParameterList& plist);
  void RegionAnalysis();
  void OutputBCs();

  // supporting members
  template<class Iterator> 
  Iterator SelectUniqueEntries(Iterator first, Iterator last);

 private:
  Teuchos::ParameterList* plist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string domain_;
  VerboseObject* vo_;
};

}  // namespace Amanzi

#endif
