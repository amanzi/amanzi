/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Abstract operator. 
*/

#ifndef AMANZI_OPERATOR_ABSTRACT_HH_
#define AMANZI_OPERATOR_ABSTRACT_HH_

#include <string>

#include "PDE_Helper.hh"

namespace Amanzi {
namespace Operators {

class Abstract : public PDE_Helper {
 public:
  Abstract(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op) :
      PDE_Helper(global_op) {
    Init_(plist);
  }

  Abstract(Teuchos::ParameterList& plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      PDE_Helper(mesh) {
    global_op_ = Teuchos::null;
    Init_(plist);
  }

  // main members 
  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K) { K_ = K; }
  void UpdateMatrices();
  
 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;

  Schema global_schema_row_, global_schema_col_;
  Schema local_schema_col_, local_schema_row_;

 private:
  void Init_(Teuchos::ParameterList& plist);

 private:
  std::string polytope_, method_, matrix_;
  bool hybridize_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

