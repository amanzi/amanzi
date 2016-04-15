/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_ELECTROMAGNETICS_MHD_HH_
#define AMANZI_OPERATOR_ELECTROMAGNETICS_MHD_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "exceptions.hh"
#include "Tensor.hh"

#include "BCs.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "OperatorElectromagnetics.hh"

namespace Amanzi {
namespace Operators {

class OperatorElectromagneticsMHD : public OperatorElectromagnetics {
 public:
  OperatorElectromagneticsMHD(const Teuchos::RCP<Operator>& global_op)
    : OperatorElectromagnetics(global_op)
  {};

  OperatorElectromagneticsMHD(Teuchos::ParameterList& plist,
                              const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : OperatorElectromagnetics(plist, mesh)
  {
    operator_type_ = OPERATOR_ELECTROMAGNETICS_MHD;
    InitElectromagneticsMHD_();
  }

  // main virtual members
  // -- creation of an operator
  virtual void UpdateMatrices();

  // -- before solving the problem
  virtual void ModifyMatrices(CompositeVector& E, CompositeVector& B);

  // -- after solving the problem
  virtual void ModifyFields(CompositeVector& E, CompositeVector& B);

 private:
  void InitElectromagneticsMHD_();

 protected:
  std::vector<WhetStone::DenseMatrix> mass_op_;
  std::vector<WhetStone::DenseMatrix> curl_op_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


