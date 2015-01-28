/*
  This is the operators component of the Amanzi code.

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete source operator.
*/

#ifndef AMANZI_OPERATOR_SOURCE_HH_
#define AMANZI_OPERATOR_SOURCE_HH_

#include "Operator.hh"

namespace Amanzi {
namespace Operators {

class OperatorSource : public Operator {
 public:
  OperatorSource() {};
  OperatorSource(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) : Operator(cvs, dummy) {};
  OperatorSource(const Operator& op) : Operator(op) {};
  ~OperatorSource() {};

  // main members
  void UpdateMatrices(const CompositeVector& src);
};

}  // namespace Operators
}  // namespace Amanzi

#endif

