/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_AUDIT_HH_
#define AMANZI_OPERATOR_AUDIT_HH_

#include "Operator.hh"


namespace Amanzi {
namespace Operators {

class OperatorAudit : public Operator {
 public:
  OperatorAudit() {};
  OperatorAudit(const Operator& op) : Operator(op) {};
  ~OperatorAudit() {};

  // main members
  int CheckSpectralBounds(int schema);
  int CheckMatrixSymmetry();
  int CheckMatrixCoercivity();

 private:
  void OrderByIncrease(int n, double* mem);
};

}  // namespace Operators
}  // namespace Amanzi

#endif
