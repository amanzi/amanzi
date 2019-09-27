/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)

  Interface class that binds discretization class with coefficient model.
*/

#ifndef AMANZI_OPERATOR_INTERFACE_WHETSTONE_HH_
#define AMANZI_OPERATOR_INTERFACE_WHETSTONE_HH_

#include "Teuchos_RCP.hpp"

#include "DenseMatrix.hh"

#include "CoefficientModel.hh"

namespace Amanzi {
namespace Operators {

class InterfaceWhetStone {
 public:
  InterfaceWhetStone() {};
  virtual ~InterfaceWhetStone() {};

  virtual void StiffnessMatrix(int c, WhetStone::DenseMatrix& Acell) = 0;
  virtual void FaceMatrixJump(int f, int c1, int c2, WhetStone::DenseMatrix& Aface) = 0;
};


template<class T, class U>
class InterfaceWhetStoneImpl : public InterfaceWhetStone {
 public:
  InterfaceWhetStoneImpl(const Teuchos::RCP<T>& mfd, const std::shared_ptr<U>& coef)
    : mfd_(mfd), coef_(coef) {};

  virtual void StiffnessMatrix(int c, WhetStone::DenseMatrix& Acell) override {
    mfd_->StiffnessMatrix(c, (*coef_->coef_)[c], Acell);
  }

  virtual void FaceMatrixJump(int f, int c1, int c2, WhetStone::DenseMatrix& Aface) override {
    mfd_->FaceMatrixJump(f, (*coef_->coef_)[c1], (*coef_->coef_)[c2], Aface);
  }

 private:
  Teuchos::RCP<T> mfd_;
  std::shared_ptr<U> coef_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


