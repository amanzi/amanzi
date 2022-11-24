/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)

  Interface class that hides details of the actual coefficient model
  using a few template classes.
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
  InterfaceWhetStone(){};
  virtual ~InterfaceWhetStone(){};

  virtual void MassMatrix(int c, WhetStone::DenseMatrix& Acell){};
  virtual void MassMatrixInverse(int c, WhetStone::DenseMatrix& Acell){};
  virtual void StiffnessMatrix(int c, WhetStone::DenseMatrix& Acell){};
  virtual void FaceMatrixJump(int f, int c1, int c2, WhetStone::DenseMatrix& Aface){};
};


template <class T, class U>
class InterfaceWhetStoneDG : public InterfaceWhetStone {
 public:
  InterfaceWhetStoneDG(const Teuchos::RCP<T>& dg, const std::shared_ptr<U>& coef)
    : dg_(dg), coef_(coef){};

  virtual void StiffnessMatrix(int c, WhetStone::DenseMatrix& Acell) override
  {
    dg_->StiffnessMatrix(c, coef_->get_coef(c), Acell);
  }

  virtual void FaceMatrixJump(int f, int c1, int c2, WhetStone::DenseMatrix& Aface) override
  {
    dg_->FaceMatrixJump(f, coef_->get_coef(c1), coef_->get_coef(c2), Aface);
  }

 private:
  Teuchos::RCP<T> dg_;
  std::shared_ptr<U> coef_;
};


template <class T, class U>
class InterfaceWhetStoneMFD : public InterfaceWhetStone {
 public:
  InterfaceWhetStoneMFD(const Teuchos::RCP<T>& mfd, const std::shared_ptr<U>& coef)
    : mfd_(mfd), coef_(coef){};

  virtual void MassMatrix(int c, WhetStone::DenseMatrix& Acell) override
  {
    mfd_->MassMatrix(c, coef_->get_coef(c), Acell);
  }

  virtual void MassMatrixInverse(int c, WhetStone::DenseMatrix& Acell) override
  {
    mfd_->MassMatrixInverse(c, coef_->get_coef(c), Acell);
  }

  virtual void StiffnessMatrix(int c, WhetStone::DenseMatrix& Acell) override
  {
    mfd_->StiffnessMatrix(c, coef_->get_coef(c), Acell);
  }

 private:
  Teuchos::RCP<T> mfd_;
  std::shared_ptr<U> coef_;
};

} // namespace Operators
} // namespace Amanzi

#endif
