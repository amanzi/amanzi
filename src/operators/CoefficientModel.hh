/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)

  Collection of classes describing coefficient models.
*/

#ifndef AMANZI_OPERATOR_COEFFICIENT_MODEL_HH_
#define AMANZI_OPERATOR_COEFFICIENT_MODEL_HH_

#include <string>

#include "Teuchos_RCP.hpp"

namespace Amanzi {
namespace Operators {

class CoefficientModel {
 public:
  CoefficientModel() {}; 
  ~CoefficientModel() {}; 

  virtual std::string name() const = 0;
};


class CoefficientConstant : public CoefficientModel {
 public:
  CoefficientConstant(const std::shared_ptr<std::vector<WhetStone::Tensor> >& coef)
    : coef_(coef) {}; 
  ~CoefficientConstant() {}; 

  virtual std::string name() const override { return "constant"; }

 public:
  std::shared_ptr<std::vector<WhetStone::Tensor> > coef_;
};


class CoefficientMatrixPolynomial : public CoefficientModel {
 public:
  CoefficientMatrixPolynomial(const std::shared_ptr<std::vector<WhetStone::MatrixPolynomial> >& coef)
    : coef_(coef) {};
  ~CoefficientMatrixPolynomial() {}; 

  virtual std::string name() const override { return "matrix polynomial"; }

 public:
  std::shared_ptr<std::vector<WhetStone::MatrixPolynomial> > coef_;
};


class CoefficientFunction : public CoefficientModel {
 public:
  CoefficientFunction(const std::shared_ptr<std::vector<WhetStone::WhetStoneFunction*> >& coef)
    : coef_(coef) {};
  ~CoefficientFunction() {}; 

  virtual std::string name() const override { return "function"; }

 public:
  std::shared_ptr<const std::vector<WhetStone::WhetStoneFunction*> > coef_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


