/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

  Collection of classes describing coefficient models.
*/

#ifndef AMANZI_OPERATOR_COEFFICIENT_MODEL_HH_
#define AMANZI_OPERATOR_COEFFICIENT_MODEL_HH_

#include <string>
#include <typeinfo>

#include "Teuchos_RCP.hpp"

namespace Amanzi {
namespace Operators {

template <typename T>
class CoefficientModel {
 public:
  CoefficientModel(){};
  CoefficientModel(const std::shared_ptr<std::vector<T>>& coef) : coef_(coef){};
  ~CoefficientModel(){};

  std::string name() { return typeid(T).name(); }
  T get_coef(int c) { return (*coef_)[c]; }

 public:
  std::shared_ptr<std::vector<T>> coef_;
};


/* ******************************************************************
* Specialization
****************************************************************** */
template <>
inline WhetStone::Tensor
CoefficientModel<WhetStone::Tensor>::get_coef(int c)
{
  WhetStone::Tensor Kc(2, 1);
  Kc(0, 0) = 1.0;
  return (coef_.get()) ? (*coef_)[c] : Kc;
}

} // namespace Operators
} // namespace Amanzi

#endif
