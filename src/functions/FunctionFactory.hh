#ifndef AMANZI_FUNCTION_FACTORY_HH_
#define AMANZI_FUNCTION_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {

class Function;  // forward declaration

class FunctionFactory {
 public:
  FunctionFactory() {}
  ~FunctionFactory() {}
  Function* Create(Teuchos::ParameterList&) const;

 private:
  Function* create_constant(Teuchos::ParameterList&) const;
  Function* create_tabular(Teuchos::ParameterList&) const;
  Function* create_smooth_step(Teuchos::ParameterList&) const;
  Function* create_polynomial(Teuchos::ParameterList&) const;
  Function* create_monomial(Teuchos::ParameterList&) const;
  Function* create_linear(Teuchos::ParameterList&) const;
  Function* create_separable(Teuchos::ParameterList&) const;
  Function* create_additive(Teuchos::ParameterList&) const;
  Function* create_multiplicative(Teuchos::ParameterList&) const;
  Function* create_composition(Teuchos::ParameterList&) const;
  Function* create_static_head(Teuchos::ParameterList&) const;
  Function* create_standard_math(Teuchos::ParameterList&) const;
  Function* create_bilinear(Teuchos::ParameterList&) const;
  Function* create_distance(Teuchos::ParameterList&) const;
  Function* create_squaredistance(Teuchos::ParameterList&) const;
  Function* create_bilinear_and_time(Teuchos::ParameterList&) const;
};

} // namespace Amanzi

#endif // AMANZI_FUNCTION_FACTORY_HH_
