/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#pragma once
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {

class Function; // forward declaration

class FunctionFactory {
 public:
  FunctionFactory() {}
  ~FunctionFactory() {}

  // parse the list and create a function
  std::unique_ptr<Function> Create(Teuchos::ParameterList& plist) const;

  // create a function of a given type
  std::unique_ptr<Function>
  Create(const std::string& func_type, Teuchos::ParameterList& func_params) const;

 private:
  std::unique_ptr<Function> create_constant(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_tabular(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_smooth_step(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_polynomial(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_monomial(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_linear(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_separable(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_additive(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_multiplicative(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_composition(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_static_head(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_standard_math(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_bilinear(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_distance(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_squaredistance(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_bilinear_and_time(Teuchos::ParameterList&) const;
  std::unique_ptr<Function> create_exprtk(Teuchos::ParameterList&) const;
};


} // namespace Amanzi
