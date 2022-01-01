/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Product evaluator computes product of multiple fields:
    Product = f1^a1 * f1^a2 * ... * fn^an
  the powers ai are either 1 or -1.
*/

#include "errors.hh"

#include "ProductEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Two constructors.
****************************************************************** */
ProductEvaluator::ProductEvaluator(Teuchos::ParameterList& plist)
  : MultiphaseBaseEvaluator(plist)
{
  my_key_ = plist.get<std::string>("my key");
  if (!(plist.isParameter("evaluator dependencies") &&
        plist.isParameter("powers"))) {
    Errors::Message msg("Product evaluator requires \"evaluator dependencies\" and \"powers\"");
    Exceptions::amanzi_throw(msg);
  }
  powers_ = plist_.get<Teuchos::Array<int> >("powers").toVector();
  field_n_.resize(powers_.size(), 0);
}


ProductEvaluator::ProductEvaluator(const ProductEvaluator& other)
  : MultiphaseBaseEvaluator(other),
    field_n_(other.field_n_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator> ProductEvaluator::Clone() const {
  return Teuchos::rcp(new ProductEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void ProductEvaluator::EvaluateField_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  auto& result_c = *result->ViewComponent("cell");
  int ncells = result_c.MyLength();

  int n(0);
  result_c.PutScalar(1.0);
  for (auto it = dependencies_.begin(); it != dependencies_.end(); ++it) {
    int m = field_n_[n];

    const auto& factor_c = *S->GetFieldData(*it)->ViewComponent("cell");
    if (powers_[n] == 1)
      for (int c = 0; c != ncells; ++c) result_c[0][c] *= factor_c[m][c];
    else if (powers_[n] == -1)
      for (int c = 0; c != ncells; ++c) result_c[0][c] /= factor_c[m][c];

    n++;
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void ProductEvaluator::EvaluateFieldPartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = result_c.MyLength();

  int n(0);
  result_c.PutScalar(1.0);
  for (auto it = dependencies_.begin(); it != dependencies_.end(); ++it) {
    int m = field_n_[n];

    const auto& factor_c = *S->GetFieldData(*it)->ViewComponent("cell");
    if (*it == wrt_key && powers_[n] == 1)
      ;  // do nothing
    else if (*it == wrt_key && powers_[n] == -1)
      for (int c = 0; c != ncells; ++c) result_c[0][c] /= -factor_c[m][c] * factor_c[m][c];
    else if (powers_[n] == 1)
      for (int c = 0; c != ncells; ++c) result_c[0][c] *= factor_c[m][c];
    else if (powers_[n] == -1)
      for (int c = 0; c != ncells; ++c) result_c[0][c] /= factor_c[m][c];

    n++;
  }
}

}  // namespace Flow
}  // namespace Amanzi
