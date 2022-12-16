/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Multiphase PK

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
ProductEvaluator::ProductEvaluator(Teuchos::ParameterList& plist) : MultiphaseBaseEvaluator(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }

  if (!(plist.isParameter("dependencies") && plist.isParameter("powers"))) {
    Errors::Message msg("Product evaluator requires \"dependencies\" and \"powers\"");
    Exceptions::amanzi_throw(msg);
  }
  powers_ = plist_.get<Teuchos::Array<int>>("powers").toVector();
  field_n_.resize(powers_.size(), 0);
}


ProductEvaluator::ProductEvaluator(const ProductEvaluator& other)
  : MultiphaseBaseEvaluator(other), field_n_(other.field_n_){};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
ProductEvaluator::Clone() const
{
  return Teuchos::rcp(new ProductEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
ProductEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = result_c.MyLength();

  int n(0);
  result_c.PutScalar(1.0);
  for (auto it = dependencies_.begin(); it != dependencies_.end(); ++it) {
    int m = field_n_[n];

    const auto& factor_c = *S.Get<CompositeVector>(it->first).ViewComponent("cell");
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
void
ProductEvaluator::EvaluatePartialDerivative_(const State& S,
                                             const Key& wrt_key,
                                             const Tag& wrt_tag,
                                             const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = result_c.MyLength();

  int n(0);
  result_c.PutScalar(1.0);
  for (auto it = dependencies_.begin(); it != dependencies_.end(); ++it) {
    int m = field_n_[n];

    const auto& factor_c = *S.Get<CompositeVector>(it->first).ViewComponent("cell");
    if (it->first == wrt_key && powers_[n] == 1)
      ; // do nothing
    else if (it->first == wrt_key && powers_[n] == -1)
      for (int c = 0; c != ncells; ++c) result_c[0][c] /= -factor_c[m][c] * factor_c[m][c];
    else if (powers_[n] == 1)
      for (int c = 0; c != ncells; ++c) result_c[0][c] *= factor_c[m][c];
    else if (powers_[n] == -1)
      for (int c = 0; c != ncells; ++c) result_c[0][c] /= factor_c[m][c];

    n++;
  }
}

} // namespace Multiphase
} // namespace Amanzi
