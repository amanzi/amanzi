/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Secondary variable field evaluator computes product of fields 
  or inverse of fields:

    eval = f1 * f2 * ... * fn) / (g1 * g2 * ... * gm)
*/

#include "FieldEvaluator.hh"
#include "MultiplicativeReciprocalEvaluator.hh"

namespace Amanzi {

/* ******************************************************************
* Two constructors.
****************************************************************** */
MultiplicativeReciprocalEvaluator::MultiplicativeReciprocalEvaluator(Teuchos::ParameterList& plist)
  : SecondaryVariableFieldEvaluator(plist)
{
  my_key_ = plist.get<std::string>("my key");
  Key domain = Keys::getDomain(my_key_);

  if (plist_.isParameter("evaluator dependencies")) {
    Errors::Message msg;
    msg << "MultiplicativeReciprocalEvaluator: \"" << my_key_ 
        << "\" must have separate (optional) lists for multiplicative and reciprocal dependencies.";
    Exceptions::amanzi_throw(msg);
  }

  if (plist_.isParameter("multiplicative dependencies")) {
    // since dependensies is a map, we need separate maps for numerator and denominator
    const auto& names = plist_.get<Teuchos::Array<std::string> >("multiplicative dependencies");
    for (const auto& name : names) {
      Key full_name = Keys::getKey(domain, name);
      dependencies_.insert(full_name);
      list0_.push_back(full_name);
    }
  }

  if (plist_.isParameter("reciprocal dependencies")) {
    const auto& names = plist_.get<Teuchos::Array<std::string> >("reciprocal dependencies");
    for (const auto& name : names) {
      Key full_name = Keys::getKey(domain, name);
      dependencies_.insert(full_name);
      list1_.push_back(full_name);
    }
  }

  if (list0_.size() + list1_.size() == 0) {
    Errors::Message msg;
    msg << "MultiplicativeReciprocalEvaluator for: \"" << my_key_ << "\" has no dependencies.";
    Exceptions::amanzi_throw(msg);
  }

  coef_ = plist_.get<double>("coefficient", 1.0);
  enforce_positivity_ = plist_.get<bool>("enforce positivity", false);
}


MultiplicativeReciprocalEvaluator::MultiplicativeReciprocalEvaluator(const MultiplicativeReciprocalEvaluator& other)
  : SecondaryVariableFieldEvaluator(other) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> MultiplicativeReciprocalEvaluator::Clone() const {
  return Teuchos::rcp(new MultiplicativeReciprocalEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void MultiplicativeReciprocalEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    auto& result_c = *result->ViewComponent(*comp);
    int ndofs = result_c.MyLength();

    result_c.PutScalar(coef_);

    for (auto it = list0_.begin(); it != list0_.end(); ++it) {
      const auto& factor_c = *S->GetFieldData(*it)->ViewComponent(*comp);
      for (int c = 0; c != ndofs; ++c) result_c[0][c] *= factor_c[0][c];
    }

    for (auto it = list1_.begin(); it != list1_.end(); ++it) {
      const auto& factor_c = *S->GetFieldData(*it)->ViewComponent(*comp);
      for (int c = 0; c != ndofs; ++c) result_c[0][c] /= factor_c[0][c];
    }

    if (enforce_positivity_) {
      for (int c = 0; c != ndofs; ++c) {
        result_c[0][c] = std::max(result_c[0][c], 0.0);
      }
    }
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void MultiplicativeReciprocalEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    auto& result_c = *result->ViewComponent(*comp);
    int ncells = result_c.MyLength();

    result_c.PutScalar(1.0);
    for (auto it = list0_.begin(); it != list0_.end(); ++it) {
      const auto& factor_c = *S->GetFieldData(*it)->ViewComponent(*comp);
      if (*it != wrt_key)
        for (int c = 0; c != ncells; ++c) result_c[0][c] *= factor_c[0][c];
    }

    for (auto it = list1_.begin(); it != list1_.end(); ++it) {
      const auto& factor_c = *S->GetFieldData(*it)->ViewComponent(*comp);
      if (*it == wrt_key)
        for (int c = 0; c != ncells; ++c) result_c[0][c] /= -factor_c[0][c] * factor_c[0][c];
      else 
        for (int c = 0; c != ncells; ++c) result_c[0][c] /= factor_c[0][c];
    }
  }
}

}  // namespace Amanzi
