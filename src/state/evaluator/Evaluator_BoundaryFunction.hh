/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Evaluates a function to provide Dirichlet data on faces.

/*
This is moving toward a function approach, but is horribly non-functional.
Effectively this just passes all the work the stored class, which includes
boht data and functions.  Eventually the functions will migrate out of the
data structure into this class...
 */


#ifndef AMANZI_EVALUATOR_BOUNDARY_FUNCTION_HH_
#define AMANZI_EVALUATOR_BOUNDARY_FUNCTION_HH_

#include "EvaluatorIndependent.hh"
#include "BoundaryFunctionFactory.hh"

namespace Amanzi {

class Evaluator_BoundaryFunction
  : public EvaluatorIndependent<Functions::BoundaryFunction,
                                Functions::BoundaryFunctionFactory> {
 public:
  explicit Evaluator_BoundaryFunction(Teuchos::ParameterList& plist)
    : EvaluatorIndependent<Functions::BoundaryFunction,
                           Functions::BoundaryFunctionFactory>(plist),
      list_name_(plist.get<std::string>("list name")),
      function_name_(plist.get<std::string>("function name", ""))
  {}

  Evaluator_BoundaryFunction(const Evaluator_BoundaryFunction& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new Evaluator_BoundaryFunction(*this));
  }

  virtual void EnsureCompatibility(State& S) override
  {
    EvaluatorIndependent<
      Functions::BoundaryFunction,
      Functions::BoundaryFunctionFactory>::EnsureCompatibility(S);
    auto& fac =
      S.Require<Functions::BoundaryFunction,
                Functions::BoundaryFunctionFactory>(my_key_, my_tag_, my_key_);
    fac.set_names(list_name_, function_name_);
    fac.set_parameterlist(plist_);
  }


 protected:
  // ---------------------------------------------------------------------------
  // Does the actual work to update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override
  {
    S.GetW<Functions::BoundaryFunction>(my_key_, my_tag_, my_key_)
      .Compute(S.time(my_tag_));
  }

 protected:
  Key domain_;
  Key list_name_;
  Key function_name_;

 private:
  static Utils::RegisteredFactory<Evaluator, Evaluator_BoundaryFunction> fac_;
};

} // namespace Amanzi

#endif
