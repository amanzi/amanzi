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
  the powers ai are either 1 or -1. If fi is a multivector, one
  of its entries is used in the product, see field_n_ variable.
*/

#ifndef AMANZI_MULTIPHASE_PRODUCT_EVALUATOR_HH_
#define AMANZI_MULTIPHASE_PRODUCT_EVALUATOR_HH_

#include "Factory.hh"

// Multiphase
#include "MultiphaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class ProductEvaluator : public MultiphaseEvaluator {
 public:
  ProductEvaluator(Teuchos::ParameterList& plist);
  ProductEvaluator(const ProductEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  // extended interface
  virtual void set_subvector(int ifield, int n, const std::string& name,
                             Teuchos::ParameterList& plist) override
  {
    field_n_[ifield] = n;
    AmanziEOS::VaporLiquidFactory factory(plist);
    vapor_liquid_ = factory.Create(name);
  }

 private:
  std::vector<int> powers_;
  std::vector<int> field_n_;

  static Utils::RegisteredFactory<Evaluator, ProductEvaluator> fac_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
