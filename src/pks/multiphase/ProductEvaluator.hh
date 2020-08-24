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

#ifndef AMANZI_MULTIPHASE_PRODUCT_EVALUATOR_HH_
#define AMANZI_MULTIPHASE_PRODUCT_EVALUATOR_HH_

// Multiphase
#include "MultiphaseBaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class ProductEvaluator : public MultiphaseBaseEvaluator {
 public:
  ProductEvaluator(Teuchos::ParameterList& plist);
  ProductEvaluator(const ProductEvaluator& other);

  // inteface functions to FieldEvaluator
  virtual Teuchos::RCP<FieldEvaluator> Clone() const override;

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) override;

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
      Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

  // extended interface
  virtual void set_subvector(int ifield, int n, double kH) override { field_n_[ifield] = n; kH_ = kH; } 

 private:
  std::vector<int> powers_;
  std::vector<int> field_n_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
