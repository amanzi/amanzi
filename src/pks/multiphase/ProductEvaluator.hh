/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Molar mobility is defined as the product of molar density and mobility.
*/

#ifndef AMANZI_MULTIPHASE_MOLAR_MOBILITY_EVALUATOR_HH_
#define AMANZI_MULTIPHASE_MOLAR_MOBILITY_EVALUATOR_HH_

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

 // interface to multiphase base class
 virtual void set_subvector(int ifield, int n) { field_n_[ifield] = n; } override;

 private:
  std::vector<int> power_;
  std::vector<int> field_n_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
