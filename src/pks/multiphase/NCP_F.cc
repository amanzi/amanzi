/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Field evaluator for noninear complimentary problem, function F.
*/

#include "NCP_F.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
NCP_F::NCP_F(Teuchos::ParameterList& plist) : MultiphaseBaseEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("my key");
  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");
  dependencies_.insert(std::string(saturation_liquid_key_));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
NCP_F::NCP_F(const NCP_F& other) : MultiphaseBaseEvaluator(other) {};


Teuchos::RCP<Evaluator> NCP_F::Clone() const {
  return Teuchos::rcp(new NCP_F(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void NCP_F::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& sl = *S.GetFieldData(saturation_liquid_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = 1.0 - sl[0][c];
  }      
}


/* ******************************************************************
* Required member: field derivative calculation.
****************************************************************** */
void NCP_F::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
}

}  // namespace Multiphase
}  // namespace Amanzi


/*
  } else if (ncp_ == "Fischer-Burmeister") {
    for (int c = 0; c < ncells_owned_; ++c) {
      double x_sum(0.0);
      for (int i = 0; i < num_primary_; ++i) x_sum += uci[i][c];

      double a = 1.0 - usi[0][c];
      double b = x_sum;
      fci[n][c] = pow(a * a + b * b, 0.5) - (a + b);
    }
  }
*/
