/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Adds one function to the second variable evaluator.
*/

#ifndef AMANZI_MULTIPHASE_BASE_EVALUATOR_HH_
#define AMANZI_MULTIPHASE_BASE_EVALUATOR_HH_

// Amanzi
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Multiphase {

class MultiphaseBaseEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  MultiphaseBaseEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist), n_(0), kH_(1.0) {
};

  // inteface functions to FieldEvaluator
  MultiphaseBaseEvaluator(const MultiphaseBaseEvaluator& other)
    : SecondaryVariableFieldEvaluator(other) { n_ = other.n_; kH_ = other.kH_; }

  // added interface
  virtual void set_subvector(int ifield, int n, double kH) { n_ = n; kH_ = kH; } 

 protected:
  int n_;
  double kH_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
