/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATORS_BC_HH_
#define AMANZI_OPERATORS_BC_HH_

#include <vector>

namespace Amanzi {
namespace Operators {

class BCs {
 public:
  BCs() {};
  BCs(std::vector<int>& bc_model, std::vector<double>& bc_value) { Init(bc_model, bc_value); }
  ~BCs() {};

  // main members
  void Init(std::vector<int>& bc_model, std::vector<double>& bc_value) {
    bc_model_ = &bc_model; 
    bc_value_ = &bc_value; 
  }
 
  // access
  const std::vector<int>& bc_model() { return *bc_model_; }
  const std::vector<double>& bc_value() { return *bc_value_; }

 private:
  std::vector<int>* bc_model_;
  std::vector<double>* bc_value_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


