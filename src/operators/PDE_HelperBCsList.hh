/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_OPERATOR_BCS_LIST_HH_
#define AMANZI_OPERATOR_BCS_LIST_HH_

#include "Teuchos_RCP.hpp"

#include "BCs.hh"

namespace Amanzi {
namespace Operators {

class PDE_HelperBCsList {
 public:
  PDE_HelperBCsList() {};
  virtual ~PDE_HelperBCsList() {};

  // boundary conditions (BC) require information on test and
  // trial spaces. For a single PDE, these BCs could be the same.
  virtual void SetBCs(const Teuchos::RCP<const BCs>& bc_trial,
                      const Teuchos::RCP<const BCs>& bc_test)
  {
    bcs_trial_.clear();
    bcs_test_.clear();

    bcs_trial_.push_back(bc_trial);
    bcs_test_.push_back(bc_test);
  }

  virtual void AddBCs(const Teuchos::RCP<const BCs>& bc_trial,
                      const Teuchos::RCP<const BCs>& bc_test)
  {
    bcs_trial_.push_back(bc_trial);
    bcs_test_.push_back(bc_test);
  }

  // const access
  const std::vector<Teuchos::RCP<const BCs>>& get_bcs_trial() { return bcs_trial_; }
  const std::vector<Teuchos::RCP<const BCs>>& get_bcs_test() { return bcs_test_; }

 protected:
  std::vector<Teuchos::RCP<const BCs>> bcs_trial_, bcs_test_;
};

} // namespace Operators
} // namespace Amanzi

#endif
