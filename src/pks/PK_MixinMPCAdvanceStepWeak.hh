/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

//! The mixin with default implementation for a weak MPC's AdvanceStep method.

/*!  

Solely implements the AdvanceStep() method.
  
*/


#ifndef AMANZI_PK_MIXIN_MPC_ADVANCE_WEAK_HH_
#define AMANZI_PK_MIXIN_MPC_ADVANCE_WEAK_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"

#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"

namespace Amanzi {

template <class Base_t>
class PK_MixinMPCAdvanceStepWeak : public Base_t {
 public:
  using Base_t::Base_t;

  // Advance PK from time tag old to time tag new
  bool AdvanceStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                   const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new);

 protected:
  using Base_t::vo_;
  using Base_t::S_;
  
};


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class Base_t>
bool
PK_MixinMPCAdvanceStepWeak<Base_t>::AdvanceStep(
    const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
    const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new)
{
  double t_old = S_->time(tag_old);
  double t_new = S_->time(tag_new);

  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << t_old
               << " t1 = " << t_new << " h = " << t_new - t_old << std::endl
               << "----------------------------------------------------------------" << std::endl;

  int i = 0;
  for (auto& pk : this->sub_pks_) {
    bool fail = pk->AdvanceStep(tag_old, soln_old->SubVector(i),
            tag_new, soln_new->SubVector(i));
    if (fail) return fail;
    ++i;
  }
  return false;
}


}  // namespace Amanzi

#endif
