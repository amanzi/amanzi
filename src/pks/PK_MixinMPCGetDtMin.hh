/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! The mixin with default implementation for a weak MPC's AdvanceStep method.

/*!

Solely implements the AdvanceStep() method.

*/

#ifndef AMANZI_PK_MIXIN_MPC_GET_DT_MIN_HH_
#define AMANZI_PK_MIXIN_MPC_GET_DT_MIN_HH_

#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"

namespace Amanzi {

template <class Base_t>
class PK_MixinMPCGetDtMin : public Base_t {
 public:
  using Base_t::Base_t;

  // Get the dt as the min of all children PKs.
  double get_dt();
};

// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class Base_t>
double
PK_MixinMPCGetDtMin<Base_t>::get_dt()
{
  double dt = 1.e80;
  for (auto& pk : this->sub_pks_) dt = std::min(dt, pk->get_dt());
  return dt;
}

} // namespace Amanzi

#endif
