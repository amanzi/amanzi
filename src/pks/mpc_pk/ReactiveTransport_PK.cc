/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Process kernel for coupling of Transport PK and Chemistry PK.

  Note, this MPC uses `"PK order`" only as a list of sub-PKs -- it ignores the
  actual order, because Transport's Initialize() (and maybe Setup()) methods
  require that _chemistry_ goes first, while this MPC requires the Transport's
  Advance() goes first (so that chemical species are equilibrated at the end of
  the step).  Since this order is so carefully controled by this MPC, we ignore
  the user's order.
*/

#include "ReactiveTransport_PK.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
ReactiveTransport_PK::ReactiveTransport_PK(Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                           const Teuchos::RCP<State>& S,
                                           const Teuchos::RCP<TreeVector>& soln)
  : Amanzi::PK_MPCAdditive<PK>(pk_tree, global_list, S, soln)
{
  int i_transport(-1);
  for (int i = 0; i < 2; ++i) {
    transport_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(sub_pks_[i]);
    if (transport_pk_ != Teuchos::null) {
      i_transport = i;
      break;
    }
  }
  AMANZI_ASSERT(i_transport >= 0);

  int i_chem = 1 - i_transport;
  chemistry_pk_ = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(sub_pks_[i_chem]);
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);

  // communicate chemistry engine to transport.
  transport_pk_->SetupChemistry(chemistry_pk_);

  // force sub_pks_ to have chemistry, then transport.  This will set the order
  // for all other methods (other than Advance, which controls the order
  // explicitly).
  sub_pks_ = { chemistry_pk_, transport_pk_ };

  // This algorithm is operator split -- we insert a "operator_split" tag so
  // that transport advances from CURRENT --> "operator_split", and chemistry
  // advances from "operator_split" --> NEXT.
  getSubPKPlist_(0)->set("concentration tag next", "operator_split");
  getSubPKPlist_(1)->set("concentration tag current", "operator_split");
  getSubPKPlist_(1)->set("primary variable password", "state");
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
ReactiveTransport_PK::get_dt()
{
  dTtran_ = transport_pk_->get_dt();
  dTchem_ = chemistry_pk_->get_dt();
  if (dTtran_ > dTchem_) dTtran_ = dTchem_;
  return dTchem_;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void
ReactiveTransport_PK::set_dt(double dt)
{
  dTtran_ = dt;
  dTchem_ = dt;
  chemistry_pk_->set_dt(dTchem_);
  transport_pk_->set_dt(dTtran_);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
ReactiveTransport_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // NOTE: this is implemented to force the _reverse_ order of Initialize() and
  // Setup() -- transport is done first to ensure that species are equilibrated
  // at the end of the timestep.
  bool fail = false;

  // First we do a transport step.
  bool pk_fail = sub_pks_[1]->AdvanceStep(t_old, t_new, reinit);
  if (pk_fail) return pk_fail;

  // Second, we do a chemistry step using a copy of the tcc vector
  pk_fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  return pk_fail;
};


// -----------------------------------------------------------------------------
// Only one needs to commit due to the shared primary unknown
//
// ETC: this seems dangerous -- what if transport needs to commit something
// else?  Is there a downside to commiting the primary variable twice?
// -----------------------------------------------------------------------------
void
ReactiveTransport_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  chemistry_pk_->CommitStep(t_old, t_new, tag);
}

} // namespace Amanzi
