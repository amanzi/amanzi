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
  // why put chemistry before transport if you intend to solve transport before
  // chemistry? --ETC
  transport_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(sub_pks_[1]);
  AMANZI_ASSERT(transport_pk_ != Teuchos::null);

  chemistry_pk_ = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(sub_pks_[0]);
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);

  // communicate chemistry engine to transport.
  transport_pk_->SetupChemistry(chemistry_pk_);

  // tell chemistry it is operator split, which means it will use TCC next
  // instead of current.  Also tell it to use Amanzi generic passwd
  getSubPKPlist_(0)->set("operator split", true);
  getSubPKPlist_(0)->set("primary variable password", "state");
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
  // NOTE: this need not be implemented except that the PKs are reversed.  Is
  // there a good reason that chemistry's Setup() is called before transport's
  // Setup(), but transport's Advance() is called before chemistry's Advance()?
  // --ETC
  bool fail = false;

  // First we do a transport step.
  bool pk_fail = transport_pk_->AdvanceStep(t_old, t_new, reinit);
  if (pk_fail) return pk_fail;

  // Second, we do a chemistry step using a copy of the tcc vector
  pk_fail = chemistry_pk_->AdvanceStep(t_old, t_new, reinit);
  return pk_fail;
};

} // namespace Amanzi
