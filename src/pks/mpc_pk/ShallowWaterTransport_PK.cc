/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Sequential coupling of shallow water and solute transport.
*/

#include <string>

#include "PK_Utils.hh"
#include "Transport_PK.hh"

#include "ShallowWaterTransport_PK.hh"

namespace Amanzi {

/* ******************************************************************
* Constructor
****************************************************************** */
ShallowWaterTransport_PK::ShallowWaterTransport_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : PK_MPCWeak(pk_tree, global_list, S, soln), cfl_(1.0), failed_steps_(0)
{
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = my_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("ShallowWaterTransport", vlist));
}


/* ******************************************************************
* Calculate the min of sub-PKs timesteps and limit it.
****************************************************************** */
double
ShallowWaterTransport_PK::get_dt()
{
  return cfl_ * PK_MPCWeak::get_dt();
}


/* ******************************************************************
* Setup of PK
****************************************************************** */
void
ShallowWaterTransport_PK::Setup()
{
  PK_MPCWeak::Setup();
}


/* ******************************************************************
* Extended treatment of time step in transport PK.
****************************************************************** */
bool
ShallowWaterTransport_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail(false);
  std::string domain = my_list_->get<std::string>("domain name");

  Key velocity_key = Keys::getKey(domain, "velocity");
  Key discharge_key = Keys::getKey(domain, "discharge");
  Key ponded_depth_key = Keys::getKey(domain, "ponded_depth");
  Key prev_ponded_depth_key = Keys::getKey(domain, "prev_ponded_depth");

  StateArchive archive(S_, vo_);
  archive.Add({ ponded_depth_key, prev_ponded_depth_key },
              { discharge_key },
              { ponded_depth_key, velocity_key },
              Tags::DEFAULT,
              name());

  double dt0 = t_new - t_old;
  fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  if (fail) {
    cfl_ /= 2.0;
    return fail;
  }

  // Typically transport has bigger time step; however, numerics may
  // brings challenges. FIXME
  double dt1 = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(sub_pks_[1])->StableTimeStep(1);
  if (dt1 < dt0 / 100) {
    archive.Restore("");

    cfl_ /= 2.0;
    failed_steps_++;

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "total failed steps: " << failed_steps_ << std::endl;
    }
    return true;
  }

  cfl_ = 1.0;
  fail = sub_pks_[1]->AdvanceStep(t_old, t_new, reinit);
  return fail;
}

} // namespace Amanzi
