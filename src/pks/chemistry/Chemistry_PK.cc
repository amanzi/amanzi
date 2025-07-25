/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry PK

  Base class for chemical process kernels.

  Handles timestep control, and generic advancing by cell.
*/

// Chemistry
#include "message.hh"
#include "Reductions.hh"
#include "PK_Helpers.hh"
#include "TimestepControllerFactory.hh"
#include "Chemistry_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

/* ******************************************************************
* Standard PK constructor
****************************************************************** */
Chemistry_PK::Chemistry_PK(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln)
  : PK_Physical_Default(pk_tree, glist, S, soln),
    PK(pk_tree, glist, S, soln),
    number_aqueous_components_(0),
    number_gaseous_components_(0),
    number_mineral_components_(0),
    dt_next_(-1.)
{
  // note, we pass in null to the factory here to make sure there is no error
  // control used, which doesn't make sense for this application.
  timestep_controller_ = createTimestepController<TreeVector>(
    name_, plist_->sublist("timestep controller"), S_, Teuchos::null, Teuchos::null);
};


double
Chemistry_PK::get_dt()
{
  if (dt_next_ < 0.) dt_next_ = timestep_controller_->getInitialTimestep();
  return dt_next_;
}


/* ******************************************************************
* Parser
******************************************************************* */
void
Chemistry_PK::parseParameterList()
{
  // with subfield names, the default width is often insufficient
  if (!plist_->sublist("verbose object").isParameter("debug cell header width"))
    plist_->sublist("verbose object").set("debug cell header width", 34);

  if (!plist_->isParameter("primary variable key suffix")) {
    plist_->set<std::string>("primary variable key suffix", "total_component_concentration");
  }
  PK_Physical_Default::parseParameterList();

  // other parameters
  saturation_tolerance_ = plist_->get<double>("saturation tolerance", 1e-14);
  tcc_tag_current_ = Tag(plist_->get<std::string>("concentration tag current", tag_current_.get()));
  tcc_tag_next_ = Tag(plist_->get<std::string>("concentration tag next", tag_next_.get()));
}


void
Chemistry_PK::Initialize()
{
  // note, this is done here to allow the default value to be the time in state
  initial_conditions_time_ = plist_->get<double>("initial conditions time", S_->get_time());

  PK_Physical_Default::Initialize();
}


/* *******************************************************************
 * Take a timestep from t_old to t_new
 ******************************************************************* */
bool
Chemistry_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << t_old << " t1 = " << t_new << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;
  db_->WriteVector("C_old",
                   S_->GetPtr<CompositeVector>(key_, tcc_tag_current_).ptr(),
                   false,
                   S_->GetRecordSet(key_).subfieldnames());

  // these assertions fail with Amanzi, does Amanzi not call State::set_time?  --ETC
  // AMANZI_ASSERT(std::abs(S_->get_time(tag_current_) - t_old) < 1.e-4);
  // AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t_new) < 1.e-4);
  State_to_Solution(Tags::NEXT, *solution_);

  // set up the substate for faster access
  S_->GetW<CompositeVector>(key_, tcc_tag_next_, passwd_) =
    S_->Get<CompositeVector>(key_, tcc_tag_current_);
  updateSubstate(tag_current_);

  // Get the number of owned (non-ghost) cells for the mesh.
  AmanziMesh::Entity_ID num_cells =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // Now loop through all the cells and advance the chemistry.
  int max_itrs(0), imax(-1);
  int convergence_failure = 0;
  for (AmanziMesh::Entity_ID cell = 0; cell != num_cells; ++cell) {
    int num_itrs = advanceSingleCell_(cell, dt);
    if (num_itrs >= 0) {
      if (max_itrs < num_itrs) {
        max_itrs = num_itrs;
        imax = cell;
      }
    } else {
      // Convergence failure. Compute the next timestep size.
      convergence_failure = 1;
      break;
    }
  }

  // broadcast and check for global failed step
  checkForError_(convergence_failure, max_itrs, imax);

  // Compute the next global timestep.
  dt_next_ = timestep_controller_->getTimestep(dt, max_itrs, !convergence_failure);
  db_->WriteVector("C_new)",
                   S_->GetPtr<CompositeVector>(key_, tag_next_).ptr(),
                   false,
                   S_->GetRecordSet(key_).subfieldnames());

  return convergence_failure;
}


/* ******************************************************************
* Advance fields after a successful step
******************************************************************* */
void
Chemistry_PK::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  copyFields_(tag_current_, tag_next);

  if (tcc_tag_next_ != tag_next) {
    assign(key_, tag_next, tcc_tag_next_, *S_);
  }

  PK_Physical_Default::CommitStep(t_old, t_new, tag_next);
}


/* *******************************************************************
* I/O or error messages
******************************************************************* */
void
Chemistry_PK::checkForError_(int& ierr, int& max_itrs, int& max_itrs_cell) const
{
  int ierr_l(ierr);
  mesh_->getComm()->MaxAll(&ierr_l, &ierr, 1);
  if (ierr) return;

  Reductions::MaxLoc itrs_l{
    (double)max_itrs, mesh_->getMap(AmanziMesh::Entity_kind::CELL, false).GID(max_itrs_cell)
  };

  Reductions::MaxLoc itrs_g = Reductions::reduceAllMaxLoc(*mesh_->getComm(), itrs_l);
  max_itrs = itrs_g.value;
  max_itrs_cell = itrs_g.gid;

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "max Newton iterations: " << (int)itrs_g.value << " in cell " << itrs_g.gid
               << std::endl;
  }
}


} // namespace AmanziChemistry
} // namespace Amanzi
