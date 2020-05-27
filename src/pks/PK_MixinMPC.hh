/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! The mixin with default implementations for a Multi-Process Coordinator (MPC)

/*!

Mixin for the default MPC class.  A multi process coordinator (MPC) is a PK
which coordinates other PKs.  Each of these coordinated PKs may be MPCs
themselves, or physical PKs.  Note this does NOT provide a full implementation
of PK -- it does not supply the AdvanceStep() method.  Therefore this class
cannot be instantiated, but must be inherited by derived classes which finish
supplying the functionality.  Instead, this provides the data structures and
methods (which may be overridden by derived classes) for managing multiple
PKs.

Most of these methods simply loop through the coordinated PKs, calling their
respective methods.

* `"PKs`" ``[parameter-list]`` A sublist including all PK specs

*/

#ifndef AMANZI_PK_MIXIN_MPC_HH_
#define AMANZI_PK_MIXIN_MPC_HH_

#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"

namespace Amanzi {

template <class Base_t, class PK_Contained_t>
class PK_MixinMPC : public Base_t {
 public:
  typedef std::vector<Teuchos::RCP<PK_Contained_t>> SubPKList;

  PK_MixinMPC(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
              const Teuchos::RCP<State>& S);

  // Uses a factory to construct child PKs.
  void ConstructChildren();

  // This should NOT be a part of the public interface, and is only here for
  // testing.  Should make this private and friend the test generation
  // function or something?  Children should only be created by
  // ConstructChildren().
  void SetChildren(const SubPKList& sub_pks) { sub_pks_ = sub_pks; }

  // PK methods
  // -- sets up sub-PKs
  void Setup();

  // -- calls all sub-PK initialize() methods
  void Initialize();

  // Returns validity of the step taken from tag_old to tag_new
  bool ValidStep(const Key& tag_old, const Key& tag_new);

  // Do work that can only be done if we know the step was successful.
  void CommitStep(const Key& tag_old, const Key& tag_new);

  // Revert a step from tag_new back to tag_old
  void FailStep(const Key& tag_old, const Key& tag_new);

  // Calculate any diagnostics at tag, currently used for visualization.
  void CalculateDiagnostics(const Key& tag);

  // Mark, as changed, any primary variable evaluator owned by this PK
  void ChangedSolutionPK(const Key& tag);

  // Generates the TreeVectorSpace that represents the solution vector.
  Teuchos::RCP<TreeVectorSpace> SolutionSpace();

  // Ensure consistency between a time integrator's view of data (TreeVector)
  // and the dag's view of data (dictionary of CompositeVectors).
  //
  // Pull data from the state into a TreeVector.  If the TreeVector has no
  // data, copy pointers.  Otherwise ensure pointers are equal.
  //
  // These almost certainly are implemented by default.
  void StateToSolution(TreeVector& soln, const Key& tag, const Key& suffix);

  // Ensure consistency between a time integrator's view of data (TreeVector)
  // and the dag's view of data (dictionary of CompositeVectors).
  //
  // Push data from a TreeVector into State.  If the State has no data,
  // require it (allowing time integrators to require intermediate steps,
  // etc).  If the TreeVector has data and State doesn't, copy pointers.  If
  // the State does have this data, ensure consistency.
  //
  // These almost certainly are implemented by default.
  void SolutionToState(const Key& tag, const Key& suffix);

  void StateToState(const Key& tag_from, const Key& tag_to);

 protected:
  // list of the PKs coupled by this MPC
  SubPKList sub_pks_;

  // plists needed to be saved for ConstructChildren
  ParameterList_ptr_type pk_tree_;
  ParameterList_ptr_type global_plist_;
};


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
PK_MixinMPC<Base_t,PK_Contained_t>::PK_MixinMPC(
              const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
              const Teuchos::RCP<State>& S)
  : pk_tree_(pk_tree),
    global_plist_(global_plist),
    Base_t(pk_tree, global_plist, S) {}
  
// -----------------------------------------------------------------------------
// Create child PKs.
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t,PK_Contained_t>::ConstructChildren()
{
  auto child_pks_list =
    this->plist_->template get<Teuchos::Array<std::string>>("PKs order");

  PKFactory fac;
  for (const std::string& child_pk_name : child_pks_list) {
    auto pk_as_pk = fac.CreatePK(child_pk_name, pk_tree_,
            global_plist_, this->S_);
    auto pk_as_contained_t = Teuchos::rcp_dynamic_cast<PK_Contained_t>(pk_as_pk);
    if (pk_as_contained_t == Teuchos::null) {
      Errors::Message msg;
      msg << "PK_MixinMPC attempted to cast PK \"" << child_pk_name
          << "\" to type \"" << typeid(PK_Contained_t).name() << "\" failed.";
      throw(msg);
    }
    sub_pks_.push_back(pk_as_contained_t);
  }

  for (auto& pk : sub_pks_) pk->ConstructChildren();
}


// -----------------------------------------------------------------------------
// Loop over sub-PKs, calling their Setup methods
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t, PK_Contained_t>::Setup()
{
  for (auto& pk : sub_pks_) pk->Setup();
  Base_t::Setup();
}

// -----------------------------------------------------------------------------
// Loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t, PK_Contained_t>::Initialize()
{
  for (auto& pk : sub_pks_) pk->Initialize();
  Base_t::Initialize();
}

// -----------------------------------------------------------------------------
// Valid implies ALL are valid
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
bool
PK_MixinMPC<Base_t, PK_Contained_t>::ValidStep(const Key& tag_old,
                                               const Key& tag_new)
{
  for (auto& pk : sub_pks_) {
    if (!pk->ValidStep(tag_old, tag_new)) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
// Call all commits
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t, PK_Contained_t>::CommitStep(const Key& tag_old,
                                                const Key& tag_new)
{
  for (auto& pk : sub_pks_) pk->CommitStep(tag_old, tag_new);
}

// -----------------------------------------------------------------------------
// Call all fails
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t, PK_Contained_t>::FailStep(const Key& tag_old,
                                              const Key& tag_new)
{
  for (auto& pk : sub_pks_) pk->FailStep(tag_old, tag_new);
}

// -----------------------------------------------------------------------------
// Calculate any diagnostics at tag
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t, PK_Contained_t>::CalculateDiagnostics(const Key& tag)
{
  for (auto& pk : sub_pks_) pk->CalculateDiagnostics(tag);
}

// -----------------------------------------------------------------------------
// Mark, as changed, any primary variable evaluator owned by this PK
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t, PK_Contained_t>::ChangedSolutionPK(const Key& tag)
{
  for (auto& pk : sub_pks_) pk->ChangedSolutionPK(tag);
}


// -----------------------------------------------------------------------------
// Generates the TreeVectorSpace that represents the solution vector.
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
Teuchos::RCP<TreeVectorSpace>
PK_MixinMPC<Base_t, PK_Contained_t>::SolutionSpace()
{
  auto space = Teuchos::rcp(new TreeVectorSpace());
  for (auto& pk : sub_pks_) space->PushBack(pk->SolutionSpace());
  return space;
}


// -----------------------------------------------------------------------------
// Ensure consistency between a time integrator's view of data (TreeVector)
// and the dag's view of data (dictionary of CompositeVectors).
//
// Pull data from the state into a TreeVector.  If the TreeVector has no
// data, copy pointers.  Otherwise ensure pointers are equal.
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t, PK_Contained_t>::StateToSolution(TreeVector& soln,
                                                     const Key& tag,
                                                     const Key& suffix)
{
  AMANZI_ASSERT(soln.size() == sub_pks_.size());
  // if (soln.size() != sub_pks_.size()) {
  //   AMANZI_ASSERT(soln.size() == 0);
  //   for (auto &pk : sub_pks_) {
  //     auto tv = Teuchos::rcp(new TreeVector());
  //     soln.PushBack(tv);
  //   }
  // }

  int i = 0;
  for (auto& pk : sub_pks_) {
    pk->StateToSolution(*soln.SubVector(i), tag, suffix);
    ++i;
  }
}

// -----------------------------------------------------------------------------
// Ensure consistency between a time integrator's view of data (TreeVector)
// and the dag's view of data (dictionary of CompositeVectors).
//
// Push data from a TreeVector into State.  If the State has no data,
// require it (allowing time integrators to require intermediate steps,
// etc).  If the TreeVector has data and State doesn't, copy pointers.  If
// the State does have this data, ensure consistency.
// -----------------------------------------------------------------------------
template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t, PK_Contained_t>::SolutionToState(const Key& tag,
                                                     const Key& suffix)
{
  for (auto& pk : sub_pks_) pk->SolutionToState(tag, suffix);
}

template <class Base_t, class PK_Contained_t>
void
PK_MixinMPC<Base_t, PK_Contained_t>::StateToState(const Key& tag_from,
                                                  const Key& tag_to)
{
  for (auto& pk : sub_pks_) pk->StateToState(tag_from, tag_to);
}

} // namespace Amanzi

#endif
