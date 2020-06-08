/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A base class with default implementations of methods for a leaf of the PK
//! tree (a conservation equation, or similar).

/*!

``PK_MixinLeaf`` is a mixin class providing some functionality for PKs which
are defined on a single mesh, and represent a single process model.  Typically
all leaves of the PK tree will inherit from ``PK_MixinLeaf``.

* `"domain`" ``[string]`` **""**, e.g. `"surface`".

  Domains and meshes are 1-to-1, and the empty string refers to the main domain
or mesh.  PKs defined on other domains must specify which domain/mesh they refer
to.

* `"primary variable`" ``[string]``

  The primary variable associated with this PK, i.e. `"pressure`",
`"temperature`", `"surface_pressure`", etc.

*/

#ifndef AMANZI_PK_MIXIN_LEAF_HH_
#define AMANZI_PK_MIXIN_LEAF_HH_

#include "Teuchos_ParameterList.hpp"

#include "Debugger.hh"
#include "EvaluatorPrimary.hh"
#include "Key.hh"
#include "PK.hh"
#include "CompositeVector.hh"

namespace Amanzi {

template <class Base_t, class Data_t = CompositeVector,
          class DataFactory_t = CompositeVectorSpace>
class PK_MixinLeaf : public Base_t {
 public:
  PK_MixinLeaf(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
               const Teuchos::RCP<State>& S);

  // require data
  void Setup();
  void Initialize();


  void ChangedSolutionPK(const Key& tag);

  void CommitStep(const Key& tag_old, const Key& tag_new);
  void FailStep(const Key& tag_old, const Key& tag_new);

  // Accessor for debugger, for use by coupling MPCs
  Teuchos::Ptr<Debugger> debugger() { return db_.ptr(); }

  Teuchos::RCP<TreeVectorSpace> SolutionSpace();
  void SolutionToState(const Key& tag, const Key& suffix);
  void StateToState(const Key& tag_from, const Key& tag_to);

  // Construct all sub-PKs.  Leaf PKs have no children, hence have nothing to
  // do
  void ConstructChildren(){};


 protected:
  // name of domain, associated mesh
  Key domain_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  // solution and evaluator
  Key key_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> db_;

  using Base_t::S_;
  using Base_t::plist_;
};

template <class Base_t, class Data_t, class DataFactory_t>
PK_MixinLeaf<Base_t, Data_t, DataFactory_t>::PK_MixinLeaf(
  const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
  const Teuchos::RCP<State>& S)
  : Base_t(pk_tree, global_plist, S)
{
  // get the domain and mesh
  domain_ = plist_->template get<std::string>("domain name");
  mesh_ = S->GetMesh(domain_);

  // get the primary variable key
  key_ = Keys::readKey(*plist_, domain_, "primary variable");

  // create a debugger
  auto db_plist = Teuchos::sublist(plist_, "debugger");
  if (!db_plist->isSublist("verbose object")) {
    db_plist->set("verbose object", plist_->sublist("verbose object"));
  }
  db_ = Teuchos::rcp(new Debugger(mesh_, this->name(), *db_plist));
};

template <class Base_t, class Data_t, class DataFactory_t>
void
PK_MixinLeaf<Base_t, Data_t, DataFactory_t>::Setup()
{
  // note that this will inherit all of its metadata from requirements made
  // elsewhere.  Don't claim ownership here yet either, as requirements may
  // need to be set still.
  S_->template Require<Data_t, DataFactory_t>(key_, "");
}

template <class Base_t, class Data_t, class DataFactory_t>
void
PK_MixinLeaf<Base_t, Data_t, DataFactory_t>::Initialize()
{
  S_->GetRecordSet(key_).Initialize(plist_->sublist("initial conditions"));
}


// Mark, as changed, any primary variable evaluator owned by this PK
template <class Base_t, class Data_t, class DataFactory_t>
void
PK_MixinLeaf<Base_t, Data_t, DataFactory_t>::ChangedSolutionPK(const Key& tag)
{
  if (S_->HasEvaluator(key_, tag)) {
    auto eval = S_->GetEvaluatorPtr(key_, tag);
    auto eval_primary =
        Teuchos::rcp_dynamic_cast<EvaluatorPrimary<Data_t, DataFactory_t>>(eval);
    AMANZI_ASSERT(eval_primary.get());
    eval_primary->SetChanged();
  }
};

template <class Base_t, class Data_t, class DataFactory_t>
void
PK_MixinLeaf<Base_t, Data_t, DataFactory_t>::CommitStep(const Key& tag_old,
                                                        const Key& tag_new)
{
  Base_t::CommitStep(tag_old, tag_new);
  this->StateToState(tag_new, tag_old); // copy new into old
}

template <class Base_t, class Data_t, class DataFactory_t>
void
PK_MixinLeaf<Base_t, Data_t, DataFactory_t>::FailStep(const Key& tag_old,
                                                      const Key& tag_new)
{
  Base_t::FailStep(tag_old, tag_new);
  this->StateToState(tag_old, tag_new); // copy old into new

  // e.g. explicit PKs don't have an eval at the next tag.
  if (S_->HasEvaluator(key_, tag_new)) { this->ChangedSolutionPK(tag_new); }
}


template <class Base_t, class Data_t, class DataFactory_t>
Teuchos::RCP<TreeVectorSpace>
PK_MixinLeaf<Base_t, Data_t, DataFactory_t>::SolutionSpace()
{
  // This sort of implies that a TreeVectorSpace can contain any DataFactory_t.
  // If this ever gets used, it will have to be thought more about.
  return Teuchos::rcp(new TreeVectorSpace());
}


template <class Base_t, class Data_t, class DataFactory_t>
void
PK_MixinLeaf<Base_t, Data_t, DataFactory_t>::SolutionToState(const Key& tag,
                                                             const Key& suffix)
{
  Key key = key_ + suffix;
  S_->template Require<Data_t, DataFactory_t>(key, tag, key);
}

template <class Base_t, class Data_t, class DataFactory_t>
void
PK_MixinLeaf<Base_t, Data_t, DataFactory_t>::StateToState(const Key& tag_from,
                                                          const Key& tag_to)
{
  S_->template GetW<Data_t>(key_, tag_to, key_) =
    S_->template Get<Data_t>(key_, tag_from);
}


//
// Specialized class for CompositeVector with PrimaryVariable evaluator
//
// Note other classes would need their own evaluator types.
// -----------------------------------------------------------------------------
template <class Base_t>
class PK_MixinLeafCompositeVector
  : public PK_MixinLeaf<Base_t, CompositeVector, CompositeVectorSpace> {
 public:
  PK_MixinLeafCompositeVector(
      const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
      const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
      const Teuchos::RCP<State>& S)
      : PK_MixinLeaf<Base_t, CompositeVector, CompositeVectorSpace>(pk_tree, global_plist, S)
  {
    // create the primary variable evaluator NOTE: this is a CompositeVector --
    // this would need to get removed/fixed for a non-CV-based leaf... --etc
    S->FEList().sublist(key_).set("evaluator type", "primary variable");
  }
        
  Teuchos::RCP<TreeVectorSpace> SolutionSpace();
  void StateToSolution(TreeVector& soln, const Key& tag, const Key& suffix);
  void StateToState(const Key& tag_from, const Key& tag_to);
  

 protected:
  using PK_MixinLeaf<Base_t, CompositeVector, CompositeVectorSpace>::key_;
};

template <class Base_t>
Teuchos::RCP<TreeVectorSpace>
PK_MixinLeafCompositeVector<Base_t>::SolutionSpace()
{
  auto space = Teuchos::rcp(new TreeVectorSpace());
  space->SetData(
    this->S_
      ->template Require<CompositeVector, CompositeVectorSpace>(this->key_, "")
      .CreateSpace());
  return space;
}


template <class Base_t>
void
PK_MixinLeafCompositeVector<Base_t>::StateToSolution(TreeVector& soln,
                                                     const Key& tag,
                                                     const Key& suffix)
{
  Key key = this->key_ + suffix;
  if (soln.Data() == Teuchos::null)
    soln.SetData(this->S_->template GetPtrW<CompositeVector>(key, tag, key));

  // AMANZI_ASSERTS for now?
  AMANZI_ASSERT(soln.Data() ==
                this->S_->template GetPtr<CompositeVector>(key, tag));
}

template <class Base_t>
void
PK_MixinLeafCompositeVector<Base_t>::StateToState(const Key& tag_from,
        const Key& tag_to)
{
  PK_MixinLeaf<Base_t,CompositeVector,CompositeVectorSpace>::StateToState(tag_from,tag_to);
  this->ChangedSolutionPK(tag_to);
}



} // namespace Amanzi

#endif
