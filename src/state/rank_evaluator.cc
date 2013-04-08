/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an unchanging MPI Comm rank.

------------------------------------------------------------------------- */

#include "rank_evaluator.hh"
#include "Mesh.hh"

namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,RankModel> RankModel::fac_("mpi_comm_rank");

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
RankModel::RankModel(Teuchos::ParameterList& plist) :
    FieldEvaluator(plist),
    computed_once_(false) {

  my_key_ = plist.get<std::string>("evaluator name", "mpi_comm_rank");
  if (plist.isParameter("mesh name")) {
    my_mesh_ = plist.get<std::string>("mesh name");
  } else if (my_key_ == std::string("mpi_comm_rank")) {
    my_mesh_ = "domain";
  } else if (my_key_.length() > std::string("mpi_comm_rank").length()) {
    my_mesh_ = my_key_.substr(0,my_key_.length() - std::string("mpi_comm_rank").length()-1);
  } else {
    ASSERT(0);
  }
}


RankModel::RankModel(const RankModel& other) :
    FieldEvaluator(other),
    my_key_(other.my_key_),
    my_mesh_(other.my_mesh_),
    computed_once_(false) {}


Teuchos::RCP<FieldEvaluator>
RankModel::Clone() const {
  return Teuchos::rcp(new RankModel(*this));
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void RankModel::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Require the field and claim ownership.
  S->RequireField(my_key_, my_key_)->SetMesh(S->GetMesh(my_mesh_))
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(false);
  S->GetField(my_key_, my_key_)->set_io_vis(true);
}


// ---------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// ---------------------------------------------------------------------------
bool RankModel::HasFieldChanged(const Teuchos::Ptr<State>& S, Key request) {
  if (!computed_once_) {
    // field DOES have to be computed at least once, even if it never changes.
    UpdateField_(S);
    computed_once_ = true;
    return true;
  }

  return false;
}


// ---------------------------------------------------------------------------
// Answers the question, Has This Field's derivative with respect to Key
// wrt_key changed since it was last requested for Field Key reqest.
// Updates the derivative if needed.
// ---------------------------------------------------------------------------
inline
bool RankModel::HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
        Key request, Key wrt_key) {
  return false;  // no derivatives, though this should never be called
}


inline
bool RankModel::IsDependency(const Teuchos::Ptr<State>& S, Key key) const {
  return false;
}

inline
bool RankModel::ProvidesKey(Key key) const { return key == my_key_; }


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void RankModel::UpdateField_(const Teuchos::Ptr<State>& S) {
  // NOTE: RankModel owns its own data.
  Teuchos::RCP<CompositeVector> rank = S->GetFieldData(my_key_, my_key_);

  // initialize from mesh
  for (int c=0; c!=rank->size("cell"); ++c) {
    (*rank)("cell",0,c) = rank->mesh()->get_comm()->MyPID();
  }
}


} // namespace
