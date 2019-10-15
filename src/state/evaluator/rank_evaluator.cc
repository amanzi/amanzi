/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "rank_evaluator.hh"
#include "Mesh.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
RankEvaluator::RankEvaluator(Teuchos::ParameterList& plist)
  : Evaluator(plist), computed_once_(false)
{
  my_key_ = plist.get<std::string>("evaluator name", "mpi_comm_rank");
  if (plist.isParameter("mesh name")) {
    my_mesh_ = plist.get<std::string>("mesh name");
  } else if (my_key_ == std::string("mpi_comm_rank")) {
    my_mesh_ = "domain";
  } else if (my_key_.length() > std::string("mpi_comm_rank").length()) {
    my_mesh_ = my_key_.substr(
      0, my_key_.length() - std::string("mpi_comm_rank").length() - 1);
  } else {
    AMANZI_ASSERT(0);
  }
}

RankEvaluator::RankEvaluator(const RankEvaluator& other)
  : Evaluator(other),
    my_key_(other.my_key_),
    my_mesh_(other.my_mesh_),
    computed_once_(false)
{}

Teuchos::RCP<Evaluator>
RankEvaluator::Clone() const
{
  return Teuchos::rcp(new RankEvaluator(*this));
}

// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void
RankEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // Require the field and claim ownership.
  S->RequireField(my_key_, my_key_)
    .SetMesh(S->GetMesh(my_mesh_))
    .SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(false);
  S->GetField(my_key_, my_key_)->set_io_vis(true);
}

// ---------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// ---------------------------------------------------------------------------
bool
RankEvaluator::Update(const Teuchos::Ptr<State>& S, Key request)
{
  if (!computed_once_) {
    // field DOES have to be computed at least once, even if it never changes.
    Update_(S);
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
inline bool
RankEvaluator::UpdateDerivative(const Teuchos::Ptr<State>& S, Key request,
                                Key wrt_key)
{
  return false; // no derivatives, though this should never be called
}

inline bool
RankEvaluator::IsDependency(const Teuchos::Ptr<State>& S, Key key) const
{
  return false;
}

inline bool
RankEvaluator::ProvidesKey(Key key) const
{
  return key == my_key_;
}

// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
RankEvaluator::Update_(const Teuchos::Ptr<State>& S)
{
  // NOTE: RankEvaluator owns its own data.
  Teuchos::RCP<CompositeVector> rank = S->GetFieldData(my_key_, my_key_);
  MultiVector_type& rank_c = *rank->ViewComponent("cell", false);

  // initialize from mesh
  for (int c = 0; c != rank_c.getLocalLength(); ++c) {
    rank_c[0][c] = rank->Comm().MyPID();
  }
}

std::string
RankEvaluator::WriteToString() const
{
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: independent" << std::endl
         << std::endl;
  return result.str();
}

} // namespace Amanzi
