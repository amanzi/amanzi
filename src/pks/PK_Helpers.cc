/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A set of helper functions for doing common things in PKs.
#include "MeshAlgorithms.hh"
#include "PK_Helpers.hh"

namespace Amanzi {

bool
aliasVector(State& S, const Key& key, const Tag& target, const Tag& alias)
{
  if (S.HasEvaluator(key, target) && !S.HasEvaluator(key, alias)) {
    S.RequireEvaluator(key, alias, true);
    S.GetRecordSetW(key).AliasRecord(target, alias);
    return true;
  }
  return false;
}


// -----------------------------------------------------------------------------
// Given a vector, apply the Dirichlet data to that vector.
// -----------------------------------------------------------------------------
void
applyDirichletBCs(const Operators::BCs& bcs, CompositeVector& u)
{
  if (u.HasComponent("face")) {
    Epetra_MultiVector& u_f = *u.ViewComponent("face", false);
    for (unsigned int f = 0; f != u_f.MyLength(); ++f) {
      if (bcs.bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) { u_f[0][f] = bcs.bc_value()[f]; }
    }

  } else if (u.HasComponent("boundary_face")) {
    Epetra_MultiVector& u_bf = *u.ViewComponent("boundary_face", false);
    const Epetra_MultiVector& u_c = *u.ViewComponent("cell", false);
    const Epetra_Map& vandalay_map =
      u.Mesh()->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false);
    const Epetra_Map& face_map = u.Mesh()->getMap(AmanziMesh::Entity_kind::FACE, false);

    for (int bf = 0; bf != u_bf.MyLength(); ++bf) {
      AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));
      if (bcs.bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_bf[0][bf] = bcs.bc_value()[f];
      }
    }
  }
}


// -----------------------------------------------------------------------------
// Given a vector and a face ID, get the value at that location.
//
// Looks in the following order:
//  -- face component
//  -- boundary Dirichlet data
//  -- boundary_face value
//  -- internal cell
// -----------------------------------------------------------------------------
double
getFaceOnBoundaryValue(AmanziMesh::Entity_ID f, const CompositeVector& u, const Operators::BCs& bcs)
{
  if (u.HasComponent("face")) {
    return (*u.ViewComponent("face", false))[0][f];
  } else if (bcs.bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
    return bcs.bc_value()[f];
    // } else if (u.HasComponent("boundary_face")) {
    //   AmanziMesh::Entity_ID bf = getFaceOnBoundaryBoundaryFace(*u.Mesh(), f);
    //   return (*u.ViewComponent("boundary_face",false))[0][bf];
  } else {
    auto c = getFaceOnBoundaryInternalCell(*u.Mesh(), f);
    return (*u.ViewComponent("cell", false))[0][c];
  }
  return -1;
}


// -----------------------------------------------------------------------------
// Get the directional int for a face that is on the boundary.
// -----------------------------------------------------------------------------
int
getBoundaryDirection(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f)
{
  auto cells = mesh.getFaceCells(f);
  AMANZI_ASSERT(cells.size() == 1);
  const auto& [faces, dirs] = mesh.getCellFacesAndDirections(cells[0]);
  return dirs[std::find(faces.begin(), faces.end(), f) - faces.begin()];
}


// -----------------------------------------------------------------------------
// Marks a primary evaluator as changed.
// -----------------------------------------------------------------------------
bool
changedEvaluatorPrimary(const Key& key, const Tag& tag, State& S, bool or_die)
{
  bool changed = false;
  Teuchos::RCP<Evaluator> eval = S.GetEvaluatorPtr(key, tag);
  Teuchos::RCP<EvaluatorPrimaryCV> eval_pv = Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(eval);
  if (eval_pv == Teuchos::null) {
    if (or_die) {
      Errors::Message msg;
      msg << "Expected primary variable evaluator for " << key << " @ " << tag.get();
      Exceptions::amanzi_throw(msg);
    }
  } else {
    eval_pv->SetChanged();
    changed = true;
  }
  return changed;
}


// -----------------------------------------------------------------------------
// Get a primary variable evaluator for a key at tag
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireEvaluatorPrimary(const Key& key, const Tag& tag, State& S, const Key& owner, bool or_die)
{
  CompositeVectorSpace& cvs = S.Require<CompositeVector, CompositeVectorSpace>(key, tag, owner);

  // first check, is there one already
  if (S.HasEvaluator(key, tag)) {
    // if so, make sure it is primary
    Teuchos::RCP<Evaluator> eval = S.GetEvaluatorPtr(key, tag);
    Teuchos::RCP<EvaluatorPrimaryCV> eval_pv = Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(eval);
    if (or_die && eval_pv == Teuchos::null) {
      Errors::Message msg;
      msg << "Expected primary variable evaluator for " << key << " @ " << tag.get();
      Exceptions::amanzi_throw(msg);
    }

  } else {
    // if not, create one, only at this tag, not to be shared across tags.  By
    // this, we mean we don't stick the "type" = "primary" back into the
    // evaluator list -- this allows "copy evaluators" e.g. "water content at the
    // old tag" to differ from the standard evalulator, e.g. "water content at
    // the new tag" which is likely a secondary variable evaluator.
    Teuchos::ParameterList plist(key);
    plist.set("evaluator type", "primary variable");
    plist.set("tag", tag.get());
    auto eval_pv = Teuchos::rcp(new EvaluatorPrimaryCV(plist));
    S.SetEvaluator(key, tag, eval_pv);
  }
  return cvs;
}


// -----------------------------------------------------------------------------
// Require a vector and a primary variable evaluator at current tag(s).
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireEvaluatorAtCurrent(const Key& key, const Tag& tag, State& S, const Key& owner)
{
  CompositeVectorSpace& cvs = S.Require<CompositeVector, CompositeVectorSpace>(key, tag);
  if (!owner.empty()) {
    Key tag_current_owner = S.GetRecord(key, tag).owner();
    if (tag_current_owner.empty()) {
      S.Require<CompositeVector, CompositeVectorSpace>(key, tag, owner);
      requireEvaluatorPrimary(key, tag, S, owner);
    }

    if (tag != Tags::CURRENT) {
      // NOTE, this differs from requireAtNext() -- the CURRENT copy is required
      // to recover from a failed step.
      S.Require<CompositeVector, CompositeVectorSpace>(key, Tags::CURRENT);
      Key current_owner = S.GetRecord(key, Tags::CURRENT).owner();
      if (current_owner.empty()) {
        requireEvaluatorAssign(key, Tags::CURRENT, S, owner);
      }
    }
  } else {
    S.RequireEvaluator(key, tag);
  }
  return cvs;
}


// -----------------------------------------------------------------------------
// Require a vector and a primary variable evaluator at next tag(s).
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireEvaluatorAtNext(const Key& key, const Tag& tag, State& S, bool managed_here, const Key& owner)
{
  CompositeVectorSpace& cvs = S.Require<CompositeVector, CompositeVectorSpace>(key, tag);
  if (!owner.empty()) {
    managed_here = true;
    requireEvaluatorPrimary(key, tag, S, owner);
  } else {
    S.RequireEvaluator(key, tag);
  }

  if (managed_here && tag != Tags::NEXT) {
    aliasVector(S, key, tag, Tags::NEXT);
  }
  return cvs;
}


// -----------------------------------------------------------------------------
// Require assignment evaluator, which allows tracking old data.
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireEvaluatorAssign(const Key& key, const Tag& tag, State& S, const Key& owner)
{
  // in the future, this will likely derive from primary instead of just being
  // primary.  This will allow confirming that the times are the same.
  return requireEvaluatorPrimary(key, tag, S, owner, false);
}


// -----------------------------------------------------------------------------
// Assign if it is an assignment evaluator.
// -----------------------------------------------------------------------------
void
assign(const Key& key, const Tag& tag_dest, const Tag& tag_source, State& S)
{
  S.GetEvaluator(key, tag_source).Update(S, Keys::getKey(key, tag_dest));
  if (S.HasEvaluator(key, tag_dest)) {
    bool changed = changedEvaluatorPrimary(key, tag_dest, S, false);
    if (changed) S.Assign(key, tag_dest, tag_source);
  }
}


void
copyMeshCoordinatesToVector(const AmanziMesh::Mesh& mesh, CompositeVector& vec)
{
  Epetra_MultiVector& nodes = *vec.ViewComponent("node", true);

  int ndim = mesh.getSpaceDimension();
  AmanziGeometry::Point nc;
  for (int i = 0; i != nodes.MyLength(); ++i) {
    nc = mesh.getNodeCoordinate(i);
    for (int j = 0; j != ndim; ++j) nodes[j][i] = nc[j];
  }
}

void
copyVectorToMeshCoordinates(const CompositeVector& vec, AmanziMesh::Mesh& mesh)
{
  const Epetra_MultiVector& nodes = *vec.ViewComponent("node", true);
  int ndim = mesh.getSpaceDimension();

  Amanzi::AmanziMesh::Entity_ID_View node_ids("node_ids", nodes.MyLength());
  Amanzi::AmanziMesh::Point_View new_positions("new_positions", nodes.MyLength());
  for (int n = 0; n != nodes.MyLength(); ++n) {
    node_ids[n] = n;
    if (mesh.getSpaceDimension() == 2) {
      new_positions[n] = Amanzi::AmanziGeometry::Point{ nodes[0][n], nodes[1][n] };
    } else {
      new_positions[n] = Amanzi::AmanziGeometry::Point{ nodes[0][n], nodes[1][n], nodes[2][n] };
    }
  }
  Amanzi::AmanziMesh::deform(mesh, node_ids, new_positions);
}

int
commMaxValLoc(const Comm_type& comm, const ValLoc& local, ValLoc& global)
{
  MpiComm_type const* mpi_comm = dynamic_cast<const MpiComm_type*>(&comm);
  const MPI_Comm& mpi_comm_raw = mpi_comm->Comm();
  return MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpi_comm_raw);
}

ValLoc
maxValLoc(const Epetra_Vector& vec)
{
  ValLoc local{ 0., 0 };
  for (int i = 0; i != vec.MyLength(); ++i) {
    if (vec[i] > local.value) {
      local.value = vec[i];
      local.gid = vec.Map().GID(i);
    }
  }
  ValLoc global{ 0., 0 };
  int ierr = commMaxValLoc(vec.Comm(), local, global);
  AMANZI_ASSERT(!ierr);
  return global;
}

} // namespace Amanzi
