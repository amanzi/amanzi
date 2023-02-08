/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Process Kernels

  Miscalleneous collection of simple non-member functions.
*/

#include "EvaluatorPrimary.hh"

#include "PK_Utils.hh"

namespace Amanzi {

/* ******************************************************************
* Deep copy of state fields.
****************************************************************** */
void
StateArchive::Add(std::vector<std::string> fields,
                  std::vector<std::string> evals,
                  std::vector<std::string> primary,
                  const Tag& tag,
                  const std::string& requestor)
{
  tag_ = tag;
  primary_ = primary;

  for (const auto& name : fields) fields_.emplace(name, S_->Get<CompositeVector>(name, tag));

  for (const auto& name : evals) {
    S_->GetEvaluator(name).Update(*S_, requestor);
    evals_.emplace(name, S_->Get<CompositeVector>(name, tag));
  }
}


/* ******************************************************************
* Deep copy of state fields.
****************************************************************** */
void
StateArchive::Restore(const std::string& passwd)
{
  for (auto it = fields_.begin(); it != fields_.end(); ++it) {
    S_->GetW<CompositeVector>(it->first, passwd) = it->second;

    if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "reverted field \"" << it->first << "\"" << std::endl;
    }
  }

  for (auto it = evals_.begin(); it != evals_.end(); ++it) {
    S_->GetW<CompositeVector>(it->first, tag_, it->first) = it->second;

    if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "reverted evaluator \"" << it->first << "\"" << std::endl;
    }
  }

  for (auto it = primary_.begin(); it != primary_.end(); ++it) {
    Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(
      S_->GetEvaluatorPtr(*it, tag_))
      ->SetChanged();

    if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "set primary solution \"" << *it << "\" to changed status" << std::endl;
    }
  }
}


/* *******************************************************************
* Copy: Evaluator (BASE) -> Field (prev_BASE)
******************************************************************* */
void
StateArchive::Swap(const std::string& passwd)
{
  for (auto it = fields_.begin(); it != fields_.end(); ++it) {
    std::string prev(it->first), next(it->first);
    auto pos = next.find("prev_");
    if (pos != std::string::npos) {
      next.erase(pos, 5);
      S_->GetW<CompositeVector>(prev, tag_, passwd) = S_->Get<CompositeVector>(next);
    }
  }
}


/* ******************************************************************
* Return a copy
****************************************************************** */
const CompositeVector&
StateArchive::get(const std::string& name)
{
  {
    auto it = fields_.find(name);
    if (it != fields_.end()) return it->second;
  }

  {
    auto it = evals_.find(name);
    if (it != evals_.end()) return it->second;
  }

  AMANZI_ASSERT(false);
}


/* ******************************************************************
* Average permeability tensor in horizontal direction.
****************************************************************** */
void
PKUtils_CalculatePermeabilityFactorInWell(const Teuchos::Ptr<State>& S,
                                          Teuchos::RCP<Epetra_Vector>& Kxy)
{
  if (!S->HasRecord("permeability", Tags::DEFAULT)) return;

  const auto& cv = S->Get<CompositeVector>("permeability", Tags::DEFAULT);
  cv.ScatterMasterToGhosted("cell");
  const auto& perm = *cv.ViewComponent("cell", true);

  int ncells_wghost = S->GetMesh()->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  int dim = perm.NumVectors();

  Kxy = Teuchos::rcp(new Epetra_Vector(S->GetMesh()->getMap(AmanziMesh::Entity_kind::CELL,true)));

  for (int c = 0; c < ncells_wghost; c++) {
    (*Kxy)[c] = 0.0;
    int idim = std::max(1, dim - 1);
    for (int i = 0; i < idim; i++) (*Kxy)[c] += perm[i][c];
    (*Kxy)[c] /= idim;
  }
}


/* ******************************************************************
* Return coordinate of mesh entity (
****************************************************************** */
AmanziGeometry::Point
PKUtils_EntityCoordinates(int id, AmanziMesh::Entity_ID kind, const AmanziMesh::Mesh& mesh)
{
  if (kind == AmanziMesh::Entity_kind::FACE) {
    return mesh.getFaceCentroid(id);
  } else if (kind == AmanziMesh::Entity_kind::CELL) {
    return mesh.getCellCentroid(id);
  } else if (kind == AmanziMesh::Entity_kind::NODE) {
    int d = mesh.getSpaceDimension();
    AmanziGeometry::Point xn(d);
    xn = mesh.getNodeCoordinate(id);
    return xn;
  } else if (kind == AmanziMesh::Entity_kind::EDGE) {
    return mesh.getEdgeCentroid(id);
  }
  return AmanziGeometry::Point();
}

} // namespace Amanzi
