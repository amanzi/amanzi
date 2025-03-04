/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK
*/

// Amanzi
#include "PDE_HelperDiscretization.hh"

// Amanzi::Flow
#include "FlowDefs.hh"
#include "ApertureModel.hh"
#include "ApertureModelEvaluator.hh"
#include "ApertureModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
ApertureModelEvaluator::ApertureModelEvaluator(Teuchos::ParameterList& plist,
                                               Teuchos::RCP<ApertureModelPartition> apm)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), apm_(apm)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("aperture key"), Tags::DEFAULT));
  }

  // dependencies
  pressure_key_ = plist_.get<std::string>("pressure key");
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));

  use_stress_ = plist.get<bool>("use stress", false);
  if (use_stress_) dependencies_.insert(std::make_pair("hydrostatic_stress", Tags::DEFAULT));
}


ApertureModelEvaluator::ApertureModelEvaluator(const ApertureModelEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    apm_(other.apm_),
    pressure_key_(other.pressure_key_){};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
ApertureModelEvaluator::Clone() const
{
  return Teuchos::rcp(new ApertureModelEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
ApertureModelEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  UpdateOverburdenPressure_(S);

  auto& a_c = *results[0]->ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  int ncells = a_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    a_c[0][c] = apm_->second[(*apm_->first)[c]]->Aperture(p_c[0][c], pov_[c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
ApertureModelEvaluator::EvaluatePartialDerivative_(const State& S,
                                                   const Key& wrt_key,
                                                   const Tag& wrt_tag,
                                                   const std::vector<CompositeVector*>& results)
{
  UpdateOverburdenPressure_(S);

  auto& a_c = *results[0]->ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  int ncells = a_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    a_c[0][c] = apm_->second[(*apm_->first)[c]]->dAperturedPressure(p_c[0][c], pov_[c]);
  }
}


/* ******************************************************************
* Project hysdostatic pressure in matrix to fracture.
****************************************************************** */
void
ApertureModelEvaluator::UpdateOverburdenPressure_(const State& S)
{
  std::string key("hydrostatic_stress");
  bool flag = S.HasRecord(key);
  Teuchos::RCP<const Epetra_MultiVector> p_c;
  if (flag) {
    S.Get<CompositeVector>(key, Tags::DEFAULT).ScatterMasterToGhosted();
    p_c = S.Get<CompositeVector>(key, Tags::DEFAULT).ViewComponent("cell");
  }

  auto mesh_parent = S.GetMesh();
  auto mesh = S.GetMesh("fracture");
  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  pov_.resize(ncells);

  for (int c = 0; c < ncells; ++c) {
    if (flag) {
      int f = mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
      auto cells = mesh_parent->getFaceCells(f);
      pov_[c] = -((*p_c)[0][cells[0]] + (*p_c)[0][cells[1]]) / 2;
    } else {
      pov_[c] = apm_->second[(*apm_->first)[c]]->OverburdenPressure();
    }
  }
}

} // namespace Flow
} // namespace Amanzi
