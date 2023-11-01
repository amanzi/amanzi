/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Heat diffusion coefficient between fracture and matrix.
*/


#include "Teuchos_ParameterList.hpp"

#include "UniqueLocalIndex.hh"

#include "HeatDiffusionMatrixFracture.hh"

namespace Amanzi {

/* *******************************************************************
* Constructor
******************************************************************* */
HeatDiffusionMatrixFracture::HeatDiffusionMatrixFracture(Teuchos::ParameterList& plist)
  : EvaluatorSecondary(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  domain_ = Keys::getDomain(my_keys_[0].first);

  conductivity_key_ = plist_.get<std::string>("thermal conductivity key");
  aperture_key_ = plist_.get<std::string>("aperture key");
  dependencies_.insert(std::make_pair(conductivity_key_, Tags::DEFAULT));
}


/* *******************************************************************
* Actual work is done here
******************************************************************* */
void
HeatDiffusionMatrixFracture::Update_(State& S)
{
  Key key = my_keys_[0].first;
  auto mesh = S.GetMesh(domain_);
  auto mesh_parent = mesh->getParentMesh();

  const auto& conductivity_c =
    *S.Get<CompositeVector>(conductivity_key_, Tags::DEFAULT).ViewComponent("cell");

  const auto& aperture_c =
    *S.Get<CompositeVector>(aperture_key_, Tags::DEFAULT).ViewComponent("cell");

  auto& result_c = *S.GetW<CompositeVector>(key, Tags::DEFAULT, key).ViewComponent("cell");
  int ncells = result_c.MyLength();

  for (int c = 0; c < ncells; ++c) {
    const AmanziGeometry::Point& xf = mesh->getCellCentroid(c);
    int f = mesh->getEntityParent(AmanziMesh::CELL, c);

    auto cells = mesh_parent->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int ndofs = cells.size();

    double dist = aperture_c[0][c] / 2;
    for (int k = 0; k < ndofs; ++k) {
      int pos = Operators::UniqueIndexFaceToCells(*mesh_parent, f, cells[k]);
      int c1 = cells[pos];
      result_c[k][c] = conductivity_c[0][c1] / dist;
    }
  }
}


/* *******************************************************************
* Actual work is done here
******************************************************************* */
void
HeatDiffusionMatrixFracture::EnsureCompatibility(State& S)
{
  Key key = my_keys_[0].first;
  Tag tag = my_keys_[0].second;
  Key domain = Keys::getDomain(key);

  S.Require<CompositeVector, CompositeVectorSpace>(key, tag, key)
    .SetMesh(S.GetMesh(domain_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 2);

  // For dependencies, all we really care is whether there is an evaluator
  // or not.
  for (const auto& dep : dependencies_) { S.RequireEvaluator(dep.first, dep.second); }

  // It would be nice to verify mesh parenting
  for (const auto& dep : dependencies_) {
    auto& dep_fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);

    Key dep_domain = Keys::getDomain(dep.first);
    if (domain == dep_domain)
      dep_fac.SetMesh(S.GetMesh(domain));
    else
      dep_fac.SetMesh(S.GetMesh(domain)->getParentMesh());
  }
}

} // namespace Amanzi
