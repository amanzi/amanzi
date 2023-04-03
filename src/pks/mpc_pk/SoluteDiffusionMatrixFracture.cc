/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Molecular diffusion coefficeint between fracture and matrix.
*/


#include "Teuchos_ParameterList.hpp"

#include "UniqueLocalIndex.hh"

#include "SoluteDiffusionMatrixFracture.hh"

namespace Amanzi {

/* *******************************************************************
* Constructor
******************************************************************* */
SoluteDiffusionMatrixFracture::SoluteDiffusionMatrixFracture(Teuchos::ParameterList& plist)
  : EvaluatorSecondary(plist), mol_diff_(0.0)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  domain_ = Keys::getDomain(my_keys_[0].first);

  saturation_key_ = plist_.get<std::string>("saturation key");
  tortuosity_key_ = plist_.get<std::string>("tortuosity key");
  porosity_key_ = plist_.get<std::string>("porosity key");
  if (plist_.isParameter("molecular diffusion")) {
    mol_diff_ = plist_.get<double>("molecular diffusion");
  }

  dependencies_.insert(std::make_pair(saturation_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));
  // dependencies_.insert(std::make_pair(tortuosity_key_, Tags::DEFAULT));  // static field
}


/* *******************************************************************
* Actual work is done here
******************************************************************* */
void
SoluteDiffusionMatrixFracture::Update_(State& S)
{
  Key key = my_keys_[0].first;
  auto mesh = S.GetMesh(domain_);
  auto mesh_parent = mesh->getParentMesh();

  const auto& sat_c = *S.Get<CompositeVector>(saturation_key_, Tags::DEFAULT).ViewComponent("cell");
  const auto& poro_c = *S.Get<CompositeVector>(porosity_key_, Tags::DEFAULT).ViewComponent("cell");
  const auto& tau_c = *S.Get<CompositeVector>(tortuosity_key_, Tags::DEFAULT).ViewComponent("cell");

  auto& result_c = *S.GetW<CompositeVector>(key, Tags::DEFAULT, key).ViewComponent("cell");
  int ncells = result_c.MyLength();

  for (int c = 0; c < ncells; ++c) {
    const AmanziGeometry::Point& xf = mesh->getCellCentroid(c);
    int f = mesh->getEntityParent(AmanziMesh::CELL, c);

    auto cells = mesh_parent->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int ndofs = cells.size();

    for (int k = 0; k < ndofs; ++k) {
      int pos = Operators::UniqueIndexFaceToCells(*mesh_parent, f, cells[k]);
      int c1 = cells[pos];

      const AmanziGeometry::Point& xc = mesh_parent->getCellCentroid(c1);
      double dist = norm(xc - xf);

      result_c[k][c] = sat_c[0][c1] * tau_c[0][c1] * poro_c[0][c1] * mol_diff_ / dist;
    }
  }
}


/* *******************************************************************
* Actual work is done here
******************************************************************* */
void
SoluteDiffusionMatrixFracture::EnsureCompatibility(State& S)
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
