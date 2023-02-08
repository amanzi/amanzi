/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Evaluator for Darcy velocity.
*/

#include "MFD3D_Diffusion.hh"
#include "UniqueLocalIndex.hh"

#include "DarcyVelocityEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* A simple constructor: create dependencies.
****************************************************************** */
DarcyVelocityEvaluator::DarcyVelocityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  AMANZI_ASSERT(my_keys_.size() > 0);
  vol_flowrate_key_ = plist.get<std::string>("volumetric flow rate key");
  dependencies_.insert(std::make_pair(vol_flowrate_key_, Tags::DEFAULT));
}


/* ******************************************************************
* A copy constructor.
****************************************************************** */
DarcyVelocityEvaluator::DarcyVelocityEvaluator(const DarcyVelocityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    vol_flowrate_key_(other.vol_flowrate_key_)
{}


/* ******************************************************************
* Clone with unclear yet purpose.
****************************************************************** */
Teuchos::RCP<Evaluator>
DarcyVelocityEvaluator::Clone() const
{
  return Teuchos::rcp(new DarcyVelocityEvaluator(*this));
}


/* ******************************************************************
* Required member function: basic algorithm.
****************************************************************** */
void
DarcyVelocityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  Key domain = plist_.get<std::string>("domain name");
  S.Get<CompositeVector>(vol_flowrate_key_).ScatterMasterToGhosted("face");

  const auto& flux = *S.Get<CompositeVector>(vol_flowrate_key_).ViewComponent("face", true);
  auto& result_c = *results[0]->ViewComponent("cell");

  const auto& fmap = *S.Get<CompositeVector>(vol_flowrate_key_).Map().Map("face", true);

  auto mesh = S.GetMesh(domain);
  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int dim = mesh->getSpaceDimension();

  WhetStone::MFD3D_Diffusion mfd(mesh);

  WhetStone::Polynomial gradient(dim, 1);
  AmanziMesh::Entity_ID_View cells;

  for (int c = 0; c < ncells_owned; c++) {
    const auto& faces = mesh->getCellFaces(c);
    int nfaces = faces.size();
    std::vector<WhetStone::Polynomial> solution(nfaces);

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      int g = fmap.FirstPointInElement(f);

      // the case of two DOFs on the face:
      if (fmap.ElementSize(f) == 2) {
        cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        if (cells.size() == 2) { g += Operators::UniqueIndexFaceToCells(*mesh, f, c); }
      }
      solution[n].Reshape(dim, 0);
      solution[n](0) = flux[0][g];
    }

    mfd.L2Cell(c, solution, solution, NULL, gradient);
    for (int i = 0; i < dim; i++) result_c[i][c] = -gradient(i + 1);
  }
}

} // namespace Flow
} // namespace Amanzi
