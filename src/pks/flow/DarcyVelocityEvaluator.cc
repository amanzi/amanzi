/*
  Flow PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Evaluator for darcy_velocity(darcy_flux)
*/

#include "MFD3D_Diffusion.hh"
#include "UniqueLocalIndex.hh"

#include "DarcyVelocityEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* A simple constructor: create dependencies.
****************************************************************** */
DarcyVelocityEvaluator::DarcyVelocityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  // hard-coded keys
  my_key_ = plist.get<std::string>("darcy velocity key");
  darcy_flux_key_ = plist.get<std::string>("darcy flux key");
  dependencies_.insert(darcy_flux_key_);
}


/* ******************************************************************
* A copy constructor.
****************************************************************** */
DarcyVelocityEvaluator::DarcyVelocityEvaluator(const DarcyVelocityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    darcy_flux_key_(other.darcy_flux_key_) {}


/* ******************************************************************
* Clone with unclear yet purpose.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> DarcyVelocityEvaluator::Clone() const {
  return Teuchos::rcp(new DarcyVelocityEvaluator(*this));
}


/* ******************************************************************
* Required member function: basic algorithm.
****************************************************************** */
void DarcyVelocityEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result) 
{
  Key domain = plist_.get<std::string>("domain name");
  S->GetFieldData(darcy_flux_key_)->ScatterMasterToGhosted("face");

  const Epetra_MultiVector& flux = *S->GetFieldData(darcy_flux_key_)->ViewComponent("face", true);
  Epetra_MultiVector& result_c = *(result->ViewComponent("cell", false));

  const auto& fmap = *S->GetFieldData(darcy_flux_key_)->Map().Map("face", true);

  auto mesh = S->GetMesh(domain);
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int dim = mesh->space_dimension();

  WhetStone::MFD3D_Diffusion mfd(mesh);

  WhetStone::Polynomial gradient(dim, 1);
  AmanziMesh::Entity_ID_List cells;

  for (int c = 0; c < ncells_owned; c++) {
    const auto& faces = mesh->cell_get_faces(c);
    int nfaces = faces.size();
    std::vector<WhetStone::Polynomial> solution(nfaces);

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      int g = fmap.FirstPointInElement(f);

      // the case of two DOFs on the face:
      if (fmap.ElementSize(f) == 2) {
        mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        if (cells.size() == 2) {
          g += Operators::UniqueIndexFaceToCells(*mesh, f, c);
        }
      }
      solution[n].Reshape(dim, 0);
      solution[n](0) = flux[0][g];
    }
  
    mfd.L2Cell(c, solution, solution, NULL, gradient);
    for (int i = 0; i < dim; i++) result_c[i][c] = -gradient(i + 1);
  }
}

}  // namespace Flow
}  // namespace Amanzi
