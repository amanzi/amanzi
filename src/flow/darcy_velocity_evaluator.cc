/*
  Evaluator for darcy_velocity(darcy_flux)

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "mfd3d_diffusion.hh"
#include "darcy_velocity_evaluator.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* A simple constructor: create dependencies.
****************************************************************** */
DarcyVelocityEvaluator::DarcyVelocityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  // hard-coded keys
  my_key_ = std::string("darcy_velocity");
  darcy_flux_key_ = std::string("darcy_flux");
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
  const Epetra_MultiVector& flux = *S->GetFieldData(darcy_flux_key_)->ViewComponent("face");
  Epetra_MultiVector& result_c = *(result->ViewComponent("cell", false));

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh();
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int dim = mesh->space_dimension();

  WhetStone::MFD3D_Diffusion mfd(mesh);

  AmanziGeometry::Point gradient(dim);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    std::vector<double> solution(nfaces);

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      solution[n] = flux[0][f] * dirs[n];
    }
  
    mfd.RecoverGradient_MassMatrix(c, solution, gradient);
    for (int i = 0; i < dim; i++) result_c[i][c] = gradient[i];
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi
