/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_UPWIND_HH_
#define AMANZI_UPWIND_HH_

#include <string>
#include <vector>

#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "Mesh.hh"
#include "VerboseObject.hh"

#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

template<class Model>
class Upwind {
 public:
  Upwind(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
         Teuchos::RCP<const Model> model)
      : mesh_(mesh), model_(model) {};
  ~Upwind() {};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void Compute(const CompositeVector& flux,
               const std::vector<int>& bc_model, const std::vector<double>& bc_value,
               const CompositeVector& field, CompositeVector& field_upwind);

 protected:
  Teuchos::RCP<VerboseObject> vo_;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const Model> model_;

  int method_, order_;  // method for upwinding
};


/* ******************************************************************
* Public Init method.
****************************************************************** */
template<class Model>
void Upwind<Model>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Upwind", plist));

  std::string name = plist.get<std::string>("method", "upwind with flux");
  if (name == "upwind with gravity") {
    method_ = Operators::OPERATOR_UPWIND_WITH_CONSTANT_VECTOR;
  } else if (name == "upwind with flux") {
    method_ = Operators::OPERATOR_UPWIND_WITH_FLUX;
  } else if (name == "arithmetic mean") {
    method_ = Operators::OPERATOR_ARITHMETIC_MEAN;
  } 

  order_ = plist.get<int>("order", 1);
}


/* ******************************************************************
* Flux-based upwind.
****************************************************************** */
template<class Model>
void Upwind<Model>::Compute(
    const CompositeVector& flux,
    const std::vector<int>& bc_model, const std::vector<double>& bc_value,
    const CompositeVector& field, CompositeVector& field_upwind)
{
  ASSERT(field.HasComponent("cell"));
  ASSERT(field_upwind.HasComponent("face"));

  Teuchos::OSTab tab = vo_->getOSTab();

  field.ScatterMasterToGhosted("cell");
  flux.ScatterMasterToGhosted("face");

  const Epetra_MultiVector& u = *flux.ViewComponent("face", true);
  const Epetra_MultiVector& fcells = *field.ViewComponent("cell", true);
  const Epetra_MultiVector& ffaces = *field.ViewComponent("face", true);

  Epetra_MultiVector& upw = *field_upwind.ViewComponent("face", true);
  upw.PutScalar(0.0);

  double umin, umax;
  u.MinValue(&umin);
  u.MaxValue(&umax);
  double tol = OPERATOR_UPWIND_RELATIVE_TOLERANCE * std::max(fabs(umin), fabs(umax));

  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      bool flag = (u[0][f] * dirs[n] < -tol);  // upwind flag
      
      if (bc_model[f] == OPERATOR_BC_NONE && fabs(u[0][f]) <= tol) { 
        upw[0][f] += fcells[0][c] / 2;  // Almost vertical face.
      } else if (bc_model[f] == OPERATOR_BC_FACE_DIRICHLET && flag) {
        upw[0][f] = model_->Value(c, bc_value[f]);
      } else if (bc_model[f] == OPERATOR_BC_FACE_NEUMANN && flag) {
        // upw[0][f] = model_->Value(c, ffaces[0][f]);
        upw[0][f] = fcells[0][c];
      } else if (!flag) {
        upw[0][f] = fcells[0][c];
      }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi

#endif

