/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_UPWIND_MFD_HH_
#define AMANZI_UPWIND_MFD_HH_

#include <string>
#include <vector>

#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "Mesh.hh"
#include "VerboseObject.hh"

#include "Upwind.hh"

namespace Amanzi {
namespace Operators {

template<class Model>
class UpwindMFD : public Upwind<Model> {
 public:
  UpwindMFD(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                 Teuchos::RCP<const Model> model)
      : Upwind<Model>(mesh, model) {};
  ~UpwindMFD() {};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void Compute(const CompositeVector& flux,
               const std::vector<int>& bc_model, const std::vector<double>& bc_value,
               const CompositeVector& field, CompositeVector& field_upwind,
               double (Model::*Value)(int, double) const);

 private:
  using Upwind<Model>::vo_;
  using Upwind<Model>::mesh_;
  using Upwind<Model>::model_;

 private:
  int method_, order_;
};


/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
template<class Model>
void UpwindMFD<Model>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("UpwindMFD", plist));

  method_ = Operators::OPERATOR_UPWIND_FLUX;

  order_ = plist.get<int>("order", 1);
}


/* ******************************************************************
* Flux-based upwind.
****************************************************************** */
template<class Model>
void UpwindMFD<Model>::Compute(
    const CompositeVector& flux,
    const std::vector<int>& bc_model, const std::vector<double>& bc_value,
    const CompositeVector& field, CompositeVector& field_upwind,
    double (Model::*Value)(int, double) const)
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
      } else if (bc_model[f] == OPERATOR_BC_DIRICHLET && flag) {
        upw[0][f] = ((*model_).*Value)(c, bc_value[f]);
      } else if (bc_model[f] == OPERATOR_BC_NEUMANN && flag) {
        // upw[0][f] = ((*model_).*Value)(c, ffaces[0][f]);
        upw[0][f] = fcells[0][c];
      } else if (bc_model[f] == OPERATOR_BC_MIXED && flag) {
        // upw[0][f] = ((*model_).*Value)(c, ffaces[0][f]);
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

