/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Upwind a cell-centered field (e.g. rel perm) using a given 
  face-based flux (e.g. Darcy flux).
*/

#ifndef AMANZI_UPWIND_FLUX_HH_
#define AMANZI_UPWIND_FLUX_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Mesh.hh"

// Operators
#include "Upwind.hh"

namespace Amanzi {
namespace Operators {

template<class Model>
class UpwindFlux : public Upwind<Model> {
 public:
  UpwindFlux(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
             Teuchos::RCP<const Model> model) :
      Upwind<Model>(mesh, model) {};
  ~UpwindFlux() {};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void Compute(const CompositeVector& flux, const CompositeVector& solution,
               const std::vector<int>& bc_model, CompositeVector& field);

 private:
  using Upwind<Model>::mesh_;
  using Upwind<Model>::model_;
  using Upwind<Model>::face_comp_;

 private:
  int method_, order_;
  double tolerance_;
};


/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
template<class Model>
void UpwindFlux<Model>::Init(Teuchos::ParameterList& plist)
{
  method_ = Operators::OPERATOR_UPWIND_FLUX;
  tolerance_ = plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
  order_ = plist.get<int>("polynomial order", 1);
}


/* ******************************************************************
* Upwind field uses flux. The result is placed in field.
* Upwinded field must be calculated on all faces of the owned cells.
****************************************************************** */
template<class Model>
void UpwindFlux<Model>::Compute(
    const CompositeVector& flux, const CompositeVector& solution,
    const std::vector<int>& bc_model, CompositeVector& field)
{
  AMANZI_ASSERT(field.HasComponent("cell"));
  AMANZI_ASSERT(field.HasComponent(face_comp_));

  field.ScatterMasterToGhosted("cell");
  flux.ScatterMasterToGhosted("face");

  const Epetra_MultiVector& flx_face = *flux.ViewComponent("face", true);
  // const Epetra_MultiVector& sol_face = *solution.ViewComponent("face", true);

  const Epetra_MultiVector& fld_cell = *field.ViewComponent("cell", true);
  const Epetra_MultiVector& fld_boundary = *field.ViewComponent("dirichlet_faces", true);
  const Epetra_Map& ext_face_map = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true);
  const Epetra_Map& face_map = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);
  Epetra_MultiVector& upw_face = *field.ViewComponent(face_comp_, true);

  double flxmin, flxmax, tol;
  flx_face.MinValue(&flxmin);
  flx_face.MaxValue(&flxmax);
  tol = tolerance_ * std::max(fabs(flxmin), fabs(flxmax));

  int nfaces_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL);
  AmanziMesh::Entity_ID_List cells;

  int c1, c2, dir;
  double kc1, kc2;
  for (int f = 0; f < nfaces_wghost; ++f) {
    mesh_->getFaceCells(f, AmanziMesh::Parallel_type::ALL, cells);
    int ncells = cells.size();

    c1 = cells[0];
    kc1 = fld_cell[0][c1];

    mesh_->getFaceNormal(f,  c1, &dir);
    bool flag = (flx_face[0][f] * dir <= -tol);  // upwind flag

    if (ncells == 2) { 
      c2 = cells[1];
      kc2 = fld_cell[0][c2];

      // We average field on almost vertical faces. 
      if (fabs(flx_face[0][f]) <= tol) { 
        double v1 = mesh_->getCellVolume(c1);
        double v2 = mesh_->getCellVolume(c2);

        double tmp = v2 / (v1 + v2);
        upw_face[0][f] = kc1 * tmp + kc2 * (1.0 - tmp); 
      } else {
        upw_face[0][f] = (flag) ? kc2 : kc1; 
      }

    // We upwind only on inflow dirichlet faces.
    } else {
      upw_face[0][f] = kc1;
      if (bc_model[f] == OPERATOR_BC_DIRICHLET && flag) {
        upw_face[0][f] = fld_boundary[0][ext_face_map.LID(face_map.GID(f))];
      }
      // if (bc_model[f] == OPERATOR_BC_NEUMANN) {
      //   upw_face[0][f] = ((*model_).*Value)(c, sol_face[0][f]);
      // }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi

#endif

