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
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Mesh.hh"

// Operators
#include "Upwind.hh"

namespace Amanzi {
namespace Operators {

class UpwindFlux : public Upwind {
 public:
  UpwindFlux(Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      Upwind(mesh) {};
  ~UpwindFlux() {};

  // main methods
  void Init(Teuchos::ParameterList& plist);

  void Compute(const CompositeVector& flux, const CompositeVector& solution,
               const Kokkos::View<int*>& bc_model, CompositeVector& field);

 private:
  using Upwind::mesh_;
  using Upwind::face_comp_;

 private:
  int method_, order_;
  double tolerance_;
};


/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
void UpwindFlux::Init(Teuchos::ParameterList& plist)
{
  method_ = Operators::OPERATOR_UPWIND_FLUX;
  tolerance_ = plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
  order_ = plist.get<int>("polynomial order", 1);
}


/* ******************************************************************
* Upwind field uses flux. The result is placed in field.
* Upwinded field must be calculated on all faces of the owned cells.
****************************************************************** */
void UpwindFlux::Compute(
    const CompositeVector& flux, const CompositeVector& solution,
    const Kokkos::View<int*>& bc_model, CompositeVector& field)
{
  AMANZI_ASSERT(field.HasComponent("cell"));
  AMANZI_ASSERT(field.HasComponent(face_comp_));

  field.ScatterMasterToGhosted("cell");
  flux.ScatterMasterToGhosted("face");

  const auto& flx_face = flux.ViewComponent("face", true);
  // const Multivector_type& sol_face = *solution.ViewComponent("face", true);

  const auto& fld_cell = field.ViewComponent("cell", true);
  const auto& fld_boundary = field.ViewComponent("dirichlet_faces", true);
  const auto& ext_face_map = mesh_->exterior_face_map(true);
  const auto& face_map = mesh_->face_map(true);
  auto upw_face = field.ViewComponent(face_comp_, true);


  double flxmin = flx_face(0,0), flxmax = flx_face(0,0), tol;
  // Compute minimum value 
  for(int i = 0 ; i < flx_face.size(); ++i){
    flxmin = std::min(flxmin,flx_face(i,0)); 
    flxmax = std::max(flxmax,flx_face(i,0)); 
  }
  //flx_face.MinValue(&flxmin);
  //flx_face.MaxValue(&flxmax);
  tol = tolerance_ * std::max(fabs(flxmin), fabs(flxmax));

  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  int c1, c2, dir;
  double kc1, kc2;
  for (int f = 0; f < nfaces_wghost; ++f) {
    AmanziMesh::Entity_ID_View cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
    int ncells = cells.size();

    c1 = cells[0];
    kc1 = fld_cell(c1,0);

    mesh_->face_normal(f, false, c1, &dir);
    bool flag = (flx_face(f,0) * dir <= -tol);  // upwind flag


    if (ncells == 2) { 
      c2 = cells[1];
      kc2 = fld_cell(c2,0);

      // We average field on almost vertical faces. 
      if (fabs(flx_face(f,0)) <= tol) { 
        double v1 = mesh_->cell_volume(c1);
        double v2 = mesh_->cell_volume(c2);

        double tmp = v2 / (v1 + v2);
        upw_face(f,0) = kc1 * tmp + kc2 * (1.0 - tmp); 
      } else {
        upw_face(f,0) = (flag) ? kc2 : kc1; 
      }

    // We upwind only on inflow dirichlet faces.
    } else {
      upw_face(f,0) = kc1;
      if (bc_model[f] == OPERATOR_BC_DIRICHLET && flag) {
        upw_face(f,0) = fld_boundary(ext_face_map->getLocalElement(face_map->getGlobalElement(f)),0);
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

