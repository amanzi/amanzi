/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Upwind a cell-centered field defined on a network of manifolds 
  using a given face-based flux.
*/

// TPLs
#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"

// Operators
#include "UniqueLocalIndex.hh"
#include "UpwindFluxManifolds.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Init method is not used.
****************************************************************** */
void
UpwindFluxManifolds::Init(Teuchos::ParameterList& plist)
{
  method_ = Operators::OPERATOR_UPWIND_FLUX_MANIFOLDS;
  tolerance_ = plist.get<double>("tolerance", OPERATOR_UPWIND_RELATIVE_TOLERANCE);
}


/* ******************************************************************
* Upwind cells -> faces.
****************************************************************** */
void
UpwindFluxManifolds::Compute(const CompositeVector& flux,
                             const std::vector<int>& bc_model,
                             CompositeVector& field)
{
  AMANZI_ASSERT(field.HasComponent("cell"));
  AMANZI_ASSERT(field.HasComponent("face"));

  flux.ScatterMasterToGhosted("face");
  field.ScatterMasterToGhosted("cell");

  const auto& flux_f = *flux.ViewComponent("face", true);
  const auto& field_c = *field.ViewComponent("cell", true);
  const auto& field_bf = *field.ViewComponent("boundary_face", true);
  auto& field_f = *field.ViewComponent("face", true);

  double flxmin, flxmax, tol;
  flux_f.MinValue(&flxmin);
  flux_f.MaxValue(&flxmax);
  tol = tolerance_ * std::max(fabs(flxmin), fabs(flxmax));

  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  AmanziMesh::Entity_ID_List cells;

  // multiple DOFs on faces require usage of block map
  const auto& fmap = *flux.ComponentMap("face", true);

  for (int f = 0; f < nfaces_wghost; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    int g = fmap.FirstPointInElement(f);

    // harmonic average over upwind cells
    if (ncells > 1) {
      double mean(0.0);
      for (int i = 0; i < ncells; ++i) {
        int c = cells[i];
        mean += field_c[0][c];
      }
      mean /= ncells;

      for (int i = 0; i < ncells; ++i) {
        field_f[0][g + i] = mean;
      }

      // upwind only on inflow Dirichlet faces
    } else {
      field_f[0][g] = field_c[0][cells[0]];
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        int bf = getFaceOnBoundaryBoundaryFace(*mesh_, f);
        field_f[0][g] = field_bf[0][bf];
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi
