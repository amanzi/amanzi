/*
This is the transport component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "mfd3d.hpp"
#include "tensor.hpp"

#include "Transport_PK.hh"


namespace Amanzi {
namespace AmanziTransport {

/* *******************************************************************
 * Calculate a dispersive tensor the from Darcy fluxes. The flux is
 * assumed to be scaled by face area.
 ****************************************************************** */
void Transport_PK::CalculateDispersionTensor()
{
  const Epetra_Vector& darcy_flux = TS_nextBIG->ref_darcy_flux();
  const Epetra_Vector& ws  = TS_nextBIG->ref_water_saturation();
  const Epetra_Vector& phi = TS_nextBIG->ref_porosity();

  AmanziMesh::Entity_ID_List nodes, faces;
  AmanziGeometry::Point velocity(dim), flux(dim);
  WhetStone::Tensor T(dim, 2);

  for (int c = 0; c < ncells_owned; c++) {
    if (dispersivity_model == TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC) {
      for (int i = 0; i < dim; i++) dispersion_tensor[c](i, i) = dispersivity_longitudinal;
    } else {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      int num_good_corners = 0;
      for (int n = 0; n < nnodes; n++) {
        int v = nodes[n];
        mesh_->node_get_cell_faces(v, c, AmanziMesh::USED, &faces);
        int nfaces = faces.size();

        for (int i = 0; i < dim; i++) {
          int f = faces[i];
          const AmanziGeometry::Point& normal = mesh_->face_normal(f);
          T.add_row(i, normal);
          flux[i] = darcy_flux[f];
        }

        T.inverse();
        velocity += T * flux;
        num_good_corners++;  // each corners is good temporary (lipnikov@lanl.gov)
      }
      velocity /= num_good_corners;

      double velocity_value = norm(velocity);
      double anisotropy = dispersivity_longitudinal - dispersivity_transverse;

      for (int i = 0; i < dim; i++) {
        dispersion_tensor[c](i, i) = dispersivity_transverse * velocity_value;
        for (int j = i; j < dim; j++) {
          double s = anisotropy * velocity[i] * velocity[j];
          if (velocity_value) s /= velocity_value;
          dispersion_tensor[c](j, i) = dispersion_tensor[c](i, j) += s;
        }
      }
    }

    double vol_phi_ws = mesh_->cell_volume(c) * phi[c] * ws[c];
    for (int i = 0; i < dim; i++) dispersion_tensor[c](i, i) *= vol_phi_ws;
  }
}


/* *******************************************************************
 * Collect time-dependent boundary data in face-based arrays.                               
 ****************************************************************** */
void Transport_PK::ExtractBoundaryConditions(const int component,
                                             std::vector<int>& bc_face_id,
                                             std::vector<double>& bc_face_value)
{
  bc_face_id.assign(nfaces_wghost, 0);

  for (int n = 0; n < bcs.size(); n++) {
    if (component == bcs_tcc_index[n]) {
      for (Amanzi::Iterator bc = bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
        int f = bc->first;
        bc_face_id[f] = TRANSPORT_BC_CONSTANT_TCC;
        bc_face_value[f] = bc->second;
      }
    }
  }
}


/* *******************************************************************
 * Calculate field values at harmonic points. For harmonic points on
 * domain boundary, we use Dirichlet boundary values.
 ****************************************************************** */
void Transport_PK::PopulateHarmonicPointsValues(int component,
                                                Teuchos::RCP<Epetra_MultiVector> tcc,
                                                std::vector<int>& bc_face_id,
                                                std::vector<double>& bc_face_values)
{
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; f++) {
    double weight;
    mfd.CalculateHarmonicPoints(f, dispersion_tensor, harmonic_points[f], weight);
    harmonic_points_weight[f] = weight;

    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    if (ncells == 2) {
      harmonic_points_value[f] = weight * (*tcc)[component][cells[0]]
                               + (1 - weight) * (*tcc)[component][cells[1]];
    } else if (bc_face_id[f] == TRANSPORT_BC_CONSTANT_TCC) {
      harmonic_points_value[f] = bc_face_values[f];
    } else {
      harmonic_points_value[f] = (*tcc)[component][cells[0]];  // ad-hoc solution (lipnikov@lanl.gov)
    }
  }
}


/* *******************************************************************
 * Calculate and add dispersive fluxes of the conservative quantatity.
 ****************************************************************** */
void Transport_PK::AddDispersiveFluxes(int component,
                                       Teuchos::RCP<Epetra_MultiVector> tcc,
                                       std::vector<int>& bc_face_id,
                                       std::vector<double>& bc_face_values,
                                       Teuchos::RCP<Epetra_MultiVector> tcc_next)
{
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List nodes, faces;
  std::vector<AmanziGeometry::Point> corner_points;
  std::vector<double> corner_values, corner_fluxes;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();
    double value = (*tcc)[component][c];

    for (int n = 0; n < nnodes; n++) {
      int v = nodes[n];
      mesh_->node_get_cell_faces(v, c, AmanziMesh::USED, &faces);
      int nfaces = faces.size();

      corner_points.clear();
      corner_values.clear();
      for (int i = 0; i < nfaces; i++) {
        int f = faces[i];
        corner_points.push_back(harmonic_points[f]);
        corner_values.push_back(harmonic_points_value[f]);
      }

      mfd.DispersionCornerFluxes(v, c, dispersion_tensor[c], 
                                 corner_points, value, corner_values, corner_fluxes);

      for (int i = 0; i < nfaces; i++) {
        int f = faces[i];
        if (bc_face_id[f] == TRANSPORT_BC_DISPERSION_FLUX) {
          corner_fluxes[i] = bc_face_values[i];
        }

        (*tcc_next)[component][c] += corner_fluxes[i];
        int c2 = mfd.cell_get_face_adj_cell(c, f);
        if (c2 >= 0) (*tcc_next)[component][c2] -= corner_fluxes[i];
      }
    }
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi



