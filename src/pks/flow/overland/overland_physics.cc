/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Gianmarco Manzini
         Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Mesh.hh"
#include "Mesh_MSTK.hh"

#include "overland.hh"

namespace Amanzi {
namespace Flow {

#if 0
}}
#endif

// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void OverlandFlow::ApplyDiffusion_(const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<CompositeVector>& g) {

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> upwind_conductivity =
    S->GetFieldData("upwind_overland_conductivity", "overland_flow");

  // update the stiffness matrix
  matrix_->CreateMFDstiffnessMatrices(*upwind_conductivity);
  matrix_->CreateMFDrhsVectors();

  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData("pres_elev");
  matrix_->ComputeNegativeResidual(*pres_elev, g);
};


// -------------------------------------------------------------
// Accumulation of internal energy term du/dt
// -------------------------------------------------------------
void OverlandFlow::AddAccumulation_(const Teuchos::RCP<CompositeVector>& g) {
  const CompositeVector & pres0        = *(S_inter_->GetFieldData("overland_pressure"));
  const CompositeVector & pres1        = *(S_next_ ->GetFieldData("overland_pressure"));
  const CompositeVector & cell_volume0 = *(S_inter_->GetFieldData("surface_cell_volume"));
  const CompositeVector & cell_volume1 = *(S_next_ ->GetFieldData("surface_cell_volume"));

  double dt = S_next_->time() - S_inter_->time();

  int c_owned = g->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    (*g)("cell",0,c) += (cell_volume1("cell",c)*pres1("cell",c)
                         - cell_volume0("cell",c)*pres0("cell",c))/dt ;
  }
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void OverlandFlow::AddLoadValue_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& g) {

  Teuchos::RCP <const CompositeVector> cell_volume = S->GetFieldData("surface_cell_volume");
  Teuchos::RCP <const CompositeVector> rain = S->GetFieldData("rainfall_rate");

  int c_owned = g->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    (*g)("cell",c) -= (*rain)("cell",c) * (*cell_volume)("cell",c);
  }
};


// -----------------------------------------------------------------------------
// Update variables, like p + z and rel perm.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateSecondaryVariables_(const Teuchos::RCP<State>& S) {
  // add elevation
  AddElevation_(S) ;

  // compute non-linear coefficient
  UpdateConductivity_(S);

  // update rainfall
  UpdateLoadValue_(S);
};


// -----------------------------------------------------------------------------
// Evaluate model for elevation + pressure
// -----------------------------------------------------------------------------
void OverlandFlow::AddElevation_(const Teuchos::RCP<State>& S) {

  const CompositeVector & pres      = *(S->GetFieldData("overland_pressure"));
  const CompositeVector & elevation = *(S->GetFieldData("elevation"));
  CompositeVector       & pres_elev = *(S->GetFieldData("pres_elev","overland_flow"));

  // add elevation to pres_elev, cell dofs
  int c_owned = pres_elev.size("cell");
  for (int c=0; c!=c_owned; ++c) {
    pres_elev("cell",c) = pres("cell",c) + elevation("cell",c) ;
  }

  // add elevation to pres_elev, face dofs
  int f_owned = pres_elev.size("face");
  for (int f=0; f!=f_owned; ++f) {
    pres_elev("face",f) = pres("face",f) + elevation("face",f) ;
  }
}


// -----------------------------------------------------------------------------
// Update the rainfall rate.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateLoadValue_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<CompositeVector> rain_rate =
    S->GetFieldData("rainfall_rate", "overland_flow");
  rain_rate_function_->Compute(S->time(), rain_rate.ptr());
};


// -----------------------------------------------------------------------------
// Update the conductivity.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateConductivity_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("overland_pressure");
  Teuchos::RCP<const CompositeVector> mann = S->GetFieldData("manning_coef");
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData("slope_magnitude");

  Teuchos::RCP<CompositeVector> cond =
    S->GetFieldData("overland_conductivity", "overland_flow");

  Conductivity_(S, *pres, *mann, *slope, cond);
};


// -----------------------------------------------------------------------------
// Update the "relative permeability".
// -----------------------------------------------------------------------------
void OverlandFlow::Conductivity_(const Teuchos::RCP<State>& S,
                                 const CompositeVector& pressure,
                                 const CompositeVector& manning_coef,
                                 const CompositeVector& slope_mag,
                                 const Teuchos::RCP<CompositeVector>& conductivity) {
  double eps = 1.e-14;
  double exponent = 1.0 + manning_exp_;

  int ncells = conductivity->size("cell");
  for (int c=0; c!=ncells; ++c ) {
    double scaling = manning_coef("cell",c) * std::sqrt(slope_mag("cell",c) + eps);
    (*conductivity)("cell",c) = std::pow(pressure("cell",c), exponent) / scaling;
  }
};


// -----------------------------------------------------------------------------
// Update the elevation of the surface mesh and its slope.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateElevationAndSlope_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<CompositeVector> elev = S->GetFieldData("elevation", "overland_flow");
  Teuchos::RCP<CompositeVector> slope = S->GetFieldData("slope_magnitude", "overland_flow");

  if (standalone_mode_) {
    // No domain mesh, so get the elevation and slope values from functions.
    elevation_function_->Compute(S->time(), elev.ptr());
    slope_function_->Compute(S->time(), slope.ptr());

  } else {
    // Get the elevation and slope values from the domain mesh.
    Teuchos::RCP<const AmanziMesh::Mesh_MSTK> domain_mesh =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S->Mesh());
                    // Note that static cast is safe here because we have
                    // already ensured it was MSTK.

    Teuchos::RCP<const AmanziMesh::Mesh_MSTK> surface_mesh =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S->Mesh("surface"));
                    // Note that static cast is safe here because we have
                    // already ensured it was MSTK.

    if (domain_mesh->cell_dimension() == 3) {
      // Set the elevation on cells by getting the corresponding face and its
      // centroid.
      int ncells = elev->size("cell");
      for (int c=0; c!=ncells; ++c) {
        // Note that a surface cell is a volume mesh's face
        AmanziMesh::Entity_ID domain_face =
          surface_mesh->entity_get_parent(AmanziMesh::CELL, c);

        // First elevation.
        AmanziGeometry::Point x = domain_mesh->face_centroid(domain_face);
        (*elev)("cell",c) = x[2];

        // Now slope.
        AmanziGeometry::Point n = domain_mesh->face_normal(domain_face);
        // -- S = || n - (n dot z) z || / | n dot z |
        (*slope)("cell",c) = std::sqrt( std::pow(n[0],2) + std::pow(n[1],2))
                                          / std::abs(n[2]);
      }

      // Set the elevation on faces by getting the corresponding nodes and
      // averaging.
      AmanziMesh::Entity_ID_List surface_nodes(2);
      AmanziMesh::Entity_ID node0, node1;
      AmanziGeometry::Point coord0(3);
      AmanziGeometry::Point coord1(3);

      int nfaces = elev->size("face");
      for (int f=0; f!=nfaces; ++f) {
        surface_mesh->face_get_nodes(f, &surface_nodes);
        node0 = surface_mesh->entity_get_parent(AmanziMesh::NODE, surface_nodes[0]);
        node1 = surface_mesh->entity_get_parent(AmanziMesh::NODE, surface_nodes[1]);

        domain_mesh->node_get_coordinates(node0, &coord0);
        domain_mesh->node_get_coordinates(node1, &coord1);
        (*elev)("face",f) = (coord0[2] + coord1[2])/2.0;
      }
    } else if (domain_mesh->cell_dimension() == 2) {
      // Set the elevation on cells by getting the corresponding cell and its
      // centroid.
      int ncells = elev->size("cell");
      for (int c=0; c!=ncells; ++c) {
        // Note that a surface cell is a surface mesh's cell.
        AmanziMesh::Entity_ID domain_cell =
          surface_mesh->entity_get_parent(AmanziMesh::CELL, c);

        // First elevation.
        AmanziGeometry::Point x = domain_mesh->cell_centroid(domain_cell);
        (*elev)("cell",c) = x[2];

        // Now slope.
        AmanziMesh::Entity_ID_List faces;
        domain_mesh->cell_get_faces(domain_cell, &faces);

        // -- Get the normals of all faces of the surface cell.
        int count = faces.size();
        std::vector<AmanziGeometry::Point> normals(count);
        for (int lcv=0; lcv!=count; ++lcv) {
          normals[lcv] = domain_mesh->face_normal(faces[lcv], false, domain_cell);
        }

        // -- Average the cross product of successive faces to get a cell normal.
        AmanziGeometry::Point cross(0.0, 0.0, 0.0);
        for (int lcv=0; lcv!=(count-1); ++lcv) {
          cross[0] += normals[lcv][1]*normals[lcv+1][2] - normals[lcv][2]*normals[lcv+1][1];
          cross[1] += normals[lcv][2]*normals[lcv+1][0] - normals[lcv][0]*normals[lcv+1][2];
          cross[2] += normals[lcv][0]*normals[lcv+1][1] - normals[lcv][1]*normals[lcv+1][0];
        }
        cross[0] += normals[count][1]*normals[0][2] - normals[count][2]*normals[0][1];
        cross[1] += normals[count][2]*normals[0][0] - normals[count][0]*normals[0][2];
        cross[2] += normals[count][0]*normals[0][1] - normals[count][1]*normals[0][0];
        cross /= count;

        // -- S = || n - (n dot z) z || / | n dot z |
        (*slope)("cell",c) = std::sqrt( std::pow(cross[0],2) + std::pow(cross[1],2))
                                          / std::abs(cross[2]);
      }

      // Set the elevation on faces by getting the corresponding face and its
      // centroid.
      int nfaces = elev->size("face");
      for (int f=0; f!=nfaces; ++f) {
        // Note that a surface face is a surface mesh's face.
        AmanziMesh::Entity_ID domain_face =
          surface_mesh->entity_get_parent(AmanziMesh::FACE, f);
        AmanziGeometry::Point x = domain_mesh->face_centroid(domain_face);

        (*elev)("face",f) = x[2];
      }
    }
  }

  elev->ScatterMasterToGhosted();
  slope->ScatterMasterToGhosted();
};


} //namespace
} //namespace
