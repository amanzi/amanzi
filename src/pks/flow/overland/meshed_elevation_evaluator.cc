/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An elevation evaluator getting values from the volumetric mesh.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
#include "Mesh.hh"
#include "Mesh_MSTK.hh"
#include "Point.hh"

#include "meshed_elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

MeshedElevationEvaluator::MeshedElevationEvaluator() :
    ElevationEvaluator() {};

MeshedElevationEvaluator::MeshedElevationEvaluator(const MeshedElevationEvaluator& other) :
    ElevationEvaluator() {};

Teuchos::RCP<FieldEvaluator>
MeshedElevationEvaluator::Clone() const {
  return Teuchos::rcp(new MeshedElevationEvaluator(*this));
}

void MeshedElevationEvaluator::EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {

  Teuchos::Ptr<CompositeVector> elev = results[0];
  Teuchos::Ptr<CompositeVector> slope = results[1];

  // Get the elevation and slope values from the domain mesh.
  Teuchos::RCP<const AmanziMesh::Mesh_MSTK> domain_mesh =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S->GetMesh());
      // Note that static cast is safe here because we have
      // already ensured it was MSTK.

  Teuchos::RCP<const AmanziMesh::Mesh_MSTK> surface_mesh =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S->GetMesh("surface"));
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

  elev->ScatterMasterToGhosted();
  slope->ScatterMasterToGhosted();

}

} //namespace
} //namespace
} //namespace
