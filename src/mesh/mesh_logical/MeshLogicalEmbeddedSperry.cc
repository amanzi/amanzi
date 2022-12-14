/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! A factory for embedding Sperry-style root meshes within a standard 3D mesh.
/*!

This factory looks to embed Sperry-style root meshes within a 3D standard
domain subsurface mesh.  It optionally builds a logical mesh (for sequential
coupling) or an embedded mesh (for global implicit coupling).

*/

#include "MeshLogicalEmbeddedSperry.hh"

namespace Amanzi {
namespace AmanziMesh {


Teuchos::RCP<MeshLogical>
MeshLogicalEmbeddedSperry::CreateLogical(const std::string& name,
                                         int n_stem,
                                         double max_rooting_depth,
                                         int n_rheizosphere_shells)
{
  name_ = name;
  n_stem_ = n_stem;
  max_rooting_depth_ = max_rooting_depth;
  n_rheizosphere_shells_ = n_rheizosphere_shells;

  sperry_ = CreateLogical_();
  return sperry_;
}


Teuchos::RCP<MeshLogical>
MeshLogicalEmbeddedSperry::CreateEmbedded(const std::string& name,
                                          int n_stem,
                                          double max_rooting_depth,
                                          int n_rheizosphere_shells,
                                          const Teuchos::RCP<AmanziMesh::Mesh>& bg_mesh)
{
  name_ = name;
  n_stem_ = n_stem;
  max_rooting_depth_ = max_rooting_depth;
  n_rheizosphere_shells_ = n_rheizosphere_shells;

  sperry_ = CreateLogical_();
  return CreateEmbedded_();
}


Teuchos::RCP<MeshLogical>
MeshLogicalEmbeddedSperry::CreateLogical_(const std::string& pftname, int column) const
{
  bg_mesh_->build_columns();
  const auto& cells_of_col = bg_mesh_->cells_of_column(column);
  const auto& faces_of_col = bg_mesh_->faces_of_column(column);

  // I never remember which is up, 0 or -1?  This checks that it is 0.
  auto ground_surface = bg_mesh_->face_centroid(faces_of_col[0]);
  AMANZI_ASSERT(ground_surface[2] > bg_mesh_->cell_centroid(cells_of_col[0])[2]);

  // find how many cells in the rooting zone
  double root_bottom_z = ground_surface[2] - max_rooting_depth_;
  int n_cells_root_zone = 0;
  for (auto f : faces_of_col) {
    if (root_bottom_z >= bg_mesh_->face_centroid(f)[2]) break;
    n_cells_root_zone++;
  }

  // add the leaf
  AmanziGemoetry::Point first(ground_surface);
  AmanziGemoetry::Point last(ground_surface);
  first[2] += 3.1;
  last[2] += 3.;

  AmanziMesh::Entity_ID_List leaf_cell, leaf_face;
  AddSegment(1,
             first,
             last,
             1.0,
             MeshLogicalFactory::LogicalTip_t::BOUNDARY,
             MeshLogicalFactory::LogicalTip_t::JUNCTION,
             Name_(pftname, column, "leaf"),
             &leaf_cell,
             &leaf_face);

  // reserve the leaf-stem face
  Entity_ID leaf_stem_face = ReserveFace();

  // add the leaf-ground segment
  first = last;
  last = ground_surface;
  AmanziMesh::Entity_ID_List stem_cells, stem_faces;
  AddSegment(n_stem,
             first,
             last,
             1.0,
             MeshLogicalFactory::LogicalTip_t::BRANCH,
             MeshLogicalFactory::LogicalTip_t::JUNCTION,
             Name_(pftname, column, "stem"),
             &stem_cells,
             &stem_faces);

  // add the leaf-stem face
  AddFace(leaf_stem_face,
          Entity_ID_List{ leaf_cell[0], stem_cells[0] },
          std::vector<double>{ 0.5, 3.0 / n_stem_ / 2.0 },
          1.0);

  // reserve the stem-root face
  Entity_ID stem_root_face = ReserveFace();

  // the transporting root needs to coincide with the soil spacing
  first = ground_surface;
  last[2] -= max_rooting_depth_;
  std::vector<double> cell_dzs;
  for (int cc = 0; cc != n_cells_root_zone; ++cc) {
    cells_dzs.push_back(bg_mesh_->face_centroid(faces_of_col[cc])[2] -
                        bg_mesh_->face_centroid(faces_of_col[cc + 1])[2]);
  }
  cells_dzs.back() = bg_mesh_->face_centroid(faces_of_col[n_cells_root_zone - 1])[2] -
                     (ground_surface[2] - max_rooting_depth_);
  AMANZI_ASSERT(cells_dzs.back() > 0.);

  // centroids are cell centroids, the bottom one is close enough
  std::vector<AmanziGemoetry::Point> centroids;
  for (int cc = 0; cc != n_cells_root_zone; ++cc) {
    centroids.push_back(bg_mesh_ - cell_centroid(cells_of_col[cc]));
  }

  std::vector<double> face_areas(n_cells_root_zone - 1, 1.);
  AmanziMesh::Entity_ID_List troot_cells, troot_faces;
  AddSegment(centroids,
             nullptr,
             cells_dzs,
             face_areas,
             bottom - top,
             MeshLogicalFactory::LogicalTip_t::BRANCH,
             MeshLogicalFactory::LogicalTip_t::JUNCTION,
             Name_(pftname, column, "troot"),
             &troot_cells,
             &troot_faces);

  // add the stem-to-troot face
  AddFace(stem_root_face,
          Entity_ID_List{ stem_cells.back(), troot_cells.front() },
          std::vector<double>{ 3.0 / n_stem_ / 2.0, cells_dzs[0] / 2. },
          1.0);

  // add absorbing roots and faces
  for (int cc = 0; cc != n_cells_root_zone; ++cc) {
    // reserve the aroot-troot face
    Entity_ID aroot_troot = ReserveFace();

    // add the aroot cell
    AmanziMesh::Entity_ID_List aroot_cell, aroot_faces;

    // aroot centroid
    begin = centroids[cc];
    end = centroids[cc];
    end[0] += 0.01;
    AddSegment(1,
               begin,
               end,
               1.0,
               MeshLogicalFactory::LogicalTip_t::BRANCH,
               MeshLogicalFactory::LogicalTip_t::JUNCTION,
               Name_(pftname, column, "aroot", cc),
               &aroot_cell,
               &aroot_face);

    // add the aroot-troot face
    AddFace(aroot_troot,
            Entity_ID_List{ troot_cells[cc], aroot_cell[0] },
            std::vector<double>{ 1., 1. },
            1.0);

    // reserve the soil-aroot face
    Entity_ID soil_aroot = ReserveFace();

    // add the rheizosphere shells
    begin = end;
    end[0] += 0.1;
    AmanziMesh::Entity_ID_List shell_cells, shell_faces;
    AddSegment(n_rheizosphere_shells,
               begin,
               end,
               1.0,
               MeshLogicalFactory::LogicalTip_t::BRANCH,
               MeshLogicalFactory::LogicalTip_t::JUNCTION,
               Name_(pftname, column, "rheizosphere", cc),
               &shell_cells,
               &shell_faces);

    // add the shell-to-aroot faces
    AddFace(soil_aroot,
            Entity_ID_List{ aroot_cell[0], shell_cells.front() },
            std::vector<double>{ 1., 1. },
            1.);

    // store info to connect rheizosphere_cell to background cell
    rheizosphere_to_bg_.emplace_back(Entity_ID_List{ shell_cells.back(), cells_in_col[cc] });
  }
};

void
MeshLogicalEmbeddedSperry::set_largest_shell_diameter(double shell_diameter)
{}


void
MeshLogicalEmbeddedSperry::set_stem_height(double height)
{
  // update centroids
}

//  ??????
// other geometric options from traits


} // namespace AmanziMesh
} // namespace Amanzi

#endif
