#include <cmath>
#include "boost/math/constants/constants.hpp"

#include "MeshLogicalFactory.hh"

#include "plant_1D_mesh.hh"

namespace Amanzi {
namespace Testing {

const double PI = boost::math::constants::pi<double>();

using namespace Amanzi::AmanziMesh;

// Creates a plant mesh
Teuchos::RCP<AmanziMesh::MeshLogical>
plantMesh(const Epetra_MpiComm* comm,
          const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
          bool include_soil) {

  
  AmanziMesh::MeshLogicalFactory fac(comm, gm);

  int n_stem = 10;
  double tree_height = 3.;
  int n_root = 10;
  double root_depth = 2.;
  double stem_cap_radius = 0.01; // effective cappillary radius
  double root_cap_radius = 0.01; // effective cappillary radius
  double root_radius = 0.1; // physical root radius
  double soil_depth = 10.;

  double aroot_radius = 0.01; // physical thickness of the a-root shell

  // mean distance from the soil cylinder to the centroid of the root
  double mean_soil_to_centroid_length = 2 * PI * std::pow(1.0,3) / 3.0
      / (2*PI*std::pow(1,3));
  
  
  // add the leafs
  AmanziGeometry::Point top(3), bottom(3);
  top.set(0.); bottom.set(0.);
  top[1] = mean_soil_to_centroid_length;
  top[2] = tree_height + 0.1;
  bottom[1] = top[1];
  bottom[2] = tree_height;
  double leaf_area = PI * std::pow(stem_cap_radius,2);
  
  std::vector<AmanziMesh::Entity_ID> leaf_cells, leaf_faces;
  double leaf_length;
  fac.AddSegment(1, top, bottom, leaf_area,
                 MeshLogicalFactory::TIP_BOUNDARY,
                 MeshLogicalFactory::TIP_DEFERRED,
                 "leaf",
                 &leaf_cells, &leaf_faces, &leaf_length);
  ASSERT(std::abs(leaf_length - 0.1) < 1.e-8);

  // add the leaf top bc
  AmanziMesh::Entity_ID_List leaf_tip;
  leaf_tip.push_back(leaf_faces[0]);
  fac.AddSet("leaf tip", "face", leaf_tip);
  
  // add the stem
  double dz_stem = tree_height / n_stem;
  top[2] = tree_height;
  bottom[2] = 0.0;
  double stem_area = PI * std::pow(stem_cap_radius, 2);

  std::vector<AmanziMesh::Entity_ID> stem_cells, stem_faces;
  double stem_length;
  fac.AddSegment(n_stem, top, bottom, stem_area,
                 MeshLogicalFactory::TIP_DEFERRED,
                 MeshLogicalFactory::TIP_DEFERRED,
                 "stem", &stem_cells, &stem_faces, &stem_length);
  ASSERT(std::abs(stem_length - dz_stem) < 1.e-8);

  // add the transporting root
  double dz_root = root_depth / n_root;
  top[2] = 0.0;
  bottom[2] = -root_depth;
  double troot_area = PI * std::pow(root_cap_radius, 2);

  std::vector<AmanziMesh::Entity_ID> troot_cells, troot_faces;
  double troot_length;
  fac.AddSegment(n_root, top, bottom, troot_area,
                 MeshLogicalFactory::TIP_DEFERRED,
                 MeshLogicalFactory::TIP_DEFERRED,
                 "troot", &troot_cells, &troot_faces, &troot_length);

  // add the absorbing roots, each their own segment
  // NOTE: to do more than one, this needs variable areas
  //
  // The absorbing root is a set of cylindrical shells around the troot
  // connection troot to soil.
  std::vector<std::vector<AmanziMesh::Entity_ID> > aroot_cells;
  aroot_cells.resize(n_root);
  int n_aroot = 1;
  AmanziMesh::Entity_ID_List aroots;
  top.set(0.); bottom.set(0.);

  for (int i=0; i!=n_root; ++i) {
    top[2] = -(i+0.5) *  dz_root;
    bottom[2] = top[2];

    top[1] = mean_soil_to_centroid_length - root_radius; // physical radius!
    bottom[1] = top[1] - aroot_radius;

    double aroot_area = troot_length * PI
        * std::pow(root_radius + 0.5*aroot_radius, 2);
    std::vector<AmanziMesh::Entity_ID> aroot_faces;
    double aroot_length;
    std::stringstream setname;
    setname << "aroot_" << i;
    fac.AddSegment(n_aroot, top, bottom, aroot_area,
                   MeshLogicalFactory::TIP_DEFERRED,
                   MeshLogicalFactory::TIP_DEFERRED,
                   setname.str(), &aroot_cells[i], &aroot_faces, &aroot_length);
    aroots.insert(aroots.end(), aroot_cells[i].begin(), aroot_cells[i].end());
  }
  fac.AddSet("aroot", "cell", aroots);

  std::vector<AmanziMesh::Entity_ID> soil_cells, soil_faces;
  if (include_soil) {
    // add the soil columns
    double dz_soil = dz_root;
    int n_soil = (int)std::round(soil_depth / dz_soil);
    ASSERT(std::abs(n_soil * dz_soil - soil_depth) < 1.e-8);
    top.set(0.);
    bottom.set(0.);
    bottom[2] = -soil_depth;
    double soil_area = 2*PI*std::pow(1,2); // 1m cylinder
    double soil_length;
    fac.AddSegment(n_soil, top, bottom, soil_area,
                   MeshLogicalFactory::TIP_BOUNDARY,
                   MeshLogicalFactory::TIP_BOUNDARY,
                   "soil", &soil_cells, &soil_faces, &soil_length);

    // add the soil bcs
    AmanziMesh::Entity_ID_List soil_bc;
    soil_bc.push_back(soil_faces[0]);
    fac.AddSet("soil surface", "face", soil_bc);
    soil_bc[0] = soil_faces[n_soil-1];
    fac.AddSet("soil bottom", "face", soil_bc);
  }
  
  // Connect the segments:
  // -- connect leaf to stem
  Entity_ID_List stem_to_leaf(2);
  AmanziGeometry::Point normal(3);
  std::vector<double> stem_to_leaf_lengths(2);

  stem_to_leaf[0] = leaf_cells[0];
  stem_to_leaf_lengths[0] = leaf_length / 2.0;

  stem_to_leaf[1] = stem_cells[0];
  stem_to_leaf_lengths[1] = stem_length / 2.0;
  normal[2] = -1.;

  int iconn = fac.AddConnection(stem_to_leaf, normal,
          stem_to_leaf_lengths, stem_area);

  // -- connect stem to troot
  Entity_ID_List stem_to_troot(2);
  std::vector<double> stem_to_troot_lengths(2);

  stem_to_troot[0] = stem_cells[n_stem-1];
  stem_to_troot_lengths[0] = stem_length / 2.0;

  stem_to_troot[1] = troot_cells[0];
  stem_to_troot_lengths[1] = troot_length / 2.0;
  normal[2] = -1.;

  iconn = fac.AddConnection(stem_to_troot, normal,
          stem_to_troot_lengths, (stem_area+troot_area) / 2.);

  // -- connect troot to aroot, aroot to soil
  for (int i=0; i!=n_root; ++i) {
    // troot to aroot
    Entity_ID_List troot_to_aroot(2);
    std::vector<double> troot_to_aroot_lengths(2);

    troot_to_aroot[0] = troot_cells[i];
    troot_to_aroot_lengths[0] = root_radius;

    troot_to_aroot[1] = aroot_cells[i][0];
    troot_to_aroot_lengths[1] = aroot_radius / 2.;
    normal.set(0.); normal[1] = 1.;
    double troot_to_aroot_area = dz_root *
        PI * std::pow(root_radius, 2);        

    iconn = fac.AddConnection(troot_to_aroot, normal,
            troot_to_aroot_lengths, troot_to_aroot_area);

    if (include_soil) {
      // aroot to soil
      Entity_ID_List aroot_to_soil(2);
      std::vector<double> aroot_to_soil_lengths(2);

      aroot_to_soil[0] = aroot_cells[i][n_aroot-1];
      aroot_to_soil_lengths[0] = aroot_radius / 2.;

      aroot_to_soil[1] = soil_cells[i];
      aroot_to_soil_lengths[1] = mean_soil_to_centroid_length
          - root_radius - aroot_radius;  // this is (more or less) the mean
      // distance to the centroid of the soil cell?  Seems reasonable...
      normal.set(0.); normal[1] = 1.;
      double aroot_to_soil_area = dz_root *
          PI * std::pow(root_radius+aroot_radius, 2);        

      iconn = fac.AddConnection(aroot_to_soil, normal,
              aroot_to_soil_lengths, aroot_to_soil_area);
    }
  }

  return fac.Create();
}


} // namespace
} // namespace
