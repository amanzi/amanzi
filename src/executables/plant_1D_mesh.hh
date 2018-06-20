//
// Set of functions that create simple plant meshes.
//
// ------------------------------------------------------------------

#ifndef MESH_LOGICAL_PLANT_MESHES_HH_
#define MESH_LOGICAL_PLANT_MESHES_HH_

#include <mpi.h>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

#include "RegionEnumerated.hh"
#include "MeshLogical.hh"
#include "MeshLogicalFactory.hh"
#include "MeshEmbeddedLogical.hh"

#include "PFT.hh"



namespace ATS {
namespace Testing {

using namespace Amanzi;

struct HydraulicMeshParams {
  HydraulicMeshParams() :
      eps_(1.e-10) {}

  HydraulicMeshParams(const std::string& name_) :
      name(name_),
      eps_(1.e-10) {}
  
  
  std::string name;             // name of the PFT/cohort
  double height;                // full-plant height [m]

  double leaf_height;           // height of the canopy leaf layer [m]  (< height)
  double leaf_area;             // total leaf area [m^2]
  double leaf_volume;           // total leaf volume [m^3]

  double stem_length;           // stem length [m]  (stem_length + leaf_height = height)
  double sapwood_area;          // total cross-sectional sapwood area at breast height [m^2]

  double troot_volume;          // transporting root volume [m^3]
  double troot_length;          // typical (radial) length of transporting root [m]

  double aroot_radius;          // typical distance from soil-root interface to absorbing root [m]
  double aroot_length;          // total length of absorbing roots [m] (typically root-soil interface area / (2*PI*aroot_radius))
  double aroot_volume;          // volume of absorbing root [m^3] (typically aroot_length * PI * aroot_radius^2)

  double eps_;
  
  bool Valid() const {
    if (height <= eps_) return false;
    if (leaf_height <= eps_) return false;
    if (leaf_volume <= eps_) return false;
    if (stem_length <= eps_) return false;
    if (std::abs(stem_length + leaf_height - height) > eps_) return false;
    if (sapwood_area <= eps_) return false;
    if (troot_volume <= eps_) return false;
    if (troot_length <= eps_) return false;
    if (aroot_radius <= eps_) return false;
    if (aroot_length <= eps_) return false;
    if (aroot_volume <= eps_) return false;
    return true;
  }
};

Teuchos::RCP<AmanziMesh::Mesh> 
plantMeshes(const Teuchos::RCP<AmanziMesh::Mesh>& subsurface_mesh,
            const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
            const AmanziMesh::Mesh& surface_mesh,
            const Epetra_MultiVector& rooting_fraction,
            const std::vector<std::vector<BGC::PFT>>& pfts,
            int n_stem, int n_rheiz_shell,
            double rheiz_shell_fraction);

void
plantMesh(AmanziMesh::MeshLogicalFactory& fac,
          const HydraulicMeshParams& params,
          int n_stem,
          int n_rheiz_shell,
          const std::vector<double> rooting_fraction,
          const std::vector<double>& soil_col_dz,
          const AmanziGeometry::Point& soil_col_centroid,
          const AmanziGeometry::Point& stem_displacement,
          const std::vector<double>& rheiz_outer_radius,
          std::vector<std::vector<int>>& soil_rheiz_conn_cells,
          std::vector<std::vector<double>>& soil_rheiz_conn_lengths,
          std::vector<AmanziGeometry::Point>& soil_rheiz_conn_area_normals,
          std::vector<double>& belowground_volume);

} // namespace
} // namespace

#endif
