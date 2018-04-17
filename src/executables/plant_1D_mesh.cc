#include <cmath>
#include "boost/math/constants/constants.hpp"
#include "boost/math/special_functions/round.hpp"

#include "MeshLogicalAudit.hh"
#include "MeshLogicalFactory.hh"

#include "plant_1D_mesh.hh"

/*

volume of:

leaf volume: leaf biomass / density
leaf length: min(0.1, 0.9*height)
leaf area: sapwood area?

stem length: (height - leaf_length) / nstem
stem volume: sapwood area * length  (does not account for taper at this point!)
stem area: sapwood area

troot volume:  belowground biomass --(?? woody carbon bg)--> troot biomass --(via wood density)--> belowground volume --(via rooting fraction)--> troot volume
troot length:  above/below = dz, to aroot = volume / (sapwood area * rooting fraction)?
troot area: accumulated rooting fraction * sapwood area (above/below),  ?? (to aroot)

aroot length: br (fine root biomass) * hydr_srl
aroot volume: aroot length * (PI * hydr_rs2 **2 )  NOTE ARE THESE CONSISTENT?  Probably not!  Maybe should be br / density
aroot area: ?? (to troot), 

rheiz shell length: (radius_outer - radius_inner)
rheiz shell area: aroot area
rheiz shell volume: length * (area_outer - area_inner)


rheiz_shell radius outer must be calculated such that volume of all total rheiz shells is less than the volume of the soil column

*/


namespace ATS {
namespace Testing {

const double PI = boost::math::constants::pi<double>();

using namespace Amanzi;

Teuchos::RCP<AmanziMesh::Mesh> 
plantMeshes(const Teuchos::RCP<AmanziMesh::Mesh>& subsurface_mesh,
            const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
            const AmanziMesh::Mesh& surface_mesh,
            const Epetra_MultiVector& rooting_fraction,
            const std::vector<std::vector<BGC::PFT>>& pfts,
            int n_stem, int n_rheiz_shell,
            double rheiz_shell_fraction)
{
  // plant 1 mesh per PFT per surface grid cell
  AmanziMesh::MeshLogicalFactory fac(subsurface_mesh->get_comm(),gm,false);

  std::vector<std::vector<int>> soil_rheiz_conn_cells;
  std::vector<std::vector<double>> soil_rheiz_conn_lengths;
  std::vector<AmanziGeometry::Point> soil_rheiz_conn_area_normals;

  int n_surf_cells = surface_mesh.num_entities(AmanziMesh::CELL,AmanziMesh::Parallel_type::OWNED);
  for (int sc=0; sc!=n_surf_cells; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh.entity_get_parent(AmanziMesh::CELL, sc);
    auto soil_col_centroid = subsurface_mesh->face_centroid(f);

    auto& cells_of_col = subsurface_mesh->cells_of_column(sc);
    auto& faces_of_col = subsurface_mesh->faces_of_column(sc);

    std::vector<double> soil_col_dz(cells_of_col.size());
    for (int cc=0; cc!=cells_of_col.size(); ++cc) {
      soil_col_dz[cc] = subsurface_mesh->face_centroid(faces_of_col[cc])[2]
                        - subsurface_mesh->face_centroid(faces_of_col[cc+1])[2];
    }
    
    std::vector<double> volume_shells(cells_of_col.size(), 0.);

    double rotation = 2 * PI / pfts[sc].size(); // rotate roots for visualization
    double surface_area = surface_mesh.cell_volume(sc);
    double shift_size = std::sqrt(surface_area/PI)/2.;

    // set the geometric parameters
    std::vector<HydraulicMeshParams> hps;
    for (auto pft : pfts[sc]) {
      HydraulicMeshParams hp(pft.pft_type+"_"+std::to_string(sc));
      hp.height = pft.height;
      hp.leaf_height = std::min(0.1, 0.9*pft.height);
      hp.stem_length = hp.height - hp.leaf_height;
      hp.leaf_area = surface_area * pft.lai;
      
      double leaf_density = -2.3231*pft.SLA / pft.carbon2biomass + 781.899;
          // from FATES, empirical regression data from leaves at Caxiuana (~ 8 spp)
      hp.leaf_volume = surface_area * pft.Bleaf * pft.carbon2biomass  / leaf_density;

      double leaf_to_sapwood_area_ratio = 1.e4 * std::pow(1e3 / pft.SLA , -2.14)*pft.height*546;
          // from BChristoffersen paper, Al:As [m^2 cm^-2] = 546 * H * LMA [g m^-2]^-2.14, converted back to m^2 m^-2.
      hp.sapwood_area = hp.leaf_area / leaf_to_sapwood_area_ratio;

      hp.troot_length = hp.height / 4.;  // made up!
      hp.troot_volume = surface_area * 0.4 * pft.Bstem * pft.carbon2biomass / pft.wood_density;
          // made up, 0.4 * stem biomass?  Need some measure of belowground vs aboveground wood
      
      double aroot_absorbing_area = 1 * hp.leaf_area; // from BChristoffersen paper, leaf area : absorbing area = 1
      hp.aroot_radius = pft.aroot_radius;
      hp.aroot_length = aroot_absorbing_area / (2 * PI * hp.aroot_radius);
      hp.aroot_volume = hp.aroot_length * PI * std::pow(hp.aroot_radius,2);

      hps.emplace_back(hp);
    }

    // set the rheizosphere shell radius so that the volume of that soil is about rheiz_shell_fraction of the total subsurface soil volume.
    std::vector<double> rheiz_outer_radius(cells_of_col.size(),0.);
    std::vector<double> total_aroot_length(cells_of_col.size(), 0.);
    if (n_rheiz_shell > 0) {
      for (int i_pft=0; i_pft!=pfts[sc].size(); ++i_pft) {
        for (int cc=0; cc!=cells_of_col.size(); ++cc) {
          total_aroot_length[cc] += rooting_fraction[i_pft][cells_of_col[cc]] * hps[i_pft].aroot_length;
        }
      }

      for (int cc=0; cc!=cells_of_col.size(); ++cc) {
        if (total_aroot_length[cc] > 0.) {
          rheiz_outer_radius[cc] = std::sqrt(subsurface_mesh->cell_volume(cells_of_col[cc]) * rheiz_shell_fraction / total_aroot_length[cc] / PI);
        }
      }
    }

    for (int i_pft=0; i_pft!=pfts[sc].size(); ++i_pft) {
      // shift the stem for vis
      AmanziGeometry::Point stem_displacement(0.,0.,0.);
      stem_displacement[0] += shift_size * std::cos(rotation * i_pft);
      stem_displacement[1] += shift_size * std::sin(rotation * i_pft);
      
      std::vector<double> col_rooting_fraction(cells_of_col.size());
      for (int cc=0; cc!=cells_of_col.size(); ++cc) {
        col_rooting_fraction[cc] = rooting_fraction[i_pft][cells_of_col[cc]];
      }

      std::vector<std::vector<int>> l_soil_rheiz_conn_cells;
      std::vector<std::vector<double>> l_soil_rheiz_conn_lengths;
      std::vector<AmanziGeometry::Point> l_soil_rheiz_conn_area_normals;
      std::vector<double> l_volume_shells(cells_of_col.size());
      
      plantMesh(fac, hps[i_pft], n_stem, n_rheiz_shell,
                col_rooting_fraction, soil_col_dz, soil_col_centroid, stem_displacement,
                rheiz_outer_radius,
                l_soil_rheiz_conn_cells, l_soil_rheiz_conn_lengths,
                l_soil_rheiz_conn_area_normals, l_volume_shells);

      for (int cc=0; cc!=cells_of_col.size(); ++cc) volume_shells[cc] += l_volume_shells[cc];

      // update the connections to point to the background mesh's cells
      // Note background mesh ID must come second so embedding can know how to reset IDs
      for (int cc=0; cc!=l_soil_rheiz_conn_cells.size(); ++cc) {
        l_soil_rheiz_conn_cells[cc][1] = cells_of_col[cc];
        l_soil_rheiz_conn_lengths[cc][1] = rheiz_outer_radius[cc];          // totally made up!
      }

      // save these for embedded mesh
      soil_rheiz_conn_cells.insert(soil_rheiz_conn_cells.end(),
              l_soil_rheiz_conn_cells.begin(), l_soil_rheiz_conn_cells.end());
      soil_rheiz_conn_lengths.insert(soil_rheiz_conn_lengths.end(),
              l_soil_rheiz_conn_lengths.begin(), l_soil_rheiz_conn_lengths.end());
      soil_rheiz_conn_area_normals.insert(soil_rheiz_conn_area_normals.end(),
              l_soil_rheiz_conn_area_normals.begin(), l_soil_rheiz_conn_area_normals.end());
    }

    // check the the volume!
    for (int cc=0; cc!=cells_of_col.size(); ++cc) {
      if (volume_shells[cc] > rheiz_shell_fraction * subsurface_mesh->cell_volume(cells_of_col[cc])) {
        std::cout << "On soil level " << cc << ", volume of shells = " << volume_shells[cc] << " but cell volume is " << subsurface_mesh->cell_volume(cells_of_col[cc]) << std::endl;
        std::cout << "  aroot_len = ";
        for (const auto& hp : hps) std::cout << " " << hp.aroot_length;
        std::cout << std::endl << "  rhiez_radius = " << rheiz_outer_radius[cc] << std::endl;
        Errors::Message msg("Plant1DMesh: volume issues -- rheizosphere shell volume is bigger than half the cell volume.  We hope to keep it << the cell volume.");
        Exceptions::amanzi_throw(msg);
      }
    }
    std::cout << "Successful with rheiz radii = ";
    for (int cc=0; total_aroot_length[cc] > 0. && cc!=cells_of_col.size(); ++cc) {
      std::cout << " " << rheiz_outer_radius[cc];
    }
    std::cout << std::endl;
  }

  // create the logical mesh
  auto logical_mesh = fac.Create();
  AmanziMesh::MeshLogicalAudit audit(logical_mesh, std::cout);
  ASSERT(!audit.Verify());
  
  // now create the embedded mesh
  auto embedded = Teuchos::rcp(new AmanziMesh::MeshEmbeddedLogical(subsurface_mesh->get_comm(),
          subsurface_mesh, logical_mesh, soil_rheiz_conn_cells, soil_rheiz_conn_lengths,
          soil_rheiz_conn_area_normals));
  embedded->set_geometric_model(gm);
  return embedded;
}
          
          
void
plantMesh(AmanziMesh::MeshLogicalFactory& fac,
          const HydraulicMeshParams& p,
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
          std::vector<double>& belowground_soil_volume)
{
  // input DBC
#ifdef ENABLE_DBC
  ASSERT(p.Valid());
  double rf_total = 0.;
  for (auto r : rooting_fraction) {
    ASSERT(r >= 0.);
    rf_total += r;
  }
  ASSERT(std::abs(1-rf_total) < 1.e-8);
  for (auto dz : soil_col_dz) {
    ASSERT(dz > 0);
  }
  ASSERT(n_stem > 0);
  ASSERT(n_rheiz_shell >= 0);
#endif

  
  // create space
  std::vector<AmanziGeometry::Point> centroids;
  std::vector<double> cell_volumes;
  std::vector<double> cell_lengths;
  std::vector<double> face_areas;
  AmanziGeometry::Point orientation(0.,0.,1.);
  AmanziGeometry::Point stem_centroid(soil_col_centroid + stem_displacement);

  // output
  belowground_soil_volume.resize(soil_col_dz.size(), 0.);
  
  // leaf
  centroids.resize(1,stem_centroid);
  centroids[0][2] = stem_centroid[2] + p.height;
  cell_volumes.resize(1, p.leaf_volume);
  cell_lengths.resize(1, p.leaf_height);
  face_areas.resize(1, p.sapwood_area);

  std::vector<AmanziGeometry::Entity_ID> leaf_cells, leaf_faces;
  fac.AddSegment(&centroids, &cell_volumes, cell_lengths, face_areas, orientation,
                 AmanziMesh::MeshLogicalFactory::LogicalTip_t::BOUNDARY,
                 AmanziMesh::MeshLogicalFactory::LogicalTip_t::JUNCTION,
                 p.name+"_leaf",
                 &leaf_cells, &leaf_faces);

  // reserve for the leaf-stem face
  int f = fac.ReserveFace();

  // add the stem
  double stem_cell_length = p.stem_length / n_stem;
  centroids = std::vector<AmanziGeometry::Point>(n_stem, stem_centroid);
  for (int i=0; i!=n_stem; ++i) {
    centroids[i][2] = stem_centroid[2] + p.stem_length - stem_cell_length * (i+0.5);
  }
  cell_lengths = std::vector<double>(n_stem, stem_cell_length);
  face_areas = std::vector<double>(n_stem-1, p.sapwood_area);
  cell_volumes = std::vector<double>(n_stem, stem_cell_length * p.sapwood_area); // NO TAPER

  std::vector<AmanziGeometry::Entity_ID> stem_cells, stem_faces;
  fac.AddSegment(&centroids, &cell_volumes, cell_lengths, face_areas, orientation,
                 AmanziMesh::MeshLogicalFactory::LogicalTip_t::BRANCH,
                 AmanziMesh::MeshLogicalFactory::LogicalTip_t::JUNCTION,
                 p.name+"_stem",
                 &stem_cells, &stem_faces);

  // add the leaf-stem face
  std::vector<int> leaf_stem_face_cells = { leaf_cells.back(), stem_cells.front() };
  std::vector<double> leaf_stem_face_lengths = { p.leaf_height / 2., stem_cell_length / 2. };
  fac.AddFace(f, leaf_stem_face_cells, orientation, leaf_stem_face_lengths, p.sapwood_area);

  
  // reserve for the stem-troot face
  f = fac.ReserveFace();

  // add the TRoot
  int n_lev_roots = 0.;
  ASSERT(rooting_fraction.size() == soil_col_dz.size());
  while (rooting_fraction[n_lev_roots] > 0.) {
    n_lev_roots++;
  }
  
  centroids.resize(n_lev_roots);
  auto tracking_centroid(stem_centroid);
  cell_lengths.resize(n_lev_roots);
  for (int i=0; i!=n_lev_roots; ++i) {
    tracking_centroid[2] -= soil_col_dz[i]/2.;
    centroids[i] = tracking_centroid;
    tracking_centroid[2] -= soil_col_dz[i]/2.;
    cell_lengths[i] = soil_col_dz[i];
  }
  face_areas.resize(n_lev_roots - 1);
  cell_volumes.resize(n_lev_roots);

  double sapwood_area_troot = p.sapwood_area;
  cell_volumes[0] = p.troot_volume * rooting_fraction[0];
  for (int i=1; i!=n_lev_roots; ++i) {
    sapwood_area_troot -= p.sapwood_area * rooting_fraction[i-1];
    face_areas[i-1] = sapwood_area_troot;
    cell_volumes[i] = p.troot_volume * rooting_fraction[i];
  }
  std::vector<AmanziGeometry::Entity_ID> troot_cells, troot_faces;
  fac.AddSegment(&centroids, &cell_volumes, cell_lengths, face_areas, orientation,
                 AmanziMesh::MeshLogicalFactory::LogicalTip_t::BRANCH,
                 AmanziMesh::MeshLogicalFactory::LogicalTip_t::JUNCTION,
                 p.name+"_troot",
                 &troot_cells, &troot_faces);

  //  for (int i=0; i!=n_lev_roots; ++i) {
  //    belowground_soil_volume[i] += cell_volumes[i];
  //  }

  // add the troot-stem face
  std::vector<int> stem_troot_face_cells = { stem_cells.back(), troot_cells.front() };
  std::vector<double> stem_troot_face_lengths = { stem_cell_length / 2., soil_col_dz[0] / 2. };
  fac.AddFace(f, stem_troot_face_cells, orientation, stem_troot_face_lengths, p.sapwood_area);

  // turn the orientation
  ASSERT(stem_displacement[2] == 0.);
  orientation = stem_displacement / AmanziGeometry::norm(stem_displacement);
  
  for (int slev=0; slev!=n_lev_roots; ++slev) {
    // reserve the aroot-troot
    f = fac.ReserveFace();

    AmanziGeometry::Point troot_centroid = centroids[slev];

    // add the aroot cell
    std::vector<AmanziGeometry::Point> l_centroids(1);
    l_centroids[0] = troot_centroid - orientation * p.troot_length * rooting_fraction[slev];
    cell_lengths = std::vector<double>(1, 0.); // does not matter since we have no faces
    face_areas.resize(0);
    cell_volumes = std::vector<double>(1, p.aroot_volume * rooting_fraction[slev]);

    std::vector<AmanziGeometry::Entity_ID> aroot_cells, aroot_faces;
    fac.AddSegment(&l_centroids, &cell_volumes, cell_lengths, face_areas, orientation,
                   AmanziMesh::MeshLogicalFactory::LogicalTip_t::BRANCH,
                   AmanziMesh::MeshLogicalFactory::LogicalTip_t::JUNCTION,
                   p.name+"_aroot_"+std::to_string(slev),
                   &aroot_cells, &aroot_faces);
    //    belowground_soil_volume[slev] += cell_volumes[0];

    // add the troot aroot face
    double aroot_troot_area = p.sapwood_area * rooting_fraction[slev];
    std::vector<int> troot_aroot_face_cells = { troot_cells[slev], aroot_cells[0] };
    std::vector<double> troot_aroot_face_lengths = { p.troot_length * rooting_fraction[slev]/2.,
                                                     p.troot_length * rooting_fraction[slev]/2. };
      // chosen to be equal to ensure equal conductance in troot and aroot,
      // assuming equal conductivity.  This is not physical length here!
    fac.AddFace(f, troot_aroot_face_cells, orientation, troot_aroot_face_lengths, aroot_troot_area);

    if (n_rheiz_shell > 0) {
      // reserve the aroot rheizosphere shell face
      f = fac.ReserveFace();
    
      // add the rheizosphere shells
      // -- set the rheizosphere shell inner/outer radii
      std::vector<double> rheiz_radii(n_rheiz_shell+1);
      rheiz_radii[0] = p.aroot_radius;
      rheiz_radii[n_rheiz_shell] = rheiz_outer_radius[slev];
      for (int i=1; i!=n_rheiz_shell; ++i) {
        rheiz_radii[i] = p.aroot_radius * std::pow(rheiz_outer_radius[slev] / p.aroot_radius, ((double) i)/n_rheiz_shell);
      }

      // -- set the centroids, volumes, areas
      l_centroids = std::vector<AmanziGeometry::Point>(n_rheiz_shell, stem_centroid);
      face_areas.resize(n_rheiz_shell - 1);
      cell_volumes.resize(n_rheiz_shell);
      cell_lengths.resize(n_rheiz_shell);

      for (int i=0; i!=n_rheiz_shell; ++i) {
        l_centroids[i] = troot_centroid - (p.troot_length * rooting_fraction[slev] + (rheiz_radii[i] + rheiz_radii[i+1])/2.) * orientation;
        cell_volumes[i] = rooting_fraction[slev] * p.aroot_length * PI * (std::pow(rheiz_radii[i+1],2) - std::pow(rheiz_radii[i],2));
        cell_lengths[i] = rheiz_radii[i+1] - rheiz_radii[i];
        
        if (i!=0) {
          face_areas[i-1] = rooting_fraction[slev]*p.aroot_length * 2 * PI * rheiz_radii[i];
        }
      }
      std::vector<AmanziGeometry::Entity_ID> rheiz_cells, rheiz_faces;
      fac.AddSegment(&l_centroids, &cell_volumes, cell_lengths, face_areas, orientation,
                     AmanziMesh::MeshLogicalFactory::LogicalTip_t::BRANCH,
                     AmanziMesh::MeshLogicalFactory::LogicalTip_t::JUNCTION,
                     p.name+"_rheizosphere_shell_"+std::to_string(slev),
                     &rheiz_cells, &rheiz_faces);
      for (int i=0; i!=n_rheiz_shell; ++i) belowground_soil_volume[slev] += cell_volumes[i];

      // add the shell-aroot connection
      std::vector<int> aroot_rheiz_face_cells = { aroot_cells.back(), rheiz_cells[0] };
      std::vector<double> aroot_rheiz_face_lengths = { rheiz_radii[0], (rheiz_radii[1] - rheiz_radii[0])/2. };
      double aroot_rheiz_area = rooting_fraction[slev] * p.aroot_length * 2 * PI * rheiz_radii[0];
      fac.AddFace(f, aroot_rheiz_face_cells, orientation, aroot_rheiz_face_lengths, aroot_rheiz_area);

      // save soil-rheiz connection info
      // Note background mesh ID must come second so embedding can know how to reset IDs
      soil_rheiz_conn_cells.emplace_back(std::vector<int>{rheiz_cells.back(), -1});
      soil_rheiz_conn_lengths.emplace_back(
          std::vector<double>{ (rheiz_radii[n_rheiz_shell] - rheiz_radii[n_rheiz_shell-1])/2., -1.0 } );
      soil_rheiz_conn_area_normals.emplace_back(rooting_fraction[slev] * p.aroot_length * 2 * PI * rheiz_radii[n_rheiz_shell]*orientation);

    } else {
      // save the soil-aroot connection info
      // Note background mesh ID must come second so embedding can know how to reset IDs
      soil_rheiz_conn_cells.emplace_back(std::vector<int>{aroot_cells.back(), -1});
      soil_rheiz_conn_lengths.emplace_back(
          std::vector<double>{ p.aroot_radius, -1.0 } );
      soil_rheiz_conn_area_normals.emplace_back(rooting_fraction[slev] * p.aroot_length * 2 * PI * p.aroot_radius * orientation);
    }
  }
}



} // namespace
} // namespace
