/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "errors.hh"
#include "Key.hh"
#include "RegionEnumerated.hh"
#include "MeshLogicalFactory.hh"

namespace Amanzi {
namespace AmanziMesh {


// Add a generic geometric segment.  Here geometric means that volumes are not
// provided, and that centroids are.  Volumes are calculated by areas times
// lengths.
void
MeshLogicalFactory::AddSegmentGeometric(
    std::vector<AmanziGeometry::Point> cell_centroids,
    std::vector<double> bisector_lengths,
    std::vector<double> areas,
    MeshLogicalFactory::LogicalTip_t first_tip_type,
    MeshLogicalFactory::LogicalTip_t last_tip_type,
    const AmanziGeometry::Point& orientation,
    const std::string& set_name,
    std::vector<Entity_ID>* cells,
    std::vector<Entity_ID>* faces)
{
  if (faces) faces->clear();
  if (cells) cells->clear();
  int n_cells = cell_centroids.size();
  if (n_cells == 0) return;

  if (n_cells == 1) {
    ASSERT(0); // FIX ME! --etc
    return;
  }

  // cell indices
  int cell_first = cell_volume_.size();
  int cell_last = cell_first + n_cells;

  // set cell ids
  std::vector<Entity_ID> new_cells(n_cells);
  for (int i=0; i!=n_cells; ++i) {
    new_cells[i] = cell_first + i;
  }
  if (cells) *cells = new_cells;

  // extend face_cell_list
  // - first face
  int f_index = 0;
  int l_index = 0;
  Entity_ID first_face = -1;
  if (first_tip_type == MeshLogicalFactory::TIP_BOUNDARY) {
    if (faces) faces->push_back(face_cell_list_.size());
    first_face = face_cell_list_.size();

    Entity_ID_List my_cells(1,cell_first);
    std::vector<double> my_lengths(1, bisector_lengths[l_index]);
    double my_area = areas[f_index];

    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(my_lengths);
    face_cell_normals_.emplace_back(-my_area * orientation);

    cell_volume_.push_back(my_area * my_lengths[0]);

    f_index++;
    l_index++;

  } else {
    cell_volume_.push_back(0.);
  }

  // - internal faces
  for (int i=cell_first; i!=cell_last-1; ++i) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(2);
    my_cells[0] = i; my_cells[1] = i+1;

    std::vector<double> my_lengths(2);
    my_lengths[0] = bisector_lengths[l_index]; my_lengths[1] = bisector_lengths[l_index+1];
    double my_area = areas[f_index];

    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(my_lengths);
    face_cell_normals_.emplace_back(my_area * orientation);

    cell_volume_.back() += my_area * my_lengths[0];
    cell_volume_.push_back(my_area * my_lengths[1]);

    l_index += 2;
    f_index++;
  }

  // -last face
  Entity_ID last_face = -1;
  if (last_tip_type == MeshLogicalFactory::TIP_BOUNDARY) {
    if (faces) faces->push_back(face_cell_list_.size());
    last_face = face_cell_list_.size();
    Entity_ID_List my_cells(1, cell_last-1);
    std::vector<double> my_lengths(1);
    my_lengths[0] = bisector_lengths[l_index];
    double my_area = areas[f_index];

    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(my_lengths);
    face_cell_normals_.emplace_back(my_area * orientation);

    cell_volume_.back() += my_area * my_lengths[0];
  }

  // Add centroids
  cell_centroids_.insert(cell_centroids_.end(), cell_centroids.begin(), cell_centroids.end());

  // create the region
  // - these are destroyed when the gm is destroyed
  auto enum_rgn = Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name, gm_->RegionSize(), "CELL", new_cells));
  gm_->AddRegion(enum_rgn);

  if (first_tip_type == MeshLogicalFactory::TIP_BOUNDARY) {
    Entity_ID_List boundary_face(1,first_face);
    auto boundary_rgn = Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name+"_first_tip", gm_->RegionSize(), "FACE", boundary_face));
    gm_->AddRegion(boundary_rgn);
  }
  
  if (last_tip_type == MeshLogicalFactory::TIP_BOUNDARY) {
    Entity_ID_List boundary_face(1,last_face);
    auto boundary_rgn = Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name+"_last_tip", gm_->RegionSize(), "FACE", boundary_face));
    gm_->AddRegion(boundary_rgn);
  }
}


void
MeshLogicalFactory::AddSegmentGeometric(Teuchos::ParameterList& plist) {
  // things we need to call the segment constructor
  auto seg_name = Keys::cleanPListName(plist.name());
  std::vector<AmanziGeometry::Point> cell_centroids;
  std::vector<double> bisector_lengths; // of size 2*ncells, the half-lengths
  std::vector<double> areas;   // ncells-1 <= size <= ncells+1, depending on tips.  the face areas
  MeshLogicalFactory::LogicalTip_t first_tip_type;
  MeshLogicalFactory::LogicalTip_t last_tip_type;

  // things we will use throughout to ensure consistency and calculate other things
  AmanziGeometry::Point first_tip;
  AmanziGeometry::Point last_tip;
  AmanziGeometry::Point orientation;
  double seg_length = -1;
  double seg_my_length = -1; // these differ when the tip is not a Boundary
  int n_cells = -1;

  //
  // Topology
  // =======================
  // -- number of cells
  n_cells = plist.get<int>("number of cells");

  // -- first tip type
  first_tip_type = GetTipType_(plist, "first tip type");

  // -- last tip type
  last_tip_type = GetTipType_(plist, "last tip type");
  
  // -- calculate number of faces, bisectors
  int n_faces = n_cells - 1; // internal faces
  int n_bisectors = 2*n_faces;
  if (first_tip_type == MeshLogicalFactory::TIP_BOUNDARY) {
    n_faces++;
    n_bisectors++;
  }
  if (last_tip_type == MeshLogicalFactory::TIP_BOUNDARY) {
    n_faces++;
    n_bisectors++;
  }

  //
  // Geometry
  // ====================
  // -- first tip location -- must come from PList
  first_tip = GetPoint_(plist, "first tip");
  seg_length = plist.get<double>("total length [m]"); // just provide this already!
  if (seg_length <= 0.) {
    Errors::Message msg;
    msg << "MeshLogicalFactory: Segment \"" << seg_name << "\": length provided must be positive.";
    Exceptions::amanzi_throw(msg);
  }
  
  // -- areas
  if (plist.isParameter("cross sectional area [m^2]")) {
    double area = plist.get<double>("cross sectional area [m^2]");
    areas.resize(n_faces, area);
  } else if (plist.isParameter("cross sectional areas [m^2]")) {
    areas = plist.get<Teuchos::Array<double> >("cross sectional areas [m^2]").toVector();
    if (areas.size() != n_faces) {
      Errors::Message msg;
      msg << "MeshLogicalFactory: Segment \"" << seg_name << "\": areas provided are of different length than segment number of faces requested.  Note this is often because of tip issues.";
      Exceptions::amanzi_throw(msg);
    }        
  } else {
    Errors::Message msg;
    msg << "MeshLogicalFactory: Segment \"" << seg_name << "\": areas must be provided through either \"cross sectional area [m^2]\" or \"cross sectional areas [m^2]\" parameters.";
    Exceptions::amanzi_throw(msg);
  }      

  // -- bisector lengths
  if (plist.isParameter("cell lengths [m]")) {
    auto cell_lengths = plist.get<Teuchos::Array<double> >("cell lengths [m]");
    if (cell_lengths.size() != n_cells) {
      Errors::Message msg;
      msg << "MeshLogicalFactory: Segment \"" << seg_name << "\": cell lengths provided are of different length than number of cells requested.";
      Exceptions::amanzi_throw(msg);
    }

    double seg_cell_lengths = std::accumulate(cell_lengths.begin(), cell_lengths.end(), 0.);
    if ((first_tip_type == MeshLogicalFactory::TIP_BOUNDARY) && (last_tip_type == MeshLogicalFactory::TIP_BOUNDARY)) {
      if (std::abs(seg_my_length - seg_length)/seg_length > 1.e-5) {
        Errors::Message msg;
        msg << "MeshLogicalFactory: Segment \"" << seg_name << "\": cell lengths provided do not sum to provided segment length.";
        Exceptions::amanzi_throw(msg);
      }
    }    

    bisector_lengths.resize(n_bisectors, -1.);
    int l_count = 0;
    if (first_tip_type == MeshLogicalFactory::TIP_BOUNDARY) {
      bisector_lengths[l_count] = cell_lengths[0]/2;
      ++l_count;
    }

    for (int c=0; c!=n_cells-1; ++c) {
      bisector_lengths[l_count] = cell_lengths[c]/2;
      bisector_lengths[l_count+1] = cell_lengths[c+1]/2;
      l_count += 2;
    }
    
    if (last_tip_type == MeshLogicalFactory::TIP_BOUNDARY) {
      bisector_lengths[l_count] = cell_lengths[n_cells-1]/2;
      ++l_count;
    }
    seg_my_length = std::accumulate(bisector_lengths.begin(), bisector_lengths.end(), 0.);

  } else if (plist.isParameter("bisector lengths [m]")) {
    bisector_lengths = plist.get<Teuchos::Array<double> >("bisector lengths [m]").toVector();
    if (bisector_lengths.size() != n_bisectors) {
      Errors::Message msg;
      msg << "MeshLogicalFactory: Segment \"" << seg_name << "\": \"bisector lengths [m]\" provided are of different length than number of half-cells requested.  This may be due to tip issues.";
      Exceptions::amanzi_throw(msg);
    }

    seg_my_length = std::accumulate(bisector_lengths.begin(), bisector_lengths.end(), 0.);

    if ((first_tip_type == MeshLogicalFactory::TIP_BOUNDARY) && (last_tip_type == MeshLogicalFactory::TIP_BOUNDARY)) {
      if (std::abs(seg_my_length - seg_length)/seg_length > 1.e-5) {
        Errors::Message msg;
        msg << "MeshLogicalFactory: Segment \"" << seg_name << "\": cell lengths provided do not sum to provided segment length.";
        Exceptions::amanzi_throw(msg);
      }
    } else if (seg_length <= seg_my_length) {
      Errors::Message msg;
      msg << "MeshLogicalFactory: Segment \"" << seg_name << "\": cell lengths provided are longer than the segment.";
      Exceptions::amanzi_throw(msg);
    }
    
  } else {
    // -- even division of provided length
    int n_total_bisectors = n_bisectors;
    if (first_tip_type == MeshLogicalFactory::TIP_BRANCH) {
      // add two for the eventual face connection branch to junction
      n_total_bisectors += 2;
    }
    if (last_tip_type == MeshLogicalFactory::TIP_BRANCH) {
      // add two for the eventual face connection branch to junction
      n_total_bisectors += 2;
    }

    bisector_lengths.resize(n_bisectors, seg_length / n_total_bisectors);
    seg_my_length = std::accumulate(bisector_lengths.begin(), bisector_lengths.end(), 0.);
  }



  // -- last tip and orientation
  //    Either provided directly, or given through first tip + (orientation * seg_length)
  if (plist.isParameter("last tip")) {
    last_tip = GetPoint_(plist, "last tip");

    orientation = last_tip - first_tip;
    if (AmanziGeometry::norm(orientation) != seg_length) {
      Errors::Message msg;
      msg << "MeshLogicalFactory: Segment \"" << seg_name
          << "\": inconsistency between provided lengths and first and last tip locations.";
      Exceptions::amanzi_throw(msg);
    }
    orientation /= seg_length;
        
  } else if (plist.isParameter("orientation")) {
    orientation = GetPoint_(plist, "orientation");
    orientation /= AmanziGeometry::norm(orientation);
        
    last_tip = first_tip + orientation * seg_length;

  } else {
    Errors::Message msg;
    msg << "MeshLogicalFactory: Segment \"" << seg_name
        << "\" must either have \"last tip\" or \"segment orientation\" and a length to define the endpoint.";
    Exceptions::amanzi_throw(msg);
  }      

  // -- cell centroids.  Need to shift for the case of branches and junctions
  ASSERT(seg_my_length >= 0.);
  ASSERT(seg_length >= seg_my_length-1.e-8);
  int n_buffers = 0;
  if (first_tip_type == MeshLogicalFactory::TIP_BRANCH) n_buffers++;
  if (last_tip_type == MeshLogicalFactory::TIP_BRANCH) n_buffers++;

  double s = n_buffers > 0 ? (seg_length - seg_my_length)/n_buffers : 0.;

  cell_centroids.resize(n_cells);
  int l_count = 0;
  if (first_tip_type == MeshLogicalFactory::TIP_BOUNDARY) {
    s += bisector_lengths[l_count];
    l_count++;
  }

  ASSERT(std::abs(AmanziGeometry::norm(orientation) - 1.) < 1.e-10);
  
  for (int c=0; c!=n_cells; ++c) {
    cell_centroids[c] = first_tip + s*orientation;
    s += bisector_lengths[l_count];
    s += bisector_lengths[l_count+1];
    l_count += 2;
  }

  //
  // Add the Segment
  // ====================
  std::vector<Entity_ID> cells;
  std::vector<Entity_ID> faces;  
  AddSegmentGeometric(cell_centroids, bisector_lengths, areas, first_tip_type, last_tip_type, orientation, seg_name,
                      &cells, &faces);
  seg_cells_[seg_name] = cells;
  seg_faces_[seg_name] = faces;

  //
  // Add Branches
  // ====================
  // if we have a branch, add the connection to the junction
  if (first_tip_type == MeshLogicalFactory::TIP_BRANCH) {
    // things we need to add the connetcion:
    Entity_ID_List branch_cells(2);
    std::vector<double> branch_bisector_lengths(2);
    double branch_area = 0.;
    // (also orientation)

    std::string junction_seg =
        plist.get<std::string>("first tip branch segment");

    auto junction_cells = seg_cells_.find(junction_seg);
    if (junction_cells == seg_cells_.end()) {
      Errors::Message msg("MeshLogicalFactory: A branch segment must be listed AFTER the junction it branches from.");
      Exceptions::amanzi_throw(msg);
    }

    branch_cells[1] = cells.front();
    branch_bisector_lengths[1] = bisector_lengths[0];
    branch_bisector_lengths[0] = AmanziGeometry::norm(cell_centroids[0] - first_tip) - branch_bisector_lengths[1];
    ASSERT(branch_bisector_lengths[0] > 0.);
  
    std::string junction_seg_tip =
        plist.get<std::string>("first tip branch segment tip");
    if (junction_seg_tip == "first") {
      branch_cells[0] = junction_cells->second.front();
    } else if (junction_seg_tip == "last") {
      branch_cells[0] = junction_cells->second.back();
    } else {
      Errors::Message msg("MeshLogicalFactory: \"branch segment tip\" must be \"first\" or \"last\"");
      Exceptions::amanzi_throw(msg);
    }

    if (plist.isParameter("cross sectional area [m^2]")) {
      branch_area = plist.get<double>("cross sectional area [m^2]");
    } else {
      branch_area = plist.get<double>("first tip branch segment area [m^2]");
    }
    
    AddConnection(branch_cells, orientation, branch_bisector_lengths, branch_area);
  }


  // if we have a branch, add the connection to the junction
  if (last_tip_type == MeshLogicalFactory::TIP_BRANCH) {
    // things we need to add the connetcion:
    Entity_ID_List branch_cells(2);
    std::vector<double> branch_bisector_lengths(2);
    double branch_area = 0.;
    // (also orientation)

    std::string junction_seg =
        plist.get<std::string>("last tip branch segment");

    auto junction_cells = seg_cells_.find(junction_seg);
    if (junction_cells == seg_cells_.end()) {
      Errors::Message msg("MeshLogicalFactory: A branch segment must be listed AFTER the junction it branches from.");
      Exceptions::amanzi_throw(msg);
    }

    branch_cells[0] = cells.back();
    branch_bisector_lengths[0] = bisector_lengths.back();
    branch_bisector_lengths[1] = AmanziGeometry::norm(last_tip - cell_centroids.back()) - branch_bisector_lengths[0];
    ASSERT(branch_bisector_lengths[1] > 0.);
  
    std::string junction_seg_tip =
        plist.get<std::string>("last tip branch segment tip");
    if (junction_seg_tip == "first") {
      branch_cells[1] = junction_cells->second.front();
    } else if (junction_seg_tip == "last") {
      branch_cells[1] = junction_cells->second.back();
    } else {
      Errors::Message msg("MeshLogicalFactory: \"branch segment tip\" must be \"first\" or \"last\"");
      Exceptions::amanzi_throw(msg);
    }

    if (plist.isParameter("cross sectional area [m^2]")) {
      branch_area = plist.get<double>("cross sectional area [m^2]");
    } else {
      branch_area = plist.get<double>("last tip branch segment area [m^2]");
    }
    
    AddConnection(branch_cells, orientation, branch_bisector_lengths, branch_area);
  }
}



// Add a vertical or horizontal segment of uniform size.
void
MeshLogicalFactory::AddSegment(int n_cells,
			       double length,
			       bool vertical,
			       double cross_section_area,
                               MeshLogicalFactory::LogicalTip_t first_tip,
                               MeshLogicalFactory::LogicalTip_t last_tip,
			       const std::string& set_name,
			       std::vector<Entity_ID>* cells,
			       std::vector<Entity_ID>* faces,
                               double* cell_length) {
  if (faces) faces->clear();
  if (cells) cells->clear();
  if (n_cells == 0) return;

  int cell_first = cell_volume_.size();
  int cell_last = cell_first + n_cells;

  // segment length
  int n_halflengths;
  int n_halflengths_before_my_first_cell;
  
  if (first_tip == MeshLogicalFactory::TIP_BOUNDARY
      || first_tip == MeshLogicalFactory::TIP_DEFERRED) {
    n_halflengths_before_my_first_cell = 1;

    if (last_tip == MeshLogicalFactory::TIP_BOUNDARY
        || last_tip == MeshLogicalFactory::TIP_DEFERRED) {
      // face-to-face:
      n_halflengths = n_cells*2;
    } else if (last_tip == MeshLogicalFactory::TIP_JUNCTION) {
      // face-to-cell
      n_halflengths = n_cells*2 - 1;
    } else {
      // face-to-cell, but not my cell
      n_halflengths = n_cells*2 + 1;
    }

  } else if (first_tip == MeshLogicalFactory::TIP_JUNCTION) {
    n_halflengths_before_my_first_cell = 0;

    if (last_tip == MeshLogicalFactory::TIP_BOUNDARY
        || last_tip == MeshLogicalFactory::TIP_DEFERRED) {
      // cell-to-face
      n_halflengths = n_cells*2 - 1;
    } else if (last_tip == MeshLogicalFactory::TIP_JUNCTION) {
      // cell-to-cell
      n_halflengths = n_cells*2 - 2;
    } else {
      // cell-to-cell, but not my cell
      n_halflengths = n_cells*2;
    }

  } else {
    n_halflengths_before_my_first_cell = 2;

    if (last_tip == MeshLogicalFactory::TIP_BOUNDARY
        || last_tip == MeshLogicalFactory::TIP_DEFERRED) {
      // cell-to-face
      n_halflengths = n_cells*2 + 1;
    } else if (last_tip == MeshLogicalFactory::TIP_JUNCTION) {
      // cell-to-cell
      n_halflengths = n_cells*2;
    } else {
      // cell-to-cell, but not my cell
      n_halflengths = n_cells*2 + 2;
    }
  }    

  double ds_halflength = length / n_halflengths;
  std::vector<double> ds_halflengths(2, ds_halflength);
  std::vector<double> ds_halflengths_boundary(1, ds_halflength);
  if (cell_length) *cell_length = 2*ds_halflength;


  AmanziGeometry::Point normal(0.,0.,0.);
  if (vertical) {
    normal[2] = cross_section_area;
  } else {
    normal[0] = cross_section_area;
  }

  // set cell ids
  std::vector<Entity_ID> new_cells(n_cells);
  for (int i=0; i!=n_cells; ++i) {
    new_cells[i] = cell_first + i;
  }
  if (cells) *cells = new_cells;
  
  // extend cell volumes
  std::vector<double> vols(n_cells, 2 * ds_halflength * cross_section_area);
  if (first_tip != MeshLogicalFactory::TIP_BOUNDARY
      && first_tip != MeshLogicalFactory::TIP_DEFERRED) {
    vols[0] -= ds_halflength * cross_section_area;
  }
  if (last_tip != MeshLogicalFactory::TIP_BOUNDARY
      && last_tip != MeshLogicalFactory::TIP_DEFERRED) {
    vols[vols.size()-1] -= ds_halflength * cross_section_area;
  }
  cell_volume_.insert(cell_volume_.end(), vols.begin(), vols.end());

  // extend face_cell_list
  // - first face
  Entity_ID first_face = -1;
  if (first_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    if (faces) faces->push_back(face_cell_list_.size());
    first_face = faces->front();
    Entity_ID_List my_cells(1,cell_first);
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths_boundary);
    face_cell_normals_.push_back(normal);
  }

  // - internal faces
  for (int i=cell_first; i!=cell_last-1; ++i) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(2);
    my_cells[0] = i; my_cells[1] = i+1;
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths);
    face_cell_normals_.push_back(normal);
  }

  // - last face
  Entity_ID last_face = -1;
  if (last_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    if (faces) faces->push_back(face_cell_list_.size());
    last_face = faces->back();
    Entity_ID_List my_cells(1,cell_last-1);
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths_boundary);
    face_cell_normals_.push_back(normal);
  }

  // create the region
  // - these are destroyed when the gm is destroyed
  Teuchos::RCP<AmanziGeometry::RegionEnumerated> enum_rgn =
      Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name, gm_->RegionSize(),
              "CELL", new_cells));
  gm_->AddRegion(enum_rgn);
  
  if (first_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    Entity_ID_List boundary_face(1,first_face);
    Teuchos::RCP<AmanziGeometry::RegionEnumerated> boundary_rgn =
        Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name+"_first_tip",
                gm_->RegionSize(), "FACE", boundary_face));
    gm_->AddRegion(boundary_rgn);
  }
  
  if (last_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    Entity_ID_List boundary_face(1,last_face);
    Teuchos::RCP<AmanziGeometry::RegionEnumerated> boundary_rgn =
        Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name+"_last_tip",
                gm_->RegionSize(), "FACE", boundary_face));
    gm_->AddRegion(boundary_rgn);
  }
  
  centroids_good_ = false;
  return;
}


// Add a segment of uniform size.
void
MeshLogicalFactory::AddSegment(int n_cells,
        const AmanziGeometry::Point& begin,
        const AmanziGeometry::Point& end,
        double cross_section_area,
        MeshLogicalFactory::LogicalTip_t first_tip,
        MeshLogicalFactory::LogicalTip_t last_tip,
        const std::string& set_name,
        std::vector<Entity_ID>* cells,
        std::vector<Entity_ID>* faces,
        double* cell_length) {
  if (faces) faces->clear();
  if (cells) cells->clear();
  if (n_cells == 0) return;

  int cell_first = cell_volume_.size();
  int cell_last = cell_first + n_cells;

  // segment length
  int n_halflengths;
  int n_halflengths_before_my_first_cell;
  
  if (first_tip == MeshLogicalFactory::TIP_BOUNDARY
      || first_tip == MeshLogicalFactory::TIP_DEFERRED) {
    n_halflengths_before_my_first_cell = 1;

    if (last_tip == MeshLogicalFactory::TIP_BOUNDARY
        || last_tip == MeshLogicalFactory::TIP_DEFERRED) {
      // face-to-face:
      n_halflengths = n_cells*2;
    } else if (last_tip == MeshLogicalFactory::TIP_JUNCTION) {
      // face-to-cell
      n_halflengths = n_cells*2 - 1;
    } else {
      // face-to-cell, but not my cell
      n_halflengths = n_cells*2 + 1;
    }

  } else if (first_tip == MeshLogicalFactory::TIP_JUNCTION) {
    n_halflengths_before_my_first_cell = 0;

    if (last_tip == MeshLogicalFactory::TIP_BOUNDARY
        || last_tip == MeshLogicalFactory::TIP_DEFERRED) {
      // cell-to-face
      n_halflengths = n_cells*2 - 1;
    } else if (last_tip == MeshLogicalFactory::TIP_JUNCTION) {
      // cell-to-cell
      n_halflengths = n_cells*2 - 2;
    } else {
      // cell-to-cell, but not my cell
      n_halflengths = n_cells*2;
    }

  } else {
    n_halflengths_before_my_first_cell = 2;

    if (last_tip == MeshLogicalFactory::TIP_BOUNDARY
      || last_tip == MeshLogicalFactory::TIP_DEFERRED) {
      // cell-to-face
      n_halflengths = n_cells*2 + 1;
    } else if (last_tip == MeshLogicalFactory::TIP_JUNCTION) {
      // cell-to-cell
      n_halflengths = n_cells*2;
    } else {
      // cell-to-cell, but not my cell
      n_halflengths = n_cells*2 + 2;
    }
  }    

  double ds_halflength = AmanziGeometry::norm(end-begin) / n_halflengths;
  std::vector<double> ds_halflengths(2, ds_halflength);
  std::vector<double> ds_halflengths_boundary(1, ds_halflength);
  if (cell_length) *cell_length = 2*ds_halflength;

  AmanziGeometry::Point normal(end-begin);
  normal *= cross_section_area / AmanziGeometry::norm(normal);  

  // set cell ids
  std::vector<Entity_ID> new_cells(n_cells);
  for (int i=0; i!=n_cells; ++i) {
    new_cells[i] = cell_first + i;
  }
  if (cells) *cells = new_cells;
  
  // extend cell volumes
  std::vector<double> vols(n_cells, 2 * ds_halflength * cross_section_area);
  if (first_tip != MeshLogicalFactory::TIP_BOUNDARY
      && first_tip != MeshLogicalFactory::TIP_DEFERRED) {
    vols[0] -= ds_halflength * cross_section_area;
  }
  if (last_tip != MeshLogicalFactory::TIP_BOUNDARY
      && last_tip != MeshLogicalFactory::TIP_DEFERRED) {
    vols[vols.size()-1] -= ds_halflength * cross_section_area;
  }
  cell_volume_.insert(cell_volume_.end(), vols.begin(), vols.end());

  // extend face_cell_list
  // - first face
  Entity_ID first_face = -1;
  if (first_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    if (faces) faces->push_back(face_cell_list_.size());
    first_face = faces->front();
    Entity_ID_List my_cells(1,cell_first);
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths_boundary);
    face_cell_normals_.push_back(-normal); // negate as it must be outward normal
  }

  // - internal faces
  for (int i=cell_first; i!=cell_last-1; ++i) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(2);
    my_cells[0] = i; my_cells[1] = i+1;
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths);
    face_cell_normals_.push_back(normal);
  }

  // -last face
  Entity_ID last_face = -1;
  if (last_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    if (faces) faces->push_back(face_cell_list_.size());
    last_face = faces->back();
    Entity_ID_List my_cells(1,cell_last-1);
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths_boundary);
    face_cell_normals_.push_back(normal);
  }

  // cell centroids
  normal /= cross_section_area;
  for (int i=0; i!=n_cells; ++i) {
    cell_centroids_.push_back(begin + (n_halflengths_before_my_first_cell + 2*i)*ds_halflength*normal);
  }

  // create the region
  // - these are destroyed when the gm is destroyed
  Teuchos::RCP<AmanziGeometry::RegionEnumerated> enum_rgn =
      Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name, gm_->RegionSize(),
              "CELL", new_cells));
  gm_->AddRegion(enum_rgn);

  if (first_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    Entity_ID_List boundary_face(1,first_face);
    Teuchos::RCP<AmanziGeometry::RegionEnumerated> boundary_rgn =
        Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name+"_first_tip",
                gm_->RegionSize(), "FACE", boundary_face));
    gm_->AddRegion(boundary_rgn);
  }
  
  if (last_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    Entity_ID_List boundary_face(1,last_face);
    Teuchos::RCP<AmanziGeometry::RegionEnumerated> boundary_rgn =
        Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name+"_last_tip",
                gm_->RegionSize(), "FACE", boundary_face));
    gm_->AddRegion(boundary_rgn);
  }
  
  return;
}


// Manually add a connection, returning the face id.
int
MeshLogicalFactory::AddConnection(const Entity_ID_List& cells,
				  const AmanziGeometry::Point& normal,
				  const std::vector<double>& lengths,
				  double area) {
  if (cells.size() != 2) {
    Errors::Message msg("MeshLogicalFactory: connection added is improperly formed -- all connections need two cells.");
    Exceptions::amanzi_throw(msg);
  }
  
  int f = face_cell_list_.size();
  face_cell_list_.push_back(cells);
  face_cell_normals_.push_back(area/AmanziGeometry::norm(normal)*normal);
  face_cell_lengths_.push_back(lengths);

  cell_volume_[cells[0]] += area * lengths[0];
  cell_volume_[cells[1]] += area * lengths[1];
  return f;
}


int
MeshLogicalFactory::AddSet(const std::string& set_name,
                           const std::string& ent,
                           const Entity_ID_List& ents) {

  // create the region
  // - these are destroyed when the gm is destroyed
  Teuchos::RCP<AmanziGeometry::RegionEnumerated> enum_rgn = 
      Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name, gm_->RegionSize(),
              ent, ents));
  gm_->AddRegion(enum_rgn);
  return gm_->RegionSize() - 1;
}


// One-stop shop -- create the whole thing from PList  
Teuchos::RCP<MeshLogical>
MeshLogicalFactory::Create(Teuchos::ParameterList& plist) {
  // Create each segment
  // - map to store metadata about previously inserted segments
  Teuchos::ParameterList& segments = plist.sublist("segments");
  for (auto& segment_it : segments) {
    // These are declared here to help document what actually needs to exist.
    // There are a lot of logical consistencies across these that need to be
    // tracked, and a lot of different ways of providing the same info.
    
    // things we need to call the constructor
    std::string seg_name = segment_it.first;
    AddSegmentGeometric(segments.sublist(seg_name));
  }    

  // add any additional sets
  Teuchos::ParameterList& sets = plist.sublist("sets");
  
  for (auto& set_it : sets) {
    std::string set_name = set_it.first;
    auto set_list = sets.get<Teuchos::Array<std::string> >(set_name);

    Entity_ID_List set_ents;
    for (auto& seg_name : set_list) {
      auto seg_cell = seg_cells_.find(seg_name);
      if (seg_cell == seg_cells_.end()) {
        Errors::Message msg("MeshLogicalFactory: set requests segment not created.");
        Exceptions::amanzi_throw(msg);
      }

      set_ents.insert(set_ents.end(), seg_cell->second.begin(), seg_cell->second.end());
    }
    AddSet(set_name, "cell", set_ents);
  }
  return Create();
}
  

// Create the mesh
Teuchos::RCP<MeshLogical>
MeshLogicalFactory::Create() {
  Teuchos::RCP<AmanziMesh::MeshLogical> mesh;
  if (centroids_good_) {
    ASSERT(cell_centroids_.size() == cell_volume_.size());
    mesh = Teuchos::rcp(new MeshLogical(comm_,
					cell_volume_,
					face_cell_list_,
					face_cell_lengths_,
					face_cell_normals_,
					&cell_centroids_));
    
  } else {
    mesh = Teuchos::rcp(new MeshLogical(comm_,
					cell_volume_,
					face_cell_list_,
					face_cell_lengths_,
					face_cell_normals_));
  }
  mesh->set_geometric_model(gm_);
  return mesh;
}


MeshLogicalFactory::LogicalTip_t
MeshLogicalFactory::GetTipType_(Teuchos::ParameterList& plist, const std::string& pname)
{
  std::string tip_type_s = plist.get<std::string>(pname);
  auto tip_type = MeshLogicalFactory::TIP_NULL;
  if (tip_type_s == "boundary") {
    tip_type = MeshLogicalFactory::TIP_BOUNDARY;
  } else if (tip_type_s == "junction") {
    tip_type = MeshLogicalFactory::TIP_JUNCTION;
  } else if (tip_type_s == "branch") {
    tip_type = MeshLogicalFactory::TIP_BRANCH;
  } else {
    Errors::Message msg("MeshLogicalFactory: invalid tip type specified.");
    Exceptions::amanzi_throw(msg);
  }
  return tip_type;
}


AmanziGeometry::Point
MeshLogicalFactory::GetPoint_(Teuchos::ParameterList& plist, const std::string& pname)
{
  auto point_a = plist.get<Teuchos::Array<double> >(pname);
  if (point_a.size() != 3) {
    Errors::Message msg("MeshLogicalFactory: logical meshes must be in 3D, so points must be of length 3.");
    Exceptions::amanzi_throw(msg);
  }
  return AmanziGeometry::Point{point_a[0], point_a[1], point_a[2]};
}






  
} // namespace AmanziMesh
} // namespace Amanzi
