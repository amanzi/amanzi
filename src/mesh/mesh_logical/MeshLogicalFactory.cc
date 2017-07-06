/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "errors.hh"
#include "RegionEnumerated.hh"
#include "MeshLogicalFactory.hh"

#include "plant_1D_mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

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
  if (first_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    if (faces) faces->push_back(face_cell_list_.size());
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
  if (last_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    if (faces) faces->push_back(face_cell_list_.size());
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

  return;
}


// Add a segment of uniform size.
void
MeshLogicalFactory::AddSegment(const AmanziGeometry::Point& begin,
        const AmanziGeometry::Point& end,
        std::vector<double> lengths,
        std::vector<double> areas,
        std::vector<double> vols,
        MeshLogicalFactory::LogicalTip_t first_tip,
        MeshLogicalFactory::LogicalTip_t last_tip,
        const std::string& set_name,
        std::vector<Entity_ID>* cells,
        std::vector<Entity_ID>* faces) {
  if (faces) faces->clear();
  if (cells) cells->clear();
  int n_cells = vols.size();
  if (n_cells == 0) return;

  ASSERT(lengths.size() == vols.size());
  ASSERT(vols.size()-1 <= areas.size() <= vols.size() + 1);
  ASSERT(std::abs(std::accumulate(lengths.begin(), lengths.end(), 0.0) - 
                  AmanziGeometry::norm(end-begin)) < 1.e-10);

  int cell_first = cell_volume_.size();
  int cell_last = cell_first + n_cells;

  AmanziGeometry::Point normal(end-begin);
  normal *= 1. / AmanziGeometry::norm(normal);  

  // set cell ids
  std::vector<Entity_ID> new_cells(n_cells);
  for (int i=0; i!=n_cells; ++i) {
    new_cells[i] = cell_first + i;
  }
  if (cells) *cells = new_cells;
  
  // extend cell volumes
  cell_volume_.insert(cell_volume_.end(), vols.begin(), vols.end());

  // extend face_cell_list
  // - first face
  int f_index = 0;
  int c_index = 0;
  if (first_tip == MeshLogicalFactory::TIP_BOUNDARY) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(1,cell_first);
    std::vector<double> my_lengths(1, lengths[c_index]/2.0);
    AmanziGeometry::Point my_normal = -std::abs(areas[f_index]) * normal; // negate as it must be outward normal

    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(my_lengths);
    face_cell_normals_.push_back(my_normal);
    f_index++;
  }

  // - internal faces
  for (int i=cell_first; i!=cell_last-1; ++i) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(2);
    my_cells[0] = i; my_cells[1] = i+1;
    std::vector<double> my_lengths(2);
    my_lengths[0] = lengths[c_index]/2.; my_lengths[0] = lengths[c_index+1]/2.;
    AmanziGeometry::Point my_normal = -std::abs(areas[f_index]) * normal; // negate as it must be outward normal

    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(my_lengths);
    face_cell_normals_.push_back(my_normal);

    c_index++;
    f_index++;
  }

  // -last face
  if (last_tip == MeshLogicalFactory::TIP_BOUNDARY) {

    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(1, cell_last-1);
    std::vector<double> my_lengths(1);
    my_lengths[0] = lengths[c_index]/2.;
    AmanziGeometry::Point my_normal = -std::abs(areas[f_index]) * normal; // negate as it must be outward normal

    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(my_lengths);
    face_cell_normals_.push_back(my_normal);
  }

  // cell centroids
  AmanziGeometry::Point centroid(begin);
  centroid += normal * lengths[0]/2.;
  for (int i=0; i!=n_cells; ++i) {
    cell_centroids_.push_back(centroid);
    if (i != n_cells-1) {
      centroid += normal * (lengths[i]/2. + lengths[i+1]/2.);
    }
  }

  // create the region
  // - these are destroyed when the gm is destroyed
  Teuchos::RCP<AmanziGeometry::RegionEnumerated> enum_rgn =
      Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name, gm_->RegionSize(),
              "CELL", new_cells));
  gm_->AddRegion(enum_rgn);
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
  // check for compiled meshes
  if (plist.isParameter("plant mesh")) {
    std::string meshname = plist.get<std::string>("plant mesh");
    if (meshname == "1D test mesh") {
      return Testing::plantMesh(comm_, gm_, true);
    }
  }
  

  // Create each segment
  // - map to store metadata about previously inserted segments
  std::map<std::string, Entity_ID_List> seg_cells;
  std::map<std::string, Entity_ID_List> seg_faces;
  std::map<std::string, double> seg_ds;

  Teuchos::ParameterList& segments = plist.sublist("segments");
  for (auto segment_it=segments.begin(); segment_it!=segments.end(); ++segment_it) {
    std::string seg_name = segment_it->first;
    if (!segments.isSublist(seg_name)) {
      Errors::Message msg("MeshLogicalFactory: Malformed \"segments\" list, all items must be sublists.");
      Exceptions::amanzi_throw(msg);
    }

    // a single segment
    Teuchos::ParameterList& segment = segments.sublist(seg_name);
    int n_cells = segment.get<int>("number of cells");
    double area = segment.get<double>("cross sectional area [m^2]");

    // -- first tip
    Teuchos::Array<double> first_tip =
      segment.get<Teuchos::Array<double> >("first tip");
    if (first_tip.size() != 3) {
      Errors::Message msg("MeshLogicalFactory: logical meshes must be in 3D, so tips must be of length 3.");
      Exceptions::amanzi_throw(msg);
    }
    AmanziGeometry::Point first(first_tip[0], first_tip[1], first_tip[2]);

    std::string first_tip_type = segment.get<std::string>("first tip type");
    MeshLogicalFactory::LogicalTip_t first_type;
    if (first_tip_type == "boundary") {
      first_type = MeshLogicalFactory::TIP_BOUNDARY;
    } else if (first_tip_type == "junction") {
      first_type = MeshLogicalFactory::TIP_JUNCTION;
    } else if (first_tip_type == "branch") {
      first_type = MeshLogicalFactory::TIP_BRANCH;
    } else {
      Errors::Message msg("MeshLogicalFactory: invalid tip type specified.");
      Exceptions::amanzi_throw(msg);
      first_type = MeshLogicalFactory::TIP_NULL;
    }
    
    // -- last tip
    AmanziGeometry::Point last(3);
    if (segment.isParameter("last tip")) {
      Teuchos::Array<double> last_tip = segment.get<Teuchos::Array<double> >("last tip");

      if (last_tip.size() != 3) {
        Errors::Message msg("MeshLogicalFactory: logical meshes must be in 3D, so tips must be of length 3.");
        Exceptions::amanzi_throw(msg);
      }

      last[0] = last_tip[0];
      last[1] = last_tip[1];
      last[2] = last_tip[2];

    } else if (segment.isParameter("segment orientation")
               && segment.isParameter("segment length [m]")) {
      Teuchos::Array<double> seg_orientation =
        segment.get<Teuchos::Array<double> >("segment orientation");

      if (seg_orientation.size() != 3) {
        Errors::Message msg("MeshLogicalFactory: logical meshes must be in 3D, so tips/orientations must be of length 3.");
        Exceptions::amanzi_throw(msg);
      }

      AmanziGeometry::Point orientation(seg_orientation[0],
                                        seg_orientation[1],
                                        seg_orientation[2]);
      orientation /= AmanziGeometry::norm(orientation);

      double seg_length = segment.get<double>("segment length [m]");
      last = first + orientation * seg_length;

    } else {
      Errors::Message msg("MeshLogicalFactory: segment must either have \"last tip\" or \"segment orientation\" and \"segment length [m]\" to define the endpoint.");
      Exceptions::amanzi_throw(msg);
    }      


    std::string last_tip_type = segment.get<std::string>("last tip type");
    MeshLogicalFactory::LogicalTip_t last_type;
    if (last_tip_type == "boundary") {
      last_type = MeshLogicalFactory::TIP_BOUNDARY;
    } else if (last_tip_type == "junction") {
      last_type = MeshLogicalFactory::TIP_JUNCTION;
    } else if (last_tip_type == "branch") {
      last_type = MeshLogicalFactory::TIP_BRANCH;
    } else {
      Errors::Message msg("MeshLogicalFactory: invalid tip type specified.");
      Exceptions::amanzi_throw(msg);
      last_type = MeshLogicalFactory::TIP_NULL;
    }
    
    // add the segment
    Entity_ID_List cells,faces;
    double ds;
    AddSegment(n_cells, first, last, area, first_type, last_type, seg_name,
               &cells, &faces, &ds);
    seg_cells[seg_name] = cells;
    seg_faces[seg_name] = faces;
    seg_ds[seg_name] = ds;

    // if we have a branch, add the connection to the junction
    if (first_type == MeshLogicalFactory::TIP_BRANCH) {
      std::string junction_seg =
        segment.get<std::string>("first tip branch segment");

      std::map<std::string, Entity_ID_List>::const_iterator junction_cells =
        seg_cells.find(junction_seg);
      if (junction_cells == seg_cells.end()) {
        Errors::Message msg("MeshLogicalFactory: A branch segment must be listed AFTER the junction it branches from.");
        Exceptions::amanzi_throw(msg);
      }

      Entity_ID_List branch_cells(2);
      branch_cells[1] = cells.front();
      std::vector<double> lengths(2, ds/2);
      
      std::string junction_seg_tip =
        segment.get<std::string>("first tip branch segment tip");
      if (junction_seg_tip == "first") {
        branch_cells[0] = junction_cells->second.front();
      } else if (junction_seg_tip == "last") {
        branch_cells[0] = junction_cells->second.back();
      } else {
        Errors::Message msg("MeshLogicalFactory: \"branch segment tip\" must be \"first\" or \"last\"");
        Exceptions::amanzi_throw(msg);
      }

      AddConnection(branch_cells, last-first, lengths, area);
    }


    if (last_type == MeshLogicalFactory::TIP_BRANCH) {
      std::string junction_seg =
        segment.get<std::string>("last tip branch segment");

      std::map<std::string, Entity_ID_List>::const_iterator junction_cells =
        seg_cells.find(junction_seg);
      if (junction_cells == seg_cells.end()) {
        Errors::Message msg("MeshLogicalFactory: A branch segment must be listed AFTER the junction is branches from.");
        Exceptions::amanzi_throw(msg);
      }

      Entity_ID_List branch_cells(2);
      branch_cells[0] = cells.back();
      std::vector<double> lengths(2, ds/2.);
      
      std::string junction_seg_tip =
        segment.get<std::string>("last tip branch segment tip");
      if (junction_seg_tip == "last") {
        branch_cells[1] = junction_cells->second.front();
      } else if (junction_seg_tip == "last") {
        branch_cells[1] = junction_cells->second.back();
      } else {
        Errors::Message msg("MeshLogicalFactory: \"branch segment tip\" must be \"last\" or \"last\"");
        Exceptions::amanzi_throw(msg);
      }

      AddConnection(branch_cells, last-first, lengths, area);
    }
  }

  // add sets
  Teuchos::ParameterList& sets = plist.sublist("sets");
  
  for (auto set_it=sets.begin(); set_it!=sets.end(); ++set_it) {
    std::string set_name = set_it->first;
    Teuchos::Array<std::string> set_list =
      sets.get<Teuchos::Array<std::string> >(set_name);

    Entity_ID_List set_ents;
    for (Teuchos::Array<std::string>::const_iterator seg_name=set_list.begin();
         seg_name!=set_list.end(); ++seg_name) {

      std::map<std::string, Entity_ID_List>::const_iterator seg_cell =
        seg_cells.find(*seg_name);
      if (seg_cell == seg_cells.end()) {
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





  
} // namespace AmanziMesh
} // namespace Amanzi
