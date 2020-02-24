/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! A factory for creating control-volume logical meshes from segments.

/*!

Builds up a logical mesh from a series of segments.  Note that logical meshes
include cell centroids, but only to make plotting easier -- no geometric
information from these are used.

Logical meshes need the following information:

- cell volume
- face area (and an orientation/normal for problems with gravity)
- length of the cell-to-face connector

This factory looks to build a logical mesh from a series of 1D segments which
may connect to other segments.  From this concept, the needed information can
be provided in a few ways, or inferred.  The primary distinction is between
cell volumes that are provided in the input and cell volumes that are inferred
by the equation:

V_c = sum_{each face in c} length of cell-to-face connector * face area.

* `"calculate cell volumes from lengths and areas`" ``[bool]`` **true** If true, calculate via the
    above formula.  Otherwise each segment must provide the `"cell volume`" or
    `"cell volumes`" parameter described below.

* `"segments`" ``[segment-spec]``  List of segment specs below:

Each segment consists of a collection of cells, and it is assumed that
face-to-cell connectors within the segment are half of the cell length.  For
each segment, cell volumes, face areas, an orientation, and the length of
cell-to-face connectors must be specified:

Cell volumes are either determined by the above formula if the `"cell volumes
provided`" option is specified, or through the following option:

* `"cell volumes [m^3]`" ``[Array(double)]`` list of volumes, [m^{2,3}], depending
   upon mesh dimension.

Cell lengths may be specified either as a list of cell lengths or as a segment
length and a number of cells (describing uniform cell length).

* `"segment length [m]`" ``[double]`` total length, [m]
* `"number of cells`" ``[int]`` number of cells

OR

* `"cell lengths [m]`" ``[Array(double)]`` List of cell lengths.

OR

* `"first tip`" ``[Array(double)]`` The segment start point.
* `"last tip`" ``[Array(double)]`` The last tip.

Note the order of checking is as above -- if both begin/end points and segment
lengths are provided, segment lengths will be used as specified, while the
begin/end points will only be used for visualization.

Face areas are specified through:

* `"cross sectional area [m^2]`" ``[double]`` cross sectional area of all faces in the segment.

OR

* `"cross sectional areas [m^2]`" ``[Array(double)]`` List of face areas, ordered from first tip to last tip.

Note that the number of faces, and hence the length of the list of face areas,
depends upon how exactly tips are to be handled.  See below.

Orientation of the faces is given by:

* `"orientation`" ``[Array(double)]`` Vector of length dimension describing
    the face normals.  This could be exteneded to be an Array of vectors, but
    currently that is not supported.  Note this is only used if gravity is
    involved.

OR

If segment begin/end points are provided, the orientation is given by
first tip - last tip.  Note that this is consistent with the view of a
logical mesh as a stream network, where first tip is the downstream outlet and
last tip is the upstream inlet, and orientation is aligned with the flow direction
(downstream).

Finally, often tips of segments may be boundaries (faces with boundary
conditions), no type (i.e. there is no face), junction (at least one, likely
more, faces will be connected to the end cell eventually), or branch (the tip
cell will add a face to a junction tip of another segment).  

* `"first tip type`" ``[string]`` Type of the first tip.  One of `"none`",
   `"boundary`", `"junction`", or `"branch`".
* `"last tip type`" ``[string]`` Type of the last tip, see above.

If either last (or respectively first) is a branch, then the following should
also be provided:

* `"last tip branch segment`" ``[string]`` Name of the segment to connect to.
* `"last tip branch segment tip`" ``[string]`` One of `"first`" or `"last`",
    the tip to connect this branch face to.

Note that segments with a branch tip MUST appear AFTER the junction tip to
which they will be connected.

This information then lets us know how many face areas to expect.  The face
area array must be at least as long as the number of cells - 1 (interior
faces).  On `"none`" and `"junction`" tips, face areas are not included, as
that face either doesn't exist (for `"none`") or will have its area specified
by any corresponding `"branch`" segments.  On `"boundary`" or `"branch`" tips,
face areas of the end face must be is included.

Finally, a note on centroids.  Centroids can be calculated by providing a
`"first tip`" and `"last tip`".  Note that there is no
expectation that these coordinates differ by the length of the segment.
Instead, cell centroids are located proportionally to the full segment legnth.
This allows meanding streams, for instance, to be visualized as straight lines
while still have the correct solution.

Alternatively, cell centroids can be manually specified:

* `"cell centroids`" ``[Array(double)]`` A list of dimension * n_cells length,
    where the first d values specify cell 0's centroid, the next d specify cell
    1's centroid, etc.

Finally, there are a few factory-global options:

Cell centroids can be inferred by specifying the FIRST segment's first tip,
and orientations and lengths for each segment that do not use first/last tip
coordinates.  Centroids can be spaced out, assuming the first tip of each
branch is a branch tip and the junction it connects to already has a last tip
coordinate.  In this case, specify:

* `"infer segment begin and end coordinates`" ``[bool]`` **false**.  When
    true, uses the above algorithm.

This option is provided in the `"logical from segments parameters`"
list.
    
As an example, consider the confluence of the Allegheny and the Monongahela,
which come together in Pittsburgh to form the Ohio.  This could be simply
modeled by three segments:

.. code-block:: xml

  <ParameterList name="mesh">
    <ParameterList name="my mesh">
      <Parameter name="mesh type" type="string" value="logical from segments" />
      <ParameterList name="logical from segments parameters">
        <Parameter name="cell volumes provided" type="bool" value="false" />
        <Parameter name="infer segment begin and end coordaintes" type="bool" value="true" />
        <ParameterList name="segments">
          <ParameterList name="ohio">
            <Parameter name="number of cells" type="int" value="10" />
            <Parameter name="segment length [m]" type="double" value="10000." />
            <Parameter name="cross sectional area [m^2]" type="double" value="300.0" />
            <Parameter name="first tip" type="Array(double)" value="{0.0,0.0}" />
            <Parameter name="first tip type" type="string" value="boundary" />
            <Parameter name="last tip type" type="string" value="junction" />
          </ParameterList>
          <ParameterList name="allegheny">
            <Parameter name="number of cells" type="int" value="10" />
            <Parameter name="segment length [m]" type="double" value="7000." />
            <Parameter name="cross sectional area [m^2]" type="double" value="120.0" />
            <Parameter name="orientation" type="Array(double)" value="{1,1}" />
            <Parameter name="first tip type" type="string" value="branch" />
            <Parameter name="first tip branch segment" type="string" value="ohio" />
            <Parameter name="first tip branch segment tip" type="string" value="last" />
            <Parameter name="last tip type" type="string" value="boundary" />
          </ParameterList>
          <ParameterList name="monongahela">
            <Parameter name="number of cells" type="int" value="10" />
            <Parameter name="segment length [m]" type="double" value="8000." />
            <Parameter name="cross sectional area [m^2]" type="double" value="180.0" />
            <Parameter name="orientation" type="Array(double)" value="{-1,1}" />
            <Parameter name="first tip type" type="string" value="branch" />
            <Parameter name="first tip branch segment" type="string" value="ohio" />
            <Parameter name="first tip branch segment tip" type="string" value="last" />
            <Parameter name="last tip type" type="string" value="boundary" />
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>
    
 */
#include <numeric>

#include "Key.hh"
#include "RegionEnumerated.hh"
#include "MeshLogicalFactory.hh"


namespace Amanzi {
namespace AmanziMesh {


// Create the mesh
Teuchos::RCP<MeshLogical>
MeshLogicalFactory::Create() const {
  Teuchos::RCP<AmanziMesh::MeshLogical> mesh;
  if (cell_centroids_.size() > 0) {
    AMANZI_ASSERT(cell_centroids_.size() == cell_volumes_.size());
    mesh = Teuchos::rcp(new MeshLogical(comm_,
            face_cell_list_,
            face_cell_dirs_,
            cell_volumes_,
            face_areas_,
            face_cell_bisectors_,
            &cell_centroids_));
    
  } else {
    mesh = Teuchos::rcp(new MeshLogical(comm_,
            face_cell_list_,
            face_cell_dirs_,
            cell_volumes_,
            face_areas_,
            face_cell_bisectors_,
            nullptr));
  }
  mesh->set_geometric_model(gm_);
  return mesh;
}


// One-stop shop -- create the whole thing from PList  
Teuchos::RCP<MeshLogical>
MeshLogicalFactory::Create(Teuchos::ParameterList& plist) {
  // set global options
  tracking_centroids_ = plist.get<bool>("infer cell centroids", false);
  if (plist.isParameter("calculate cell volumes from lengths and areas"))
    calculated_volume_ = plist.get<bool>("calculate cell volumes from lengths and areas");
  if (plist.isParameter("cell volumes provided")) {
    Errors::Message msg("MeshLogicalFactory: deprecated parameter \"cell volumes provided\", please use \"calculate cell volumes from lengths and areas\".");
    Exceptions::amanzi_throw(msg);
  }

  // Create each segment
  // - map to store metadata about previously inserted segments
  Teuchos::ParameterList& segments = plist.sublist("segments");
  for (auto& segment_it : segments) {
    // These are declared here to help document what actually needs to exist.
    // There are a lot of logical consistencies across these that need to be
    // tracked, and a lot of different ways of providing the same info.
    
    // things we need to call the constructor
    std::string seg_name = segment_it.first;
    AddSegment(segments.sublist(seg_name));
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


void
MeshLogicalFactory::AddSegment(
    int n_cells,
    const AmanziGeometry::Point& begin,
    const AmanziGeometry::Point& end,
    double face_area,
    MeshLogicalFactory::LogicalTip_t first_tip_type,
    MeshLogicalFactory::LogicalTip_t last_tip_type,
    std::string const& name,
    std::vector<Entity_ID> *const cells,
    std::vector<Entity_ID> *const faces) {
  AMANZI_ASSERT(calculated_volume_);
    
  // orientation
  auto orientation = begin - end;
  double seg_length = AmanziGeometry::norm(orientation);
  orientation /= seg_length;

  // lengths
  double cell_length = seg_length / n_cells;
  std::vector<double> lengths(n_cells, cell_length);

  // how many areas?
  int n_faces = n_cells - 1;
  if (first_tip_type == LogicalTip_t::BOUNDARY) n_faces++;
  if (last_tip_type == LogicalTip_t::BOUNDARY) n_faces++;
  std::vector<double> face_areas(n_faces, face_area);

  // centroids
  std::vector<AmanziGeometry::Point> centroids(n_cells);
  auto my_centroid = begin - orientation * cell_length/2.;
  centroids[0] = my_centroid;
  for (int c=0; c!=n_cells-1; ++c) {
    my_centroid -= orientation * cell_length;
    centroids[c+1] = my_centroid;
  }

  AddSegment(&centroids, nullptr, lengths, face_areas, orientation,
             first_tip_type, last_tip_type, name, cells, faces);
}



// Add a segment
//
// Centroids and cell volumes are optional arguments.
//
// Cells and faces are optional return values containing the list of
// entities in the new segment.
//
// Master add a segment -- all the others call this one!
void
MeshLogicalFactory::AddSegment(
      std::vector<AmanziGeometry::Point> const * const cell_centroids,
      std::vector<double> const *const cell_volumes,
      std::vector<double> const& cell_lengths,
      std::vector<double> const& face_areas,
      AmanziGeometry::Point const& orientation,
      MeshLogicalFactory::LogicalTip_t first_tip_type,
      MeshLogicalFactory::LogicalTip_t last_tip_type,
      std::string const& seg_name,
      std::vector<Entity_ID> *const cells,
      std::vector<Entity_ID> *const faces) {

  // manage output
  if (faces) faces->clear();
  if (cells) cells->clear();

  // number of new entities to be added
  int n_cells = cell_lengths.size();
  if (n_cells == 0) return;
  int n_faces = face_areas.size();

  // check for consistency in sizes
  if (cell_centroids) AMANZI_ASSERT(cell_centroids->size() == n_cells);
  if (cell_volumes) AMANZI_ASSERT(cell_volumes->size() == n_cells);

  // check for the number of faces relative to number of cells
  int n_faces_expected = n_cells - 1; // interior faces
  if (first_tip_type == LogicalTip_t::BOUNDARY) {
    // need the area, we add the face
    n_faces_expected++;
  }
  if (last_tip_type == LogicalTip_t::BOUNDARY) {
    // need the area, we add the face
    n_faces_expected++;
  }
  AMANZI_ASSERT(n_faces_expected == n_faces);

  // check for expected existence of cell volumes
  AMANZI_ASSERT((cell_volumes == nullptr) == calculated_volume_);
  
  // set the new cell lids
  int cell_first = cell_volumes_.size();
  std::vector<Entity_ID> new_cells(n_cells);
  for (int i=0; i!=n_cells; ++i) {
    new_cells[i] = cell_first + i;
  }
  if (cells) *cells = new_cells;

  // set face ids
  int face_first = face_cell_list_.size();
  std::vector<Entity_ID> new_faces(n_faces);
  for (int i=0; i!=n_faces; ++i) {
    new_faces[i] = face_first + i;
  }
  if (faces) *faces = new_faces;

  // insert or make space for cell volumes
  if (cell_volumes) {
    cell_volumes_.insert(cell_volumes_.end(), cell_volumes->begin(), cell_volumes->end());
  } else {
    cell_volumes_.resize(cell_volumes_.size() + n_cells, 0.);
  }

  // insert face areas
  face_areas_.insert(face_areas_.end(), face_areas.begin(), face_areas.end());

  // insert cell lengths
  cell_lengths_.insert(cell_lengths_.end(), cell_lengths.begin(), cell_lengths.end());

  // ensure orientation is unit normal
  AmanziGeometry::Point my_orientation = orientation / AmanziGeometry::norm(orientation);
  
  // add the first face
  int i_face = 0;
  if (first_tip_type == LogicalTip_t::BOUNDARY) {
    face_cell_list_.emplace_back(std::vector<int>{new_cells[0]});
    face_cell_bisectors_.emplace_back(std::vector<AmanziGeometry::Point>{ cell_lengths.front()/2. * orientation });
    face_cell_dirs_.emplace_back(std::vector<int>{1});

    if (calculated_volume_) cell_volumes_[new_cells[0]] += face_areas[i_face] * cell_lengths[0]/2;

    i_face++;
  }

  // add the interior faces
  for (int j=1; j!=n_cells; ++j) {
    face_cell_list_.emplace_back(std::vector<int>{new_cells[j-1], new_cells[j]});
    face_cell_bisectors_.emplace_back(std::vector<AmanziGeometry::Point>{ -cell_lengths[j-1]/2 * orientation,
                                                                           cell_lengths[j]/2 * orientation });
    face_cell_dirs_.emplace_back(std::vector<int>{-1,1});

    if (calculated_volume_) {
      cell_volumes_[new_cells[j-1]] += face_areas[i_face] * cell_lengths[j-1]/2;
      cell_volumes_[new_cells[j]] += face_areas[i_face] * cell_lengths[j]/2;
    }
    i_face++;
  }

  // potentially add the last face
  if (last_tip_type == LogicalTip_t::BOUNDARY) {
    face_cell_list_.emplace_back(std::vector<int>{new_cells.back()});
    face_cell_bisectors_.emplace_back(std::vector<AmanziGeometry::Point>{ -cell_lengths.back()/2 * orientation});
    face_cell_dirs_.emplace_back(std::vector<int>{-1});

    if (calculated_volume_) {
      cell_volumes_[new_cells.back()] += face_areas[i_face] * cell_lengths.back()/2;
    }
    i_face++;
  }

  // add in centroids
  if (cell_centroids) {
    cell_centroids_.insert(cell_centroids_.end(), cell_centroids->begin(), cell_centroids->end());
  }

  // add sets
  AddSet(seg_name, "CELL", new_cells);
  if (first_tip_type == LogicalTip_t::BOUNDARY) AddSet(seg_name+"_first_tip", "FACE", Entity_ID_List(1,new_faces[0]));
  if (last_tip_type == LogicalTip_t::BOUNDARY) AddSet(seg_name+"_last_tip", "FACE", Entity_ID_List(1,new_faces.back()));

}


// Add segment from sublist
void
MeshLogicalFactory::AddSegment(Teuchos::ParameterList& plist) {

  // need the following info to call AddSegment
  std::vector<AmanziGeometry::Point> cell_centroids;
  std::vector<double> cell_volumes;
  std::vector<double> cell_lengths;
  std::vector<double> face_areas;
  AmanziGeometry::Point orientation;

  auto seg_name = Keys::cleanPListName(plist.name());

  // number and size of cells, total segment length
  double seg_length = -1.0;
  int n_cells = -1;
  
  if (plist.isParameter("cell lengths [m]")) {
    cell_lengths = plist.get<Teuchos::Array<double>>("cell lengths [m]").toVector();
    seg_length = std::accumulate(cell_lengths.begin(), cell_lengths.end(), 0.);
    n_cells = cell_lengths.size();

  } else if (plist.isParameter("segment length [m]")
             && plist.isParameter("number of cells")) {
    seg_length = plist.get<double>("segment length [m]");
    n_cells = plist.get<int>("number of cells");
    cell_lengths.resize(n_cells, seg_length / n_cells);

  } else if (plist.isParameter("first tip")
             && plist.isParameter("last tip")
             && plist.isParameter("number of cells")) {
    auto begin = GetPoint_(plist, "first tip");
    auto end = GetPoint_(plist, "last tip");
    seg_length = AmanziGeometry::norm(end - begin);
    n_cells = plist.get<int>("number of cells");
    cell_lengths.resize(n_cells, seg_length / n_cells);
  } else {
    Errors::Message msg;
    msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): unable to get number of cells and segment length.  See documentation.";
    Exceptions::amanzi_throw(msg);
    
  }

  // potentially get cell volumes
  if (!calculated_volume_) {
    if (plist.isParameter("cell volume [m^3]")) {
      double cv = plist.get<double>("cell volume [m^3]");
      cell_volumes.resize(n_cells, cv);
    } else if (plist.isParameter("cell volumes [m^3]")) {
      cell_volumes = plist.get<Teuchos::Array<double>>("cell volumes [m^3]").toVector();
      if (cell_volumes.size() != n_cells) {
        Errors::Message msg;
        msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): number of cells inconsistent with number of cell volumes provided.";
        Exceptions::amanzi_throw(msg);
      }
    } else {
      Errors::Message msg;
      msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): unable to get cell volumes for mesh with non-calculated volumes.  "
          << "See documentation";
      Exceptions::amanzi_throw(msg);
    }
  }

  // Determine tip types, and through this, the expected number of faces
  int n_faces_total = n_cells + 1;
  int n_faces_mine = n_cells - 1;
  int n_faces_to_addsegment = n_cells - 1;
  auto first_tip_type = GetTipType_(plist, "first tip type");
  if (first_tip_type == LogicalTip_t::BOUNDARY) {
    n_faces_mine++;
    n_faces_to_addsegment++;
  } else if (first_tip_type == LogicalTip_t::BRANCH) {
    n_faces_mine++;
  }

  auto last_tip_type = GetTipType_(plist, "last tip type");
  if (last_tip_type == LogicalTip_t::BOUNDARY) {
    n_faces_mine++;
    n_faces_to_addsegment++;
  } else if (last_tip_type == LogicalTip_t::BRANCH) {
    n_faces_mine++;
  }

  // Get face areas
  // -- mine include BRANCH tip face areas
  std::vector<double> face_areas_mine;
  if (plist.isParameter("cross sectional area [m^2]")) {
    double face_area = plist.get<double>("cross sectional area [m^2]");
    face_areas_mine.resize(n_faces_mine, face_area);
  } else if (plist.isParameter("cross sectional areas [m^2]")) {
    face_areas_mine = plist.get<Teuchos::Array<double>>("cross sectional areas [m^2]").toVector();
    if (face_areas_mine.size() != n_faces_mine) {
      Errors::Message msg;
      msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): Incorrect number of cross sectional areas provided.";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg;
    msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): cross sectional areas not specified.";
    Exceptions::amanzi_throw(msg);
  }

  // -- face_areas do not include BRANCH tip face areas
  auto start_it = face_areas_mine.begin();
  if (first_tip_type == LogicalTip_t::BRANCH) start_it++;
  auto end_it = face_areas_mine.end();
  if (last_tip_type == LogicalTip_t::BRANCH) end_it--;
  face_areas.insert(face_areas.end(), start_it, end_it);

  // get the orientation
  if (plist.isParameter("orientation")) {
    orientation = GetPoint_(plist, "orientation");
    orientation /= AmanziGeometry::norm(orientation);
  } else if (plist.isParameter("first tip")
             && plist.isParameter("last tip")) {
    auto begin = GetPoint_(plist, "first tip");
    auto end = GetPoint_(plist, "last tip");
    orientation = begin - end;
    orientation /= AmanziGeometry::norm(orientation);
  } else {
    Errors::Message msg;
    msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): Orientation not specified.";
    Exceptions::amanzi_throw(msg);
  }

  // get centroids
  if (plist.isParameter("cell centroids")) {
    auto centroids_raw = plist.get<Teuchos::Array<double>>("cell centroids").toVector();
    if (centroids_raw.size() != n_cells * 3) {
      Errors::Message msg;
      msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): Specified cell centroid vector must be of size 3 * n_cells.";
      Exceptions::amanzi_throw(msg);
    }
    for (int i=0; i!=n_cells; ++i) {
      cell_centroids.emplace_back(AmanziGeometry::Point(centroids_raw[3*i], centroids_raw[3*i+1], centroids_raw[3*i+2]));
    }
    

  } else if (plist.isParameter("first tip")
             && plist.isParameter("last tip")) {
    auto begin = GetPoint_(plist, "first tip");
    auto end = GetPoint_(plist, "last tip");
    double provided_len = AmanziGeometry::norm(end - begin);
    auto ds = (end - begin) / seg_length;

    cell_centroids.resize(n_cells);
    auto my_centroid = begin + ds * cell_lengths[0]/2;
    cell_centroids[0] = my_centroid;
    for (int c=0; c!=n_cells-1; ++c) {
      my_centroid += ds * (cell_lengths[c]/2 + cell_lengths[c+1]/2);
      cell_centroids[c+1] = my_centroid;
    }
  } else if (tracking_centroids_) {
    AmanziGeometry::Point begin;
    if (first_tip_type == LogicalTip_t::BOUNDARY) {
      begin = GetPoint_(plist, "first tip");
    } else {
      if (first_tip_type != LogicalTip_t::BRANCH) {
        Errors::Message msg;
        msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): In tracking centroids mode, all \"first tip type\" values must be BOUNDARY or BRANCH";
        Exceptions::amanzi_throw(msg);
      }
      std::string branch_from = plist.get<std::string>("first tip branch segment");
      begin = tracking_end_points_[branch_from];
    }

    auto end = begin - orientation * seg_length;
    auto ds = -orientation;

    cell_centroids.resize(n_cells);
    auto my_centroid = begin + ds * cell_lengths[0]/2;
    cell_centroids[0] = my_centroid;
    for (int c=0; c!=n_cells-1; ++c) {
      my_centroid += ds * (cell_lengths[c]/2 + cell_lengths[c+1]/2);
      cell_centroids[c+1] = my_centroid;
    }
    tracking_end_points_[seg_name] = end;
  }

  // Now start doing the add
  
  // -- if a branch, reserve space for the first face.  This keeps faces in some reasonable
  //    order internally, but isn't really necessary.
  std::vector<int> new_cells, new_faces;
  if (first_tip_type == LogicalTip_t::BRANCH) {
    std::string branch_from = plist.get<std::string>("first tip branch segment");
    std::string branch_from_tip = plist.get<std::string>("first tip branch segment tip");

    auto branch_from_seg = seg_cells_.find(branch_from);
    if (branch_from_seg == seg_cells_.end()) {
      Errors::Message msg;
      msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): branches from segment \"" << branch_from << "\" but this segment has not yet been specified.  Junctions must come before Branches!";
      Exceptions::amanzi_throw(msg);
    }

    // note cell_volumes_.size() will be the LID of the first cell of this segment!
    std::vector<int> cells = {-1, (int) cell_volumes_.size() };
    std::vector<AmanziGeometry::Point> bisectors = { AmanziGeometry::Point(), orientation * cell_lengths[0]/2. };
    std::vector<int> dirs = { 0, 1 };
    
    if (branch_from_tip == "first") {
      cells[0] = seg_cells_[branch_from].front();
      bisectors[0] = seg_orientations_[branch_from] * cell_lengths_[cells[0]]/2.;
      dirs[0] = 1;
    } else if (branch_from_tip == "last") {
      cells[0] = seg_cells_[branch_from].back();
      bisectors[0] = -seg_orientations_[branch_from] * cell_lengths_[cells[0]]/2.;
      dirs[0] = -1;
    } else {
      Errors::Message msg;
      msg << "MeshLogicalFactory (segment \"" << seg_name << "\"): \"first tip branch segment tip\" must be \"first\" or \"last\"";
      Exceptions::amanzi_throw(msg);
    }

    // add the first face
    Entity_ID f = ReserveFace();

    // -- add the interior faces
    auto centroids_p = cell_centroids.size() > 0 ? &cell_centroids : nullptr;
    auto volumes_p = cell_volumes.size() > 0 ? &cell_volumes : nullptr;

    AddSegment(centroids_p, volumes_p, cell_lengths, face_areas, orientation,
               first_tip_type, last_tip_type, seg_name, &new_cells, &new_faces);

    AddFace(f, cells, bisectors, dirs, face_areas_mine[0]);

    
  } else {
    // -- Just add the interior faces
    auto centroids_p = cell_centroids.size() > 0 ? &cell_centroids : nullptr;
    auto volumes_p = cell_volumes.size() > 0 ? &cell_volumes : nullptr;

    AddSegment(centroids_p, volumes_p, cell_lengths, face_areas, orientation,
               first_tip_type, last_tip_type, seg_name, &new_cells, &new_faces);
  }  

  seg_cells_[seg_name] = new_cells;
  seg_faces_[seg_name] = new_faces;
  seg_orientations_[seg_name] = orientation;

}
  

// Reserve a slot for a face (likely to be added via AddFace!)
int
MeshLogicalFactory::ReserveFace() {
  int f = face_cell_list_.size();
  face_cell_list_.emplace_back(Entity_ID_List());
  face_cell_bisectors_.emplace_back(std::vector<AmanziGeometry::Point>());
  face_cell_dirs_.emplace_back(std::vector<int>());
  face_areas_.push_back(-1.0);
  return f;
}


// Manually add a connection, returning the face id.
int
MeshLogicalFactory::AddFace(int f,
                            const Entity_ID_List& cells,
                            const std::vector<AmanziGeometry::Point>& bisectors,
                            const std::vector<int>& dirs,
                            double area) {
  if (f < 0) f = ReserveFace();
  
  if (cells.size() != 2) {
    Errors::Message msg("MeshLogicalFactory: connection added is improperly formed -- all connections need two cells.");
    Exceptions::amanzi_throw(msg);
  }
  
  face_cell_list_[f] = cells;
  face_cell_bisectors_[f] = bisectors;
  face_cell_dirs_[f] = dirs;
  face_areas_[f] = area;

  if (calculated_volume_) {
    AMANZI_ASSERT(cells[0] < cell_volumes_.size());
    AMANZI_ASSERT(cells[1] < cell_volumes_.size());
    cell_volumes_[cells[0]] += area * AmanziGeometry::norm(bisectors[0]);
    cell_volumes_[cells[1]] += area * AmanziGeometry::norm(bisectors[1]);
  }
  return f;
}

int
MeshLogicalFactory::AddSet(const std::string& set_name,
                           const std::string& ent,
                           const Entity_ID_List& ents) {

  // create the region
  // - these are destroyed when the gm is destroyed
  Teuchos::RCP<AmanziGeometry::RegionEnumerated> enum_rgn = 
      Teuchos::rcp(new AmanziGeometry::RegionEnumerated(set_name, gm_->size(),
              ent, ents));
  gm_->AddRegion(enum_rgn);
  return gm_->size() - 1;
}
  

MeshLogicalFactory::LogicalTip_t
MeshLogicalFactory::GetTipType_(Teuchos::ParameterList& plist, const std::string& pname)
{
  std::string tip_type_s = plist.get<std::string>(pname);
  LogicalTip_t tip_type = LogicalTip_t::NONE;
  if (tip_type_s == "boundary") {
    tip_type = LogicalTip_t::BOUNDARY;
  } else if (tip_type_s == "junction") {
    tip_type = LogicalTip_t::JUNCTION;
  } else if (tip_type_s == "branch") {
    tip_type = LogicalTip_t::BRANCH;
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
